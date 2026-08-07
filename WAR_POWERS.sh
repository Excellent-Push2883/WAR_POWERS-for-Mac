#!/usr/bin/env bash
# =====================================================================
#  WAR_POWERS  -  Gene deletion / duplication screen  (v2.8)
#  GRCh38 / hg38.  Screening only, NOT clinical. Confirm hits in IGV.
#
#  Fixes vs the original:
#   * flank-contamination: uses the CLEANER (higher) flank OR a genome-wide
#     median baseline, so a deletion that runs into a flank is still caught.
#   * sex chromosomes: X and Y scored against their own baselines.
#   * mitochondria: MT-* genes scored vs whole-MT median (heteroplasmy note).
#   * contig naming: auto-detects bare "7" vs "chr7" BAMs.
#   * v2.6: SENSITIVITY-FIRST. Gene-depth-vs-baseline is the sole caller;
#     flanks are reported for context but never gate/veto. Threshold 0.75.
# =====================================================================
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COORDS="$SCRIPT_DIR/gencode_gene_coords.GRCh38.nochr.tsv"
GENES="$SCRIPT_DIR/genes_clean.txt"
GTF_GZ="$SCRIPT_DIR/gencode.v44.annotation.gtf.gz"
RAW="$SCRIPT_DIR/.wp_raw.tsv"; FLK="$SCRIPT_DIR/.wp_flanks.tsv"; BODY="$SCRIPT_DIR/.wp_body.tsv"
OUT_ALL="$SCRIPT_DIR/deletion_all_data.tsv"; OUT_HITS="$SCRIPT_DIR/deletion_hits.tsv"
FLANK=50000
# ---- tunable thresholds (edit these to taste) ----
DEL_HET_RATIO=0.70   # gene depth < this * baseline  -> POSSIBLE het/partial deletion
DUP_RATIO=1.35       # gene depth > this * baseline  -> POSSIBLE duplication
FULL_DEL_ABS=5       # gene depth < this many x       -> LIKELY full deletion
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
HDR=$'gene\tchrom\tstart\tend\tgene_depth\tleft_depth\tright_depth\tavg_depth\tgene_to_avg_depth\tcall'

echo ""
echo "=================================================="
echo "  Executing WAR_POWERS  (deletion/duplication screen, v2.8)"
echo "=================================================="
echo ""
read -rp "Press Enter to start, or Ctrl+C to cancel: " _ || true
echo ""

# ---- Step 1: prerequisites ----
echo "[1/6] Checking prerequisites..."
if ! command -v samtools >/dev/null 2>&1; then
  echo "  ERROR: samtools is not installed."
  case "${OSTYPE:-}" in darwin*) echo "  Install with:  brew install samtools";; *) echo "  Install with:  sudo apt install samtools";; esac
  exit 1
fi
echo "  samtools: $(samtools --version 2>&1 | head -1)"
shopt -s nullglob; BAMS=("$SCRIPT_DIR"/*.bam); shopt -u nullglob
if [ "${#BAMS[@]}" -eq 0 ]; then echo "  ERROR: no .bam file in this folder."; exit 1; fi
if [ "${#BAMS[@]}" -gt 1 ]; then echo "  ERROR: keep only ONE .bam file in this folder."; exit 1; fi
BAM="${BAMS[0]}"; echo "  BAM: $(basename "$BAM")"
if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then echo "  Indexing BAM (a few minutes)..."; samtools index "$BAM"; fi
MTNAME="$(samtools idxstats "$BAM" 2>/dev/null | awk '$1=="MT"||$1=="M"||$1=="chrM"||$1=="chrMT"{print $1; exit}')"
echo "  MT contig: ${MTNAME:-none}"
FIRSTSQ="$(samtools view -H "$BAM" 2>/dev/null | awk -F'\t' '/^@SQ/{for(i=1;i<=NF;i++) if(substr($i,1,3)=="SN:"){print substr($i,4); exit}}')"
case "$FIRSTSQ" in chr*) BAMCHR="chr";; *) BAMCHR="";; esac
echo "  BAM contig naming: first='${FIRSTSQ:-?}'  query-prefix='${BAMCHR}'"

# ---- Step 2: GENCODE annotation ----
echo ""; echo "[2/6] GENCODE v44 annotation..."
if [ ! -f "$COORDS" ]; then
  if [ ! -f "$GTF_GZ" ]; then
    echo "  Downloading GENCODE v44 (~50 MB)..."
    if command -v curl >/dev/null 2>&1; then curl -L --progress-bar -o "$GTF_GZ" "$GENCODE_URL"
    elif command -v wget >/dev/null 2>&1; then wget -q --show-progress -O "$GTF_GZ" "$GENCODE_URL"
    else echo "  ERROR: need curl or wget."; exit 1; fi
  fi
  echo "  Building gene coordinate table..."
  gzip -dc "$GTF_GZ" | awk 'BEGIN{OFS="\t"} $3=="gene"{for(i=1;i<=NF;i++) if($i=="gene_name"){g=$(i+1); gsub(/[";]/,"",g); c=$1; sub(/^chr/,"",c); if(c=="M")c="MT"; print g,c,$4,$5; break}}' > "$COORDS"
  echo "  Coordinate table: $(wc -l < "$COORDS" | tr -d ' ') genes"
else echo "  Coordinate table exists - skipping."; fi

# ---- Step 3: gene list ----
echo ""; echo "[3/6] Gene list..."
if [ ! -f "$GENES" ]; then echo "  ERROR: genes_clean.txt not found in this folder."; exit 1; fi
total="$(grep -c '[^[:space:]]' "$GENES")"; echo "  $total genes"

# ---- Step 4: measure depth (ONE samtools bedcov pass over all gene + flank intervals) ----
echo ""; echo "[4/6] Measuring depth (single-pass coverage over all genes + flanks)..."
: > "$RAW"; : > "$FLK"
REGIONS="$SCRIPT_DIR/.wp_regions.bed"; BEDCOV="$SCRIPT_DIR/.wp_bedcov.tsv"; IDX="$SCRIPT_DIR/.wp_idx.tsv"
samtools idxstats "$BAM" 2>/dev/null | awk '{print $1"\t"$2}' > "$IDX"
# build a single BED: gene body (G), left flank (L), right flank (R); clamped to contig bounds
awk -v FLANK="$FLANK" -v BAMCHR="$BAMCHR" -v MTNAME="$MTNAME" '
  FNR==1{fno++}
  fno==1{ len[$1]=$2+0; next }
  fno==2{ if(!($1 in coord)) coord[$1]=$2 SUBSEP $3 SUBSEP $4; next }
  fno==3{ g=$1; sub(/\r$/,"",g); if(g=="")next; if(!(g in coord))next;
    split(coord[g],a,SUBSEP); c=a[1]; s=a[2]+0; e=a[3]+0;
    if(c=="MT"){ bc=MTNAME; if(bc=="")next; L=len[bc]; es=s; ee=e; if(es<1)es=1; if(L>0&&ee>L)ee=L;
      if(ee>=es)printf "%s\t%d\t%d\t%s|G\n",bc,es-1,ee,g; next }
    bc=(BAMCHR=="")?c:(BAMCHR c); L=len[bc];
    gs=s; ge=e; if(gs<1)gs=1; if(L>0&&ge>L)ge=L; if(ge>=gs)printf "%s\t%d\t%d\t%s|G\n",bc,gs-1,ge,g;
    lsf=s-FLANK; lef=s-1; if(lsf<1)lsf=1; if(lef>=lsf)printf "%s\t%d\t%d\t%s|L\n",bc,lsf-1,lef,g;
    rsf=e+1; ref=e+FLANK; if(L>0&&ref>L)ref=L; if(ref>=rsf)printf "%s\t%d\t%d\t%s|R\n",bc,rsf-1,ref,g;
  }' "$IDX" "$COORDS" "$GENES" > "$REGIONS"
nreg="$(wc -l < "$REGIONS" | tr -d ' ')"
echo "  $nreg intervals -> one samtools bedcov pass (this is the slow part; a few minutes)..."
samtools bedcov "$REGIONS" "$BAM" 2>/dev/null > "$BEDCOV"
# assemble RAW (gene chrom start end gd ld rd) + FLK (chrom flank_mean) from the summed coverage
awk -v OFS='\t' -v RAW="$RAW" -v FLK="$FLK" '
  FNR==1{fno++}
  fno==1{ split($4,x,"|"); w=$3-$2; mm[x[1] SUBSEP x[2]]=(w>0)?($NF/w):"NA"; next }
  fno==2{ if(!($1 in coord)) coord[$1]=$2 SUBSEP $3 SUBSEP $4; next }
  fno==3{ g=$1; sub(/\r$/,"",g); if(g=="")next;
    if(!(g in coord)){ print g,"NA","NA","NA","NA","NA","NA" > RAW; next }
    split(coord[g],a,SUBSEP); c=a[1]; s=a[2]; e=a[3];
    gd=((g SUBSEP "G") in mm)?sprintf("%.4f",mm[g SUBSEP "G"]):"NA";
    ld=((g SUBSEP "L") in mm)?sprintf("%.4f",mm[g SUBSEP "L"]):"NA";
    rd=((g SUBSEP "R") in mm)?sprintf("%.4f",mm[g SUBSEP "R"]):"NA";
    print g,c,s,e,gd,ld,rd > RAW;
    if(c!="MT"){ if(ld!="NA")print c,ld > FLK; if(rd!="NA")print c,rd > FLK }
  }' "$BEDCOV" "$COORDS" "$GENES"
rm -f "$REGIONS" "$BEDCOV" "$IDX"

# ---- Step 5: baselines + calls ----
echo ""; echo "[5/6] Per-chromosome baselines + calls..."
median(){ sort -n | awk '{a[NR]=$1} END{if(NR==0)print "NA"; else if(NR%2)print a[(NR+1)/2]; else printf "%.4f",(a[NR/2]+a[NR/2+1])/2}'; }
AUTO="$(awk -F'\t' '$1 ~ /^[0-9]+$/{print $2}' "$FLK" | median)"
XMED="$(awk -F'\t' '$1=="X"{print $2}' "$FLK" | median)"
YMED="$(awk -F'\t' '$1=="Y"{print $2}' "$FLK" | median)"
MTMED="NA"; if [ -n "$MTNAME" ]; then MTMED="$(samtools depth -a -r "$MTNAME" "$BAM" 2>/dev/null | awk '{print $3}' | median)"; fi; [ -z "$MTMED" ] && MTMED="NA"
echo "  baselines:  autosomal=${AUTO}x  X=${XMED}x  Y=${YMED}x  MT=${MTMED}x"
if [ "$AUTO" = "NA" ]; then
  echo "  WARNING: autosomal baseline is NA -> depth was empty for every gene."
  echo "           Check the BAM contig naming above and that the BAM is complete/indexed."
fi

awk -v AUTO="$AUTO" -v XMED="$XMED" -v YMED="$YMED" -v MTMED="$MTMED" -v DELHET="$DEL_HET_RATIO" -v DUP="$DUP_RATIO" -v FULLABS="$FULL_DEL_ABS" 'BEGIN{OFS="\t"}
{
  gene=$1;chrom=$2;start=$3;end=$4;gd=$5;l=$6;r=$7;
  if(chrom=="NA"){print gene,"NA","NA","NA","NA","NA","NA","NA","NA","NOT_FOUND_IN_GENCODE"; next}
  if(chrom=="MT"){
    if(MTMED=="NA"||MTMED+0<=0||gd=="NA"){rg="NA";call="MT_NOT_ASSESSED"}
    else{rg=sprintf("%.4f",gd/MTMED); call=((gd+0)/(MTMED+0)<0.5)?"MT_POSSIBLE_HETEROPLASMIC_DELETION_REVIEW":"MT_DEPTH_NORMAL"}
    print gene,chrom,start,end,gd,"NA","NA",MTMED,rg,call; next
  }
  if(chrom ~ /^[0-9]+$/)G=AUTO; else if(chrom=="X")G=XMED; else if(chrom=="Y")G=YMED; else G="NA";
  # left/right flank depths are shown for context only; the call uses gene vs avg depth
  rg=(G!="NA"&&G+0>0&&gd!="NA")?sprintf("%.4f",gd/G):"NA";
  # sole caller = gene depth vs the (per-chromosome) genome baseline
  usable=(G!="NA"&&G+0>=10);
  if(gd=="NA") call="NO_CALL";
  else if(!usable) call="NO_CALL";
  else if(gd+0<FULLABS+0) call="LIKELY_FULL_DELETION_OR_MAJOR_COPY_LOSS";
  else if(rg!="NA"&&rg+0<DELHET+0) call="POSSIBLE_PARTIAL_OR_HET_DELETION";
  else if(rg!="NA"&&rg+0>DUP+0) call="POSSIBLE_DUPLICATION_OR_HIGH_COPY";
  else call="NORMAL";
  print gene,chrom,start,end,gd,l,r,G,rg,call;
}' "$RAW" > "$BODY"

printf '%s\n' "$HDR" > "$OUT_ALL"; cat "$BODY" >> "$OUT_ALL"
printf '%s\n' "$HDR" > "$OUT_HITS"; grep 'DELETION' "$BODY" >> "$OUT_HITS" || true
rm -f "$RAW" "$FLK" "$BODY"

# ---- Step 6: results ----
echo ""; echo "[6/6] Done."
printf "  All genes : %s  (%s rows)\n" "$OUT_ALL" "$(( $(wc -l < "$OUT_ALL") - 1 ))"
printf "  Hits only : %s  (%s hits)\n" "$OUT_HITS" "$(( $(wc -l < "$OUT_HITS") - 1 ))"
echo ""; echo "Flagged (nuclear deletions + MT review):"; echo ""
if command -v column >/dev/null 2>&1; then column -t -s "$(printf '\t')" "$OUT_HITS"; else cat "$OUT_HITS"; fi
echo ""
echo "Reminder: screening only, not clinical. Confirm any hit in IGV / with the"
echo "depth + loss-of-heterozygosity check before treating it as real."
