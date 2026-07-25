#!/usr/bin/env bash

# ─────────────────────────────────────────────────────────────────
# WAR_POWERS — Gene Deletion / Duplication Screen
# Reference: GRCh38 / hg38, no "chr" prefix (sequencing.com)
# ─────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COORDS="$SCRIPT_DIR/gencode_gene_coords.GRCh38.nochr.tsv"
GENES="$SCRIPT_DIR/genes_clean.txt"
GTF_GZ="$SCRIPT_DIR/gencode.v44.annotation.gtf.gz"
OUT_ALL="$SCRIPT_DIR/deletion_all_data.tsv"
OUT_HITS="$SCRIPT_DIR/deletion_hits.tsv"
FLANK=50000
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"

echo ""
echo "=================================================="
echo "  Executing WAR_POWERS  (gene deletion / duplication screen)"
echo "  Reference : GRCh38 / hg38 (no chr prefix)"
echo "  Output    : $SCRIPT_DIR"
echo "=================================================="
echo ""
echo "This script will:"
echo "  1. Verify samtools is installed"
echo "  2. Detect your BAM file automatically"
echo "  3. Download GENCODE v44 annotation if needed (~50 MB)"
echo "  4. Build a gene coordinate table"
echo "  5. Screen all genes in genes_clean.txt for deletions and duplications"
echo "  6. Write deletion_all_data.tsv and deletion_hits.tsv"
echo ""
echo "Runtime estimate: 15-60 minutes depending on your machine."
echo ""
read -rp "Press Enter to start, or Ctrl+C to cancel: "
echo ""

# ── Step 1: Prerequisites ──────────────────────────────────────
echo "[1/5] Checking prerequisites..."

if ! command -v samtools &>/dev/null; then
    echo ""
    echo "  ERROR: samtools is not installed."
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "  Install with:  brew install samtools"
    else
        echo "  Install with:  sudo apt install samtools"
    fi
    echo ""
    exit 1
fi
echo "  samtools: $(samtools --version 2>&1 | head -1)"

# Auto-detect BAM file
BAM_FILES=("$SCRIPT_DIR"/*.bam)
if [ ! -f "${BAM_FILES[0]}" ]; then
    echo ""
    echo "  ERROR: No .bam file found in $SCRIPT_DIR"
    echo "  Place your BAM file in the same folder as this script."
    echo ""
    exit 1
elif [ ${#BAM_FILES[@]} -gt 1 ]; then
    echo ""
    echo "  ERROR: Multiple .bam files found in $SCRIPT_DIR:"
    for f in "${BAM_FILES[@]}"; do
        echo "    $(basename "$f")"
    done
    echo "  Please keep only one .bam file in this folder."
    echo ""
    exit 1
fi
BAM="${BAM_FILES[0]}"
echo "  BAM found: $(basename "$BAM")"

if [ ! -f "${BAM}.bai" ]; then
    echo "  No index found — building BAM index (this may take a few minutes)..."
    samtools index "$BAM"
fi
echo "  BAM index OK"

# ── Step 2: GENCODE annotation ─────────────────────────────────
echo ""
echo "[2/5] GENCODE v44 annotation..."

if [ ! -f "$COORDS" ]; then
    if [ ! -f "$GTF_GZ" ]; then
        echo "  Downloading gencode.v44.annotation.gtf.gz..."
        if command -v curl &>/dev/null; then
            curl -L --progress-bar -o "$GTF_GZ" "$GENCODE_URL"
        elif command -v wget &>/dev/null; then
            wget -q --show-progress -O "$GTF_GZ" "$GENCODE_URL"
        else
            echo "  ERROR: Neither curl nor wget found. Install one and retry."
            exit 1
        fi
    else
        echo "  GTF already present, skipping download."
    fi

    echo "  Parsing GTF into gene coordinate table..."
    gzip -dc "$GTF_GZ" | awk 'BEGIN{OFS="\t"} $3=="gene" {
        for(i=1;i<=NF;i++) {
            if($i=="gene_name") {
                gene=$(i+1)
                gsub(/[";]/, "", gene)
                chrom=$1
                sub(/^chr/, "", chrom)
                if(chrom=="M") chrom="MT"
                print gene, chrom, $4, $5
                break
            }
        }
    }' > "$COORDS"
    echo "  Coordinate table written: $(wc -l < "$COORDS" | tr -d ' ') genes"
else
    echo "  Coordinate table already exists — skipping download/parse."
fi

# ── Step 3: Gene list ──────────────────────────────────────────
echo ""
echo "[3/5] Checking gene list..."

if [ ! -f "$GENES" ]; then
    echo "  ERROR: genes_clean.txt not found in $SCRIPT_DIR"
    echo "  Make sure genes_clean.txt is in the same folder as this script."
    exit 1
fi

GENE_COUNT=$(grep -c '[^[:space:]]' "$GENES")
echo "  $GENE_COUNT genes loaded from genes_clean.txt"

# ── Step 4: Deletion screen ────────────────────────────────────
echo ""
echo "[4/5] Running deletion screen ($GENE_COUNT genes × 3 depth calls each)..."
echo "      Progress will print as genes are processed."
echo ""

avgdepth() {
    samtools depth -a -r "$1" "$BAM" 2>/dev/null | awk '
    {sum+=$3; n++}
    END {
        if(n>0) printf "%.4f", sum/n
        else    printf "NA"
    }'
}

# Write headers
printf "gene\tchrom\tstart\tend\tgene_depth\tleft_depth\tright_depth\tflank_mean\tgene_to_flank_mean\tgene_to_left\tgene_to_right\tcall\n" \
    > "$OUT_ALL"
printf "gene\tchrom\tstart\tend\tgene_depth\tleft_depth\tright_depth\tflank_mean\tgene_to_flank_mean\tgene_to_left\tgene_to_right\tcall\n" \
    > "$OUT_HITS"

total=$(wc -l < "$GENES" | tr -d ' ')
count=0

while IFS= read -r gene || [ -n "$gene" ]; do
    # Strip carriage returns (in case of CRLF)
    gene="${gene%$'\r'}"
    [ -z "$gene" ] && continue

    count=$((count + 1))
    printf "[%d/%d] %s\n" "$count" "$total" "$gene"

    hit=$(awk -v g="$gene" '$1==g {print $0; exit}' "$COORDS")

    if [ -z "$hit" ]; then
        printf "%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNOT_FOUND_IN_GENCODE\n" \
            "$gene" >> "$OUT_ALL"
        continue
    fi

    chrom=$(printf "%s" "$hit" | awk '{print $2}')
    start=$(printf "%s" "$hit" | awk '{print $3}')
    end=$(printf "%s"   "$hit" | awk '{print $4}')

    left_start=$((start - FLANK))
    left_end=$((start - 1))
    right_start=$((end + 1))
    right_end=$((end + FLANK))
    [ "$left_start" -lt 1 ] && left_start=1

    gene_d=$(avgdepth  "${chrom}:${start}-${end}")
    left_d=$(avgdepth  "${chrom}:${left_start}-${left_end}")
    right_d=$(avgdepth "${chrom}:${right_start}-${right_end}")

    flank_mean=$(awk -v l="$left_d" -v r="$right_d" 'BEGIN{
        if(l=="NA" && r=="NA") print "NA"
        else if(l=="NA")       print r
        else if(r=="NA")       print l
        else printf "%.4f", (l+r)/2
    }')

    ratio_mean=$(awk -v g="$gene_d" -v f="$flank_mean" 'BEGIN{
        if(g=="NA"||f=="NA"||f==0) print "NA"
        else printf "%.4f", g/f
    }')

    ratio_left=$(awk -v g="$gene_d" -v l="$left_d" 'BEGIN{
        if(g=="NA"||l=="NA"||l==0) print "NA"
        else printf "%.4f", g/l
    }')

    ratio_right=$(awk -v g="$gene_d" -v r="$right_d" 'BEGIN{
        if(g=="NA"||r=="NA"||r==0) print "NA"
        else printf "%.4f", g/r
    }')

    call=$(awk -v gd="$gene_d" -v ld="$left_d" -v rd="$right_d" \
               -v fm="$flank_mean" -v rl="$ratio_left" -v rr="$ratio_right" 'BEGIN{
        if(gd=="NA"||fm=="NA") {
            print "NO_CALL"
        }
        # Gene near-zero, at least one flank has real coverage
        else if(gd < 5 && ((ld!="NA" && ld >= 10) || (rd!="NA" && rd >= 10))) {
            print "LIKELY_FULL_DELETION_OR_MAJOR_COPY_LOSS"
        }
        # Gene lower than BOTH flanks (both flanks must be valid and high)
        else if(ld!="NA" && rd!="NA" && ld >= 10 && rd >= 10 && rl < 0.75 && rr < 0.75) {
            print "POSSIBLE_PARTIAL_OR_HET_DELETION"
        }
        # Gene higher than BOTH flanks
        else if(ld!="NA" && rd!="NA" && rl > 1.35 && rr > 1.35) {
            print "POSSIBLE_DUPLICATION_OR_HIGH_COPY"
        }
        else {
            print "NORMAL"
        }
    }')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$gene" "$chrom" "$start" "$end" \
        "$gene_d" "$left_d" "$right_d" \
        "$flank_mean" "$ratio_mean" "$ratio_left" "$ratio_right" \
        "$call" >> "$OUT_ALL"

    if [[ "$call" == "LIKELY_FULL_DELETION_OR_MAJOR_COPY_LOSS" \
       || "$call" == "POSSIBLE_PARTIAL_OR_HET_DELETION" ]]; then
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$gene" "$chrom" "$start" "$end" \
            "$gene_d" "$left_d" "$right_d" \
            "$flank_mean" "$ratio_mean" "$ratio_left" "$ratio_right" \
            "$call" >> "$OUT_HITS"
    fi

done < "$GENES"

# ── Step 5: Results ────────────────────────────────────────────
echo ""
echo "[5/5] Done."
echo ""
echo "Output files:"
printf "  All genes  : %s  (%s data rows)\n" \
    "$OUT_ALL" "$(( $(wc -l < "$OUT_ALL") - 1 ))"
printf "  Hits only  : %s  (%s hits)\n" \
    "$OUT_HITS" "$(( $(wc -l < "$OUT_HITS") - 1 ))"
echo ""
echo "Deletion / copy-loss hits:"
echo ""
if command -v column &>/dev/null; then
    column -t -s $'\t' "$OUT_HITS"
else
    cat "$OUT_HITS"
fi
echo ""
