#!/bin/bash
cd "$(dirname "$0")"
if ! command -v samtools >/dev/null 2>&1; then
  echo "samtools not installed. Install Homebrew then:  brew install samtools"
  echo "Press Enter to close."; read _; exit 1
fi
bash "./WAR_POWERS.sh"
echo ""; echo "Finished. Press Enter to close."; read _
