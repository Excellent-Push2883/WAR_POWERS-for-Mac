#!/bin/bash
# WAR_POWERS - macOS launcher.  Double-click this file to run.
cd "$(dirname "$0")"
echo "=================================================="
echo "  Executing WAR_POWERS  (macOS)"
echo "=================================================="
if ! command -v samtools >/dev/null 2>&1; then
  echo ""; echo "  samtools is not installed yet."
  if ! command -v brew >/dev/null 2>&1; then
    echo "  1) Install Homebrew (one time) - paste this in Terminal:"
    echo '       /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"'
    echo "  2) Then:  brew install samtools"
  else
    echo "  Install it with:  brew install samtools"
  fi
  echo ""; echo "  Then double-click WAR_POWERS.command again."
  echo ""; echo "Press Enter to close."; read _; exit 1
fi
bash "./WAR_POWERS.sh"
echo ""; echo "Finished. Press Enter to close."; read _
