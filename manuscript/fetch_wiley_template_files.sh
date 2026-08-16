#!/usr/bin/env bash
set -euo pipefail

# Convenience retrieval of non-font Wiley NJDv5 support files from a public mirror.
# For submission, compare them with the current template distributed by Wiley.
base="https://raw.githubusercontent.com/ramiromagno/wiley-njd/main/_extensions/wiley-njd/wiley-njd-v5"
curl -fL "${base}/WileyNJDv5.cls" -o WileyNJDv5.cls
curl -fL "${base}/NJDapacite.sty" -o NJDapacite.sty
curl -fL "${base}/wileyNJD-APA.bst" -o WileyNJD-APA.bst
printf '%s\n' "Fetched WileyNJDv5.cls, NJDapacite.sty and WileyNJD-APA.bst."
printf '%s\n' "Font files are not included. Obtain the current Wiley NJDv5 package from Wiley."
