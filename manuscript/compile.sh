#!/usr/bin/env bash
set -euo pipefail

compile_one() {
  local base="$1"
  xelatex -interaction=nonstopmode -halt-on-error "${base}.tex"
  if grep -q '\\citation' "${base}.aux"; then
    if command -v bibtex.original >/dev/null 2>&1; then
      bibtex.original "$base"
    else
      bibtex "$base"
    fi
  fi
  xelatex -interaction=nonstopmode -halt-on-error "${base}.tex"
  xelatex -interaction=nonstopmode -halt-on-error "${base}.tex"
}

compile_one sjs_profile_quantile_main
compile_one sjs_profile_quantile_supplement
