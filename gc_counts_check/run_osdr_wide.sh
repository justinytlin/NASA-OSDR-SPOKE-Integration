#!/bin/zsh
# OSDR-wide run: every mouse bulk RNA-seq study with GeneLab counts and a control group.
#   scripts/survey_osdr.py     -> data/_catalog/survey.tsv (studies, factors, files)   [also downloads ISA metadata]
#   scripts/build_manifest.py  -> data/_catalog/manifest.json (control samples per study and tissue)
# then, per batch (A, B = ground control; C = analog-study controls): projection + one anatomy check per entry,
# and finally the first 14 datasets again in compact form so every dataset has scores + the uniform-weight profile.
cd "$(dirname "$0")"
for b in ${=1:-B C A}; do
  names=$(python3 -c "import json; print(','.join(json.load(open('data/_catalog/batches.json'))['$b']))")
  echo "=== batch $b $(date)"
  python3 scripts/gc_counts_to_psev.py --manifest data/_catalog/manifest.json --names "$names" --compact --skip-done --threads 7 --chunk 64 2>&1 | grep --line-buffered -v -i warn
  python3 scripts/check_manifest.py --names "$names"
done
echo "=== batch legacy $(date)"
python3 scripts/gc_counts_to_psev.py --compact --skip-done --threads 7 --chunk 64 2>&1 | grep --line-buffered -v -i warn
echo "=== all batches done $(date)"
