#!/bin/zsh
# All 104 control sets on the 2025 'upregulates-only' PSEV build, in three batches, each followed by its checks.
cd "$(dirname "$0")"
TAG=2025upAuG; B=/Volumes/justinytlin/SPOKE_psev/gene_psev_spoke2025-03_upAuG_only
for b in A B C; do
  names=$(python3 -c "import json; print(','.join(json.load(open('data/_catalog/batches_$TAG.json'))['$b']))")
  echo "=== batch $b $(date)"
  python3 scripts/gc_counts_to_psev.py --manifest data/_catalog/manifest_$TAG.json --names "$names" --compact --skip-done \
      --zip "$B" --spoke "$B" --threads 6 --chunk 64 2>&1 | grep --line-buffered -v -i warn | grep --line-buffered -E "datasets to run|^done|Traceback|Error|pass2 group 19"
  find /Volumes/justinytlin/SPOKE_psev -name "._*" -delete 2>/dev/null
  python3 scripts/check_build.py --tag $TAG --build "$B" --names "$names" 2>&1 | grep -v -i warn
done
echo "=== all batches done $(date)"
