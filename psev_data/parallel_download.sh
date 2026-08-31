#!/bin/bash
# Parallel chunked download of gene_psev.zip from Zenodo (resumable).
set -u
cd "$(dirname "$0")"
URL="https://zenodo.org/records/4404618/files/gene_psev.zip?download=1"
TOTAL=25643336977
N=8
CHUNK=$(( (TOTAL + N - 1) / N ))

# seed part_0 from any partial single-stream download
if [ -f gene_psev.zip ] && [ ! -f part_0 ] && [ "$(stat -f%z gene_psev.zip)" -lt "$TOTAL" ]; then
  mv gene_psev.zip part_0
fi

pids=()
for i in $(seq 0 $((N-1))); do
  start=$(( i * CHUNK ))
  end=$(( start + CHUNK - 1 )); [ $end -ge $TOTAL ] && end=$(( TOTAL - 1 ))
  want=$(( end - start + 1 ))
  (
    while :; do
      have=0; [ -f part_$i ] && have=$(stat -f%z part_$i)
      [ "$have" -ge "$want" ] && break
      curl -sS --retry 3 --retry-delay 5 -r $((start+have))-$end -o - "$URL" >> part_$i
    done
  ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done

# verify and assemble
for i in $(seq 0 $((N-1))); do
  start=$(( i * CHUNK ))
  end=$(( start + CHUNK - 1 )); [ $end -ge $TOTAL ] && end=$(( TOTAL - 1 ))
  want=$(( end - start + 1 ))
  have=$(stat -f%z part_$i)
  if [ "$have" -ne "$want" ]; then echo "PART $i SIZE MISMATCH $have != $want"; exit 1; fi
done
cat part_0 part_1 part_2 part_3 part_4 part_5 part_6 part_7 > gene_psev.zip
rm -f part_*
have=$(stat -f%z gene_psev.zip)
[ "$have" -eq "$TOTAL" ] || { echo "FINAL SIZE MISMATCH $have"; exit 1; }
python3 -c "import zipfile; zipfile.ZipFile('gene_psev.zip').testzip(); print('ZIP_OK')" || exit 1
echo DOWNLOAD_COMPLETE
