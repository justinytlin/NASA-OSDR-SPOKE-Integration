#!/bin/zsh
# Download the GeneLab RSEM unnormalized counts and ISA metadata used by the tissue check into data/.
# OSD-564's counts are shared with the replication route (psev_data/), so they are copied if present.
set -e
cd "$(dirname "$0")"; mkdir -p data; cd data
B="https://osdr.nasa.gov/geode-py/ws/studies"
dl() { curl -sS -L -o "$2" "$B/$1/download?source=datamanager&file=$2" && echo "ok $2"; }
for n in 244 245 246 288 289 137 161 162 163 168 194 270 352; do
  dl OSD-$n GLDS-${n}_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv &
done
dl OSD-564 GLDS-569_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv &
for n in 244 245 246 288 289 564 137 162 168 194 270 352; do dl OSD-$n OSD-${n}_metadata_OSD-${n}-ISA.zip & done
for n in 161 163; do dl OSD-$n OSD-${n}_metadata_GLDS-${n}-ISA.zip & done
wait
for z in OSD-*-ISA.zip; do n=${z%%_*}; mkdir -p "isa/$n"; unzip -qo "$z" -d "isa/$n"; done
echo "done: $(ls *.csv | wc -l) counts files"
