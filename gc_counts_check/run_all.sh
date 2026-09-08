#!/bin/zsh
# Reproduce everything: counts -> 2019-PSEV rank tables (needs psev_data/gene_psev.zip), then one
# anatomy check per dataset (each dataset judged on its own; nothing is compared across datasets).
set -e
cd "$(dirname "$0")"
python3 scripts/gc_counts_to_psev.py --threads 4 | tee psev_out/run.log   # all 14 datasets in one pass (~1 h); --only RR3 for the RR-3 set
for d in "OSD-564 hippocampus" "OSD-244 thymus" "OSD-289 thymus" "OSD-246 spleen" "OSD-288 spleen" "OSD-245 liver" \
         "OSD-137 liver" "OSD-161 adrenal gland" "OSD-162 eye" "OSD-163 kidney" "OSD-168 liver" "OSD-194 retina" "OSD-270 heart" "OSD-352 brain"; do
  python3 scripts/check_anatomy.py --ranks "psev_out/${d%% *}_GC_counts_ranks_and_rank_by_type_for_meta0_1.tsv" --tissue "${d#* }" --name "${d%% *}"
done
