#!/bin/bash -l
#SBATCH -A sscm013902 -p short --constraint=Cascadelake -J e2e-v9 -c 8 --mem=64G -t 03:00:00 --array=1-3%1
#SBATCH -o /user/work/fh6520/showcase/e2e/logs/e2e-v9-%A_%a.out
# fastMR arms (C1, C8) on lib8 (fastMR fix/clump-flag-preread 671c85f, CompreSSoR d158622) + traits_v8 stores, 1x1 to 50x50.
# Replicates run one at a time. The TwoSampleMR / VCF arms are reused from the v7 runs (same node type).
# The Cascadelake feature also matches Gold 6226 (2.7 GHz) nodes; all v7 cells ran on Gold 6226R, so refuse anything else.
REQUIRE_CPU=${REQUIRE_CPU:-"Gold 6226R"}
grep -m1 'model name' /proc/cpuinfo | grep -q "$REQUIRE_CPU" || { echo "wrong CPU on $(hostname): $(grep -m1 'model name' /proc/cpuinfo)"; exit 3; }
P=/user/work/fh6520/showcase/prep; E=/user/work/fh6520/showcase/e2e
[ -d $P/traits50_v8 ] || bash $E/make_traits50b.sh $P/traits25_v8 $P/traits50_v8
export LIB_CS=/user/work/fh6520/showcase/final/lib8 LIB_FASTMR=/user/work/fh6520/showcase/final/lib8 ARMS="C1 C8" CAP=1500 \
  FASTMR_SHA=671c85f83f2f753a17e4cc79045bca3f31a11634 CS_SHA=d158622004a6f9812f65cf02cd914f4cda7d24d4
for spec in "1x1 10x10:traits_v8" "25x25:traits25_v8" "50x50:traits50_v8"; do
  SIZES=${spec%%:*} TRAITS=$P/${spec##*:} RESDIR=$E/results_v9_${spec##*:} bash $E/e2e.sbatch
done
