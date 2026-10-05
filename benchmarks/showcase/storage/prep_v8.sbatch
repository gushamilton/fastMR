#!/bin/bash -l
#SBATCH -A sscm013902 -p short --constraint=Cascadelake -J stor-v8prep -c 8 --mem=64G -t 01:00:00
#SBATCH -o /user/work/fh6520/showcase/storage/logs/prep-v8-%j.out
# Re-write stores/finngen.cpr with lib7 (prep + fidelity). The v7 store is kept as stores/finngen_v7.cpr for the
# byte comparison, and the v7 sizes / prep / fidelity CSVs are copied to *_v7.csv first (prep.R overwrites them).
module load languages/R/4.5.1 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export STORAGE_ROOT=/user/work/fh6520/showcase/storage STORAGE_SCRIPTS=/user/work/fh6520/showcase/storage/scripts LIB_CS=/user/work/fh6520/showcase/final/lib7
R=$STORAGE_ROOT/results; St=$STORAGE_ROOT/stores
for f in sizes prep_writes_cpr fidelity_exactness fidelity_errors fidelity_flips fidelity_pvalue_order fidelity_cpr_nonfinite_rows; do
  [ -e $R/$f.csv ] && [ ! -e $R/${f}_v7.csv ] && cp $R/$f.csv $R/${f}_v7.csv; done
[ -e $St/finngen_v7.cpr ] || mv $St/finngen.cpr $St/finngen_v7.cpr
PREP_FORMATS=cpr Rscript $STORAGE_SCRIPTS/prep.R
for f in $(ls $St/finngen.cpr); do printf '%s,%s,%s,' $f $(stat -c %s $St/finngen_v7.cpr/$f 2>/dev/null || echo NA) $(stat -c %s $St/finngen.cpr/$f)
  cmp -s $St/finngen_v7.cpr/$f $St/finngen.cpr/$f && echo identical || echo differs; done | sed '1i file,bytes_v7,bytes_v8,content' > $R/store_bytes_v7_vs_v8.csv
cat $R/store_bytes_v7_vs_v8.csv
