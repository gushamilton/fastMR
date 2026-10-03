#!/bin/bash
#SBATCH --job-name=sc-e2e
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --constraint=Cascadelake
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:30:00
#SBATCH --array=1-3
#SBATCH --output=/user/work/fh6520/showcase/e2e/logs/e2e-%A_%a.out
# One replicate per array task.  Arms run in randomised order (seeded by the task id) inside the task.
# Override with e.g. --export=ALL,SIZES="1x1",RESDIR=...,LIB_CS=...,LIB_FASTMR=...,LIB_TSMR=...
set -uo pipefail
module load languages/R/4.5.1 apps/plink2/2.00a68LM apps/plink1.9/1.90-b77 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
E2E=/user/work/fh6520/showcase/e2e; REP=${SLURM_ARRAY_TASK_ID:-1}
SIZES=${SIZES:-"1x1 10x10"}; CAP=${CAP:-1500}
RESDIR=${RESDIR:-$E2E/results}/rep$REP; mkdir -p $RESDIR
export PLINK1=$(which plink) PLINK2=$(which plink2) SEED=$((1000+REP))
CELLS=$RESDIR/cells.csv; echo "rep,size,arm,threads,status,censored,wall_s,user_s,sys_s,maxrss_kb,host" > $CELLS
declare -A CENS
for size in $SIZES; do
  E=${size%x*}; O=${size#*x}
  arms="A B C1 C8"; [ "$size" != "1x1" ] && arms="A A8 B C1 C8"
  [ -n "${ARMS:-}" ] && arms="$ARMS"
  for arm in $(printf '%s\n' $arms | shuf --random-source=<(yes $((REP*7+${#size})))); do
    d=$RESDIR/${size}_$arm; mkdir -p $d
    if [ "${CENS[$arm]:-0}" = 1 ]; then
      echo "$REP,$size,$arm,NA,skipped,predicted,NA,NA,NA,NA,$(hostname)" >> $CELLS; continue; fi
    case $arm in A) a=A; th=8;; A8) a=A8; th=8;; B) a=B; th=8;; C1) a=C; th=1;; C8) a=C; th=8;; esac
    /usr/bin/time -v -o $d/time.txt timeout $CAP Rscript $E2E/arm.R $a $E ${O} $th $d > $d/log.txt 2>&1
    st=$?
    w=$(awk -F': ' '/Elapsed \(wall/{n=split($NF,p,":"); s=0; for(i=1;i<=n;i++) s=s*60+p[i]; print s}' $d/time.txt)
    u=$(awk -F': ' '/User time/{print $2}' $d/time.txt); sy=$(awk -F': ' '/System time/{print $2}' $d/time.txt)
    rss=$(awk -F': ' '/Maximum resident/{print $2}' $d/time.txt)
    cen=0; [ $st = 124 ] && { cen=1; CENS[$arm]=1; }
    echo "$REP,$size,$arm,$th,$st,$cen,$w,$u,$sy,$rss,$(hostname)" >> $CELLS
  done
done
echo finished
