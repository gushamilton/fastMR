#!/bin/bash -l
#SBATCH -A sscm013902 -p test --constraint=Cascadelake -J e2e-ab -c 8 --mem=64G -t 00:40:00 --array=1-3%1
#SBATCH -o /user/work/fh6520/showcase/e2e/logs/e2e-ab-%A_%a.out
# Same-node A/B control: the fastMR C8 arm on the v7 build (lib6 + traits*_v7) and the v8 build (lib7 + traits*_v8),
# interleaved in random order within one job, so node load and filesystem state are shared by both builds.
# BUILDS (default "v7 v8") may add v9 = lib8 (fastMR fix/clump-flag-preread 671c85f + CompreSSoR d158622);
# AB_DIR names the results dir (default results_ab).
set -uo pipefail
REQUIRE_CPU=${REQUIRE_CPU:-"Gold 6226R"}
grep -m1 'model name' /proc/cpuinfo | grep -q "$REQUIRE_CPU" || { echo "wrong CPU on $(hostname)"; exit 3; }
module load languages/R/4.5.1 apps/plink2/2.00a68LM apps/plink1.9/1.90-b77
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
E=/user/work/fh6520/showcase/e2e; P=/user/work/fh6520/showcase/prep; F=/user/work/fh6520/showcase/final
REP=${SLURM_ARRAY_TASK_ID:-1}; BUILDS=${BUILDS:-"v7 v8"}; RES=$E/${AB_DIR:-results_ab}/rep$REP; mkdir -p $RES
export PLINK1=$(which plink) PLINK2=$(which plink2) SEED=$((1000+REP))
CELLS=$RES/cells.csv; echo "rep,size,build,arm,threads,status,wall_s,user_s,sys_s,maxrss_kb,host" > $CELLS
for round in 1 2; do
 for size in 10x10 25x25; do
  set=traits; [ $size = 25x25 ] && set=traits25
  for b in $(printf '%s\n' $BUILDS | shuf --random-source=<(yes $((REP*13+round)))); do
    case $b in
      v7) lib=$F/lib6 tg=v7 fsha=fb62057 csha=f91bb8d;;
      v8) lib=$F/lib7 tg=v8 fsha=543cae282d0320500af75728c2c71550da75d32e csha=d158622004a6f9812f65cf02cd914f4cda7d24d4;;
      v9) lib=$F/lib8 tg=v8 fsha=671c85f83f2f753a17e4cc79045bca3f31a11634 csha=d158622004a6f9812f65cf02cd914f4cda7d24d4;;
    esac
    d=$RES/${size}_C8_${b}_r$round; mkdir -p $d
    FASTMR_SHA=$fsha CS_SHA=$csha LIB_CS=$lib LIB_FASTMR=$lib TRAITS=$P/${set}_$tg /usr/bin/time -v -o $d/time.txt timeout 1500 Rscript $E/arm.R C ${size%x*} ${size#*x} 8 $d > $d/log.txt 2>&1
    st=$?
    w=$(awk -F': ' '/Elapsed \(wall/{n=split($NF,p,":"); s=0; for(i=1;i<=n;i++) s=s*60+p[i]; print s}' $d/time.txt)
    u=$(awk -F': ' '/User time/{print $2}' $d/time.txt); sy=$(awk -F': ' '/System time/{print $2}' $d/time.txt)
    rss=$(awk -F': ' '/Maximum resident/{print $2}' $d/time.txt)
    echo "$REP,$size,$b,C8,8,$st,$w,$u,$sy,$rss,$(hostname)" >> $CELLS
  done
 done
done
echo finished
