export OMP_NUM_THREADS=1

if [ ! -f ../data/output/mean_gap_filled.txt ];
then
  echo "$PWD/gapfill.sh --> starting job"
  bsub -q mpi_large \
       -n 1 \
       -m io[3,4] \
       -M 9 \
       -R "span[hosts=1]" \
       -J hpmf_gapfill \
       -oo ../data/output/gapfill_lsf_std.out \
       -eo ../data/output/gapfill_lsf_std.err \
       Rscript gapfill.R
else
  echo "$PWD/../data/output/mean_gap_filled.txt already exist"
fi
