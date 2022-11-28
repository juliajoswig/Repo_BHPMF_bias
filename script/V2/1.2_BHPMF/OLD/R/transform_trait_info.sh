#!/bin/sh

# this script is a workaround because I forgot to transform the input data
# this should be fixed in future versions

alias R='/usr/local/apps/R/R-3.1.1/bin/R'
alias Rscript='/usr/local/apps/R/R-3.1.1/bin/Rscript'

bsub -q mpi_large \
     -m io1 \
     -M 4 \
     -n 1 \
     -oo ../data/output/transform_trait_info_lsf_std.out \
     -eo ../data/output/transform_trait_info_lsf_std.err \
     Rscript transform_trait_info.R
