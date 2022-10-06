#!/bin/bash
#BSUB -q epistasis_long
#BSUB -J interpret_dtox_compound_select_gene_comptoxai_tox21_implementation
#BSUB -n 1
#BSUB -o src/run/interpret_dtox_compound_select_gene_comptoxai_tox21_implementation.%J.out
#BSUB -e src/run/interpret_dtox_compound_select_gene_comptoxai_tox21_implementation.%J.error
#BSUB -N

module unload python
module load python/3.7

./src/run/interpret_dtox_compound_select_gene_comptoxai_tox21_implementation1.sh &
