#!/bin/bash
#BSUB -q epistasis_long
#BSUB -J dtox_compound_select_gene_comptoxai_binding_tox21_null
#BSUB -n 20
#BSUB -o src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null.%J.out
#BSUB -e src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null.%J.error
#BSUB -N

module unload python
module load python/3.7

./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null1.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null2.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null3.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null4.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null5.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null6.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null7.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null8.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null9.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null10.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null11.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null12.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null13.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null14.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null15.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null16.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null17.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null18.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null19.sh &
./src/run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null20.sh &
