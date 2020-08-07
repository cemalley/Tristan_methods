#!/bin/sh

module load cellranger/2.1.0

#mkfastq--------------------------

cellranger mkfastq --id=IS006-bcl \
--run=/data/NCATS_ifx/data/Chromium/IS006_180522_J00139_0320_AHV3TVBBXX \
--csv=/data/NCATS_ifx/data/Chromium/sample-csvs/IS006.csv \
--localcores=$SLURM_CPUS_PER_TASK \
--localmem=34


#mkfastq submitted with: sbatch --mem=35g --cpus-per-task=12 --time=48:00:00 scripts/cellranger_pipeline_052918.sh
#2292020 - timeout
# 2344439 - ran out of disk space
# 2371272 - success
#resubmitted june 19 2018 3492496

#count-----------------------------

#the following submitted as three separate shell scripts

#!/bin/sh
module load cellranger/2.1.0
cellranger count --id=IS006-bcl-count-iPSC_Fresh \
--fastqs=/data/NCATS_ifx/data/Chromium/IS006-bcl/outs/fastq_path \
--sample=iPSC_Fresh \
--transcriptome=/fdb/cellranger/refdata-cellranger-GRCh38-1.2.0 \
--localcores=32 --localmem=120 \
--maxjobs=12

#!/bin/sh
module load cellranger/2.1.0
cellranger count --id=IS006-bcl-count-iPSC_Cryo \
--fastqs=/data/NCATS_ifx/data/Chromium/IS006-bcl/outs/fastq_path \
--sample=iPSC_Cryo \
--transcriptome=/fdb/cellranger/refdata-cellranger-GRCh38-1.2.0 \
--localcores=32 --localmem=120 \
--maxjobs=12

#!/bin/sh
module load cellranger/2.1.0
cellranger count --id=IS006-bcl-count-iPSC_MtOH \
--fastqs=/data/NCATS_ifx/data/Chromium/IS006-bcl/outs/fastq_path \
--sample=iPSC_MtOH \
--transcriptome=/fdb/cellranger/refdata-cellranger-GRCh38-1.2.0 \
--localcores=32 --localmem=120 \
--maxjobs=12

sbatch --cpus-per-task=32 --exclusive --mem=120g --constraint=x2650 --time=48:00:00 scripts/cellranger_pipeline_count_iPSC_Fresh_053118.sh
sbatch --cpus-per-task=32 --exclusive --mem=120g --constraint=x2650 --time=48:00:00 scripts/cellranger_pipeline_count_iPSC_Cryo_053118.sh
sbatch --cpus-per-task=32 --exclusive --mem=120g --constraint=x2650 --time=48:00:00 scripts/cellranger_pipeline_count_iPSC_MtOH_053118.sh

#cellranger count batch job IDs
#2405913
#2405914
#2405915 - all three successful

# cellranger aggr

#/data/NCATS_ifx/data/Chromium/sample-csvs/IS006_aggr.csv:
library_id,molecule_h5
iPSC_Fresh,/data/NCATS_ifx/data/Chromium/IS006-bcl-count-iPSC_Fresh/outs/molecule_info.h5
iPSC_Cryo,/data/NCATS_ifx/data/Chromium/IS006-bcl-count-iPSC_Cryo/outs/molecule_info.h5
iPSC_MtOH,/data/NCATS_ifx/data/Chromium/IS006-bcl-count-iPSC_MtOH/outs/molecule_info.h5

#!/bin/bash

module load cellranger/2.1.0

cellranger aggr --id=aggr-IS006-bcl \
--csv=/data/NCATS_ifx/data/Chromium/sample-csvs/IS006_aggr.csv \
--normalize=mapped \
--localcores=32 \
--localmem=120 \
--maxjobs=12

sbatch --cpus-per-task=32 --exclusive --mem=120g --constraint=x2650 --time=24:00:00 scripts/cellranger_pipeline_aggr_060118.sh
#2485876


# convert matrix to dense csv file - done in sinteractive --mem=10g
module load cellranger
#usage: cellranger mat2csv sample123/outs/filtered_gene_bc_matrices sample123.csv

module load cellranger
cellranger mat2csv /data/NCATS_ifx/data/Chromium/IS006-bcl-count-iPSC_Fresh/outs/filtered_gene_bc_matrices /data/NCATS_ifx/iPSC/IS006-iPSC-Fresh-Cryo-MtOH/iPSC_Fresh_dense_expression_matrix.csv

module load cellranger
cellranger mat2csv /data/NCATS_ifx/data/Chromium/IS006-bcl-count-iPSC_Cryo/outs/filtered_gene_bc_matrices /data/NCATS_ifx/iPSC/IS006-iPSC-Fresh-Cryo-MtOH/iPSC_Cryo_dense_expression_matrix.csv

module load cellranger
cellranger mat2csv /data/NCATS_ifx/data/Chromium/IS006-bcl-count-iPSC_MtOH/outs/filtered_gene_bc_matrices /data/NCATS_ifx/iPSC/IS006-iPSC-Fresh-Cryo-MtOH/iPSC_MtOH_dense_expression_matrix.csv

module load cellranger
cellranger mat2csv /data/NCATS_ifx/data/Chromium/aggr-IS006-bcl/outs/filtered_gene_bc_matrices_mex /data/NCATS_ifx/iPSC/IS006-iPSC-Fresh-Cryo-MtOH/iPSC_Fresh_Cryo_MtOH_dense_expression_matrix.csv




#
