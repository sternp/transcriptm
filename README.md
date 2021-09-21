# transcriptM
Under development

## installation

mamba env create -n transcriptM -f /home/sternesp/lustre-microbiome/users/sternesp/transcriptM/transcriptM.yaml

conda activate transcriptM

cd /home/sternesp/lustre-microbiome/users/sternesp/transcriptM/

pip install -e .

optional, change the version of conda used:
conda env config vars set CONDA_ENV_PATH=/lustre/work-lustre/microbiome/users/sternesp/conda/envs/

##notes
Set your conda_prefix in transcriptM/template_config.yaml
(Optional) set your other parameters in transcriptM/template_config.yaml

##read QC
You may need to download database to filter contaminating reads from rRNA genes, the human genome...etc

kneaddata_database --download human_genome bowtie2 /lustre/work-lustre/microbiome/db/human_grch37_bowtie2
kneaddata_database --download ribosomal_RNA bowtie2 /lustre/work-lustre/microbiome/db/silva_rRNA_bowtie2

These are already downloaded on the CMR server. If you wish to filter against these databases then they need to be specified by transcriptM as below. You can specify as many of DBs as you want:
transcriptM \
--qc -db /lustre/work-lustre/microbiome/db/human_grch37_bowtie2 -db /lustre/work-lustre/microbiome/db/silva_rRNA_bowtie2 \



