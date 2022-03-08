# TranscriptM
Under development. A metatranscriptome pipeline with gDNA detection and removal.

## Installation
```
git clone git@github.com:sternp/transcriptm.git

cd transcriptm

mamba env create -n transcriptm -f /home/sternesp/lustre-microbiome/users/sternesp/transcriptm/transcriptm.yaml

conda activate transcriptm

cd /home/sternesp/lustre-microbiome/users/sternesp/transcriptm/

pip install -e .
```

## Notes
#read QC
You may need to download database to filter contaminating reads from rRNA genes, the human genome...etc

```
kneaddata_database --download human_genome bowtie2 /lustre/work-lustre/microbiome/db/human_grch37_bowtie2
kneaddata_database --download ribosomal_RNA bowtie2 /lustre/work-lustre/microbiome/db/silva_rRNA_bowtie2
kneaddata_database --download human_transcriptome bowtie2 /lustre/work-lustre/microbiome/db/human_transcriptome_bowtie2
```

These are already downloaded on the CMR server. If you wish to filter against these databases then they need to be specified by TranscriptM as below. You can specify as many of DBs as you want:

```
transcriptm \
-db /lustre/work-lustre/microbiome/db/human_grch37_bowtie2 -db /lustre/work-lustre/microbiome/db/silva_rRNA_bowtie2 \
```
