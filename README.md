# TranscriptM

A  metatranscriptome bioinformatics pipeline including metagenome contamination correction



## Installation (conda)
```
git clone git@github.com:sternp/transcriptm.git

cd transcriptm

conda env create -n transcriptm -f transcriptm.yaml

conda activate transcriptm

pip install -e .
```


## Usage (QC, mapping, gDNA decontamination, and read counting pipeline)
```
transcriptm count
    -1 reads_R1.fastq.gz \
    -2 reads_R2.fastq.gz \
    -n 24 \
    --ref combined_reference.fna \
    --gff combined_reference.gff \
     -m 128 \
    -db /dir/to/bowtie2_1 /dir/to/bowtie2_2 \
    -o output_directory

Specifying -g will concatenate a directory of .fna genomes into a single ref sequence and annotate with prokka (time intensive)
Alternatively, you can use pre-contructed files via --ref <combined_reference.fna> and --gff <combined_reference.gff>

Please note, currently the required format for the contig's FASTA headers are <genomeID>_<number>.
For example: >Ardenticatenaceae-ID1234_00001, >Ardenticatenaceae-ID1234_00002...etc
```

## Notes
You may need to download databases to filter contaminating reads from rRNA genes, the human genome...etc. Otherwise you can make your own Bowtie2-formatted database. You can specify multiple Bowtie2 databases in the transcriptm command, i.e. ```-db /dir/to/db1 /dir/to/db2```

Pre-made databases can be downloaded like so:
```
conda activate transcriptm

kneaddata_database --download human_genome bowtie2 /work/microbiome/db/human_grch37_bowtie2
```
Available pre-made databases:
```
KneadData Databases ( database : build = location )
human_genome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz
human_genome : bmtagger = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_BMTagger_v0.1.tar.gz
human_transcriptome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz
ribosomal_RNA : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.2.tar.gz
mouse_C57BL : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz

```
