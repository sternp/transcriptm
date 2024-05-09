import subprocess
import os

gDNA_single_strand = int(snakemake.config["gDNA"])/2
max_memory = str(snakemake.config["max_memory"]) + "G"

#get genes on reverse strand - if no contam then no reads should be present on fwd strand
command = (
    "gff2bed --max-mem %s < annotate/combined_reference.gff > detect_contam/temp.bed && "
    "awk '$8==\"CDS\"' detect_contam/temp.bed > detect_contam/temp.CDS.bed && "
    "awk '$6 == \"-\"' detect_contam/temp.CDS.bed > detect_contam/CDS.rev.bed" % (max_memory)
)
try:
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")




#prepare file with overlapping ORFs. Avoid these regions when filtering contam.
command = (
    "intersectBed -a detect_contam/temp.CDS.bed -b detect_contam/temp.CDS.bed -S > detect_contam/out && "
    "cat detect_contam/out | awk '{{print $1 \"\\t\" $2 \"\t\" $3 \"\t\" $8}}' | sort -u > detect_contam/overlappedORFs.bed ; rm detect_contam/out"
)
try:
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")



#split the reads mapping to forward strand
command = (
    "samtools view -b -f 128 -F 16 coverm_filter/combined_reference_filtered.bam > detect_contam/fwd1.bam && "
    "samtools view -b -f 64 -F 32 coverm_filter/combined_reference_filtered.bam > detect_contam/fwd2.bam && "
    "samtools merge -f detect_contam/fwd.bam detect_contam/fwd1.bam detect_contam/fwd2.bam && "
    "samtools index detect_contam/fwd.bam"
)
try:
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")




#identify genomes with antisense reads. calculates median antisense coverage across all genes in a genome (fwd strand)...
#...then extracts all contigs per genome for use in strandCheckR
command = (
    "bedtools coverage -hist -b detect_contam/fwd.bam -a detect_contam/CDS.rev.bed > detect_contam/temp1 && "
    "sed -i '/^all/d' detect_contam/temp1 && "
    "awk -F\"\t\" '{{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $11 * $12 / $13}}' detect_contam/temp1 > detect_contam/temp2 && "
    "cat detect_contam/temp2 | datamash -s -g1,2,3 sum 4 > detect_contam/temp3 && "  # sum the values from bedtools hist to get the total coverage
    "cat detect_contam/temp3 | awk '{{print $1}}' | rev | cut -d '_' -f2- | rev > detect_contam/temp3_genomeID && "  # get genome name, trims the last '_' separated field which represents contig number.
    "paste detect_contam/temp3_genomeID detect_contam/temp3 > detect_contam/temp4 && "
    "sort -k1,1 detect_contam/temp4 > detect_contam/temp4_sort && "
    "cat detect_contam/temp4_sort | datamash -s -g1 median 5 > detect_contam/temp5 && "
    "awk -F\"\t\" '$2 >= %s {{print $1}}' detect_contam/temp5 > detect_contam/contaminated_genome_list && "
    "sort -k1,1 detect_contam/contaminated_genome_list > detect_contam/contaminated_genome_list_sort && "
    "join detect_contam/contaminated_genome_list_sort detect_contam/temp4_sort > detect_contam/temp6 && "
    "cat detect_contam/temp6 | awk '{{print $2}}' | sort -u > detect_contam/contaminated_contigs_list" % (gDNA_single_strand)
)
try:
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")



#split the contaminated contigs into separate bam for strandCheckR
command = (
	"samtools view coverm_filter/combined_reference_filtered.bam | awk '{{print $3}}' | sort -u > detect_contam/all_contigs && "
	"comm -3 detect_contam/all_contigs detect_contam/contaminated_contigs_list > detect_contam/non-contaminated_contigs_list && "
	"sed -i \'s/^\t//\' detect_contam/non-contaminated_contigs_list"
)
try:
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")


#create beds for splitting
command = (
	"samtools faidx annotate/combined_reference.fna && "
	"awk -v FS=\"\t\" -v OFS=\"\t\" '{{print $1 FS \"0\" FS ($2-1)}}' annotate/combined_reference.fna.fai > annotate/combined_reference.bed && "
	"fgrep -w -f detect_contam/contaminated_contigs_list annotate/combined_reference.bed > detect_contam/contaminated_contigs.bed && "
	"fgrep -w -f detect_contam/non-contaminated_contigs_list annotate/combined_reference.bed > detect_contam/non-contaminated_contigs.bed"
)
try:
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
#    print(f"An error occurred: {e}")
	":"
command = (
	"samtools view -bh -L detect_contam/contaminated_contigs.bed coverm_filter/combined_reference_filtered.bam > detect_contam/contam_contigs.bam && "
	"samtools view -bh -L detect_contam/non-contaminated_contigs.bed coverm_filter/combined_reference_filtered.bam > detect_contam/non-contam_contigs.bam && "
	"samtools index detect_contam/contam_contigs.bam && "
	"samtools index detect_contam/non-contam_contigs.bam"
)
try:
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")



#finish
command = (
	"touch detect_contam/done"
)
try:
    subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
    print(f"An error occurred: {e}")
