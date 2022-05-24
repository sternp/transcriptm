#check for dependencies and install
if("BiocManager" %in% rownames(installed.packages()) == FALSE) {install.packages("BiocManager", repos='http://cran.us.r-project.org',quiet = TRUE)}
if("regioneR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("regioneR")}
if("ballgown" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ballgown")}

#load libraries
suppressMessages(library(strandCheckR))
suppressMessages(library(ballgown))
suppressMessages(library(GenomicRanges))
suppressMessages(library(regioneR))

##############################################################################

#read overlapping ORF regions
overlap_list <- toGRanges(snakemake@params[[5]])

#read list of genomes with contamination
contam_list = scan(snakemake@params[[4]],what="", sep="\n")

#filter, if contam exists.
if (length(contam_list) > 0){
table=filterDNA(file = snakemake@params[[1]],
  destination = snakemake@params[[2]],
  threshold = 0.9,
  winWidth = 200,
  winStep = 50,
  getWin = TRUE,
  mustKeepRanges = overlap_list)
} else {
  system('echo "No gDNA detected, skipping filtering step..."')
  system('ln -s $(pwd)/coverm_filter/combined_reference_filtered.bam $(pwd)/final_bam/final.bam')
}

#finish
system("touch filter_contam/done")
if (length(contam_list) > 0){
system("mv $(pwd)/out.stat $(pwd)/decontamination_stats.txt")
}
