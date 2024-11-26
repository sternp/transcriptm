#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

os.system('rm -r final_table')
os.mkdir('final_table')

infile = open(snakemake.params.in_table)
next(infile)
rates=[]
for line in infile :
    columns = line.rstrip("\n").split("\t")
    length = float(columns[5])
    count = float(columns[12])
    rate = count / length
    rates.append(rate)
infile.close()
sum_rate = sum(rates)


# calculate TPM for the real count table
infile = open(snakemake.params.in_table)
next(infile)
with open(snakemake.params.out_table, "w") as outfile:
    outfile.write("gene_id" + "\t" + "contig" + "\t" + "start" + "\t" + "end" + "\t" + "strand"  + "\t" + "length"  + "\t" + "gene"  + "\t" + "product" + "\t" "db_xref" + "\t" "Dbxref" + "\t" + "inference" + "\t" + "UniProtKB" + "\t" + "raw_count" + "\t" + "rpk"  + "\t" + "tpm" + "\n")
    for line in infile :
        columns = line.rstrip("\n").split("\t")
        length = float(columns[5])
        count = float(columns[12])
        rate = count / length
        tpm = (rate / sum_rate) *1000000
        rpk = count / (length / 1000)
        print(str(columns[0]) + "\t" +
              str(columns[1]) + "\t" +
              str(columns[2]) + "\t" +
              str(columns[3]) + "\t" +
              str(columns[4]) + "\t" +
              str(columns[5]) + "\t" +
              str(columns[6]) + "\t" +
              str(columns[7]) + "\t" +
              str(columns[8]) + "\t" +
              str(columns[9]) + "\t" +
              str(columns[10]) + "\t" +
              str(columns[11]) + "\t" +
              str(columns[12]) + "\t" +
              str(rpk) + "\t" +
              str(tpm), file=outfile)
infile.close()
outfile.close()


# calculate TPM for the intergenic regions
infile = open(snakemake.params.in_table_intergenic)
next(infile)
rates=[]
for line in infile :
    columns = line.rstrip("\n").split("\t")
    length = float(columns[5])
    count = float(columns[6])
    rate = count / length
    rates.append(rate)
infile.close()
sum_rate = sum(rates)

infile = open(snakemake.params.in_table_intergenic)
next(infile)
with open(snakemake.params.out_table_intergenic, "w") as outfile:
    outfile.write("gene_id" + "\t" + "contig" + "\t" + "start" + "\t" + "end" + "\t" + "strand"  + "\t" + "length"  + "\t" + "raw_count"  + "\t" + "rpk" + "\t" + "tpm" + "\n")
    for line in infile :
        columns = line.rstrip("\n").split("\t")
        length = float(columns[5])
        count = float(columns[6])
        rate = count / length
        tpm = (rate / sum_rate) *1000000
        rpk = count / (length / 1000)
        print(str(columns[0]) + "\t" +
              str(columns[1]) + "\t" +
              str(columns[2]) + "\t" +
              str(columns[3]) + "\t" +
              str(columns[4]) + "\t" +
              str(columns[5]) + "\t" +
              str(columns[6]) + "\t" +
              str(rpk) + "\t" +
              str(tpm), file=outfile)
infile.close()
outfile.close()



############################################################################################3
# Reapeat for the forward strand - only for non-Illumina RNAseq

infile = open(snakemake.params.in_table_fwd)
next(infile)
rates=[]
for line in infile :
    columns = line.rstrip("\n").split("\t")
    length = float(columns[5])
    count = float(columns[12])
    rate = count / length
    rates.append(rate)
infile.close()
sum_rate = sum(rates)


# calculate TPM for the real count table
infile = open(snakemake.params.in_table_fwd)
next(infile)
with open(snakemake.params.out_table_fwd, "w") as outfile:
    outfile.write("gene_id" + "\t" + "contig" + "\t" + "start" + "\t" + "end" + "\t" + "strand"  + "\t" + "length"  + "\t" + "gene"  + "\t" + "product" + "\t" "db_xref" + "\t" "Dbxref" + "\t" + "inference" + "\t" + "UniProtKB" + "\t" + "raw_count" + "\t" + "rpk"  + "\t" + "tpm" + "\n")
    for line in infile :
        columns = line.rstrip("\n").split("\t")
        length = float(columns[5])
        count = float(columns[12])
        rate = count / length
        tpm = (rate / sum_rate) *1000000
        rpk = count / (length / 1000)
        print(str(columns[0]) + "\t" +
              str(columns[1]) + "\t" +
              str(columns[2]) + "\t" +
              str(columns[3]) + "\t" +
              str(columns[4]) + "\t" +
              str(columns[5]) + "\t" +
              str(columns[6]) + "\t" +
              str(columns[7]) + "\t" +
              str(columns[8]) + "\t" +
              str(columns[9]) + "\t" +
              str(columns[10]) + "\t" +
              str(columns[11]) + "\t" +
              str(columns[12]) + "\t" +
              str(rpk) + "\t" +
              str(tpm), file=outfile)
infile.close()
outfile.close()


# calculate TPM for the intergenic regions
infile = open(snakemake.params.in_table_fwd_intergenic)
next(infile)
rates=[]
for line in infile :
    columns = line.rstrip("\n").split("\t")
    length = float(columns[5])
    count = float(columns[6])
    rate = count / length
    rates.append(rate)
infile.close()
sum_rate = sum(rates)

infile = open(snakemake.params.in_table_fwd_intergenic)
next(infile)
with open(snakemake.params.out_table_fwd_intergenic, "w") as outfile:
    outfile.write("gene_id" + "\t" + "contig" + "\t" + "start" + "\t" + "end" + "\t" + "strand"  + "\t" + "length"  + "\t" + "raw_count"  + "\t" + "rpk" + "\t" + "tpm" + "\n")
    for line in infile :
        columns = line.rstrip("\n").split("\t")
        length = float(columns[5])
        count = float(columns[6])
        rate = count / length
        tpm = (rate / sum_rate) *1000000
        rpk = count / (length / 1000)
        print(str(columns[0]) + "\t" +
              str(columns[1]) + "\t" +
              str(columns[2]) + "\t" +
              str(columns[3]) + "\t" +
              str(columns[4]) + "\t" +
              str(columns[5]) + "\t" +
              str(columns[6]) + "\t" +
              str(rpk) + "\t" +
              str(tpm), file=outfile)
infile.close()
outfile.close()






os.system('touch final_table/done')
#os.system('ln -s final_table/rawcount_tpm_table.txt rawcount_tpm_table.txt')
