import subprocess
import os

if snakemake.config["genome_dir"] != "none" and snakemake.config["ref"] == "none" and snakemake.config["gff"] == 'none':
	genome_dir_path = os.path.join(snakemake.config["genome_dir"], "*." + snakemake.config["fasta_extension"])
	gff_dir_path = os.path.join(snakemake.config["output"], "annotate","prokka_out")
	gff_file_path = os.path.join(gff_dir_path,"*.gff")
	subprocess.Popen(
		"cat %s > annotate/combined_reference.fna && "
		"prokka --metagenome --cpus %s --kingdom Bacteria annotate/combined_reference.fna --outdir %s --force && "
		"cp %s annotate/combined_reference.gff && "
		"touch annotate/done" %
		(genome_dir_path, snakemake.threads, gff_dir_path, gff_file_path), shell=True).wait()

elif snakemake.config["genome_dir"] == "none" and snakemake.config["ref"] != "none" and snakemake.config["gff"] != 'none':
	fna_file_path = snakemake.config["ref"]
	gff_file_path = snakemake.config["gff"]
	subprocess.Popen(
		"cp %s annotate/combined_reference.fna && "
		"cp %s annotate/combined_reference.gff && "
		"touch annotate/done" %
		(fna_file_path, gff_file_path), shell=True).wait()

elif snakemake.config["genome_dir"] == "none" and snakemake.config["ref"] == "none" and snakemake.config["gff"] == 'none':
	subprocess.Popen(
		"echo ERROR: Provide either 1\) directory of contigs/genomes to annotate \(-g\), or 2\) a reference sequence file \(--ref\) \(.fna\) and an annotation file \(--gff\) \(.gff\)", shell=True).wait()

elif snakemake.config["genome_dir"] != "none" and snakemake.config["ref"] != "none" or snakemake.config["gff"] != 'none':
	subprocess.Popen(
		"echo ERROR: Provide either 1\) directory of contigs/genomes to annotate \(-g\), or 2\) a reference sequence file \(--ref\) \(.fna\) and an annotation file \(--gff\) \(.gff\)", shell=True).wait()
