import subprocess
import os

genome_dir = snakemake.config["genome_dir"]
fna  = snakemake.config["ref"]
gff = snakemake.config["gff"]
genome_dir_path = os.path.join(snakemake.config["genome_dir"], "*." + snakemake.config["fasta_extension"])
kingdom = snakemake.config["kingdom"]

if genome_dir != "none" and fna == "none" and gff == 'none':
	subprocess.Popen(
		"echo 'Running annotation (time consuming), please standby...' && "
		"cat %s > annotate/combined_reference.fna && "
		"ls %s | parallel -j %s prokka --metagenome --cpus 1 --kingdom %s {} --outdir annotate/{/.}.prokka --force && "
#		"find . -name \'*.gff\' -exec cat {} + > annotate/out ; mv annotate/out annotate/combined_reference.gff && "
		"find . -type f -name \'*.gff\' -exec gffread -F -o- {} + > annotate/combined_reference.gff && "
		"sed -i \'s/Parent=/ID=/\' annotate/combined_reference.gff && "
		"touch annotate/done" %
		(genome_dir_path, genome_dir_path, snakemake.threads, kingdom), shell=True).wait()


elif genome_dir == "none" and fna != "none" and gff != 'none':
	subprocess.Popen(
		"ln -s %s annotate/combined_reference.fna && "
		"ln -s %s annotate/combined_reference.gff && "
		"touch annotate/done" %
		(fna, gff), shell=True).wait()


elif genome_dir == "none" and fna == "none" and gff == 'none':
	subprocess.Popen(
		"echo ERROR: Provide either 1\) directory of contigs/genomes to annotate \(-g\), or 2\) a reference sequence file \(--ref\) \(.fna\) and an annotation file \(--gff\) \(.gff\)", shell=True).wait()


elif genome_dir != "none" and fna != "none" or gff != 'none':
	subprocess.Popen(
		"echo ERROR: Provide either 1\) directory of contigs/genomes to annotate \(-g\), or 2\) a reference sequence file \(--ref\) \(.fna\) and an annotation file \(--gff\) \(.gff\)", shell=True).wait()
