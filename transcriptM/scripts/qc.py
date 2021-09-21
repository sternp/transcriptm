import subprocess
import os

short_reads_1 = snakemake.config["short_reads_1"]
short_reads_2 = snakemake.config["short_reads_2"]
output_dir_path = os.path.join(snakemake.config["output"], "qc")
max_memory = str(snakemake.config["max_memory"]) + 'g' # kneaddata takes bytes as input
human_db = snakemake.config["human_db"]
silva_db = snakemake.config["silva_db"]
trimmomatic = snakemake.config["trimmomatic"]
sequencer_source = snakemake.config["sequencer_source"]

if snakemake.config["human_db"] != "none" and snakemake.config["silva_db"] != "none":
    subprocess.Popen(
        """
        kneaddata \
        -i %s \
        -i %s \
        -o %s \
        -t %s \
        -p %s \
        --max-memory %s \
        -db %s \
        -db %s \
        --trimmomatic-options %s \
        --sequencer-source %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --cat-final-output \
        --run-trim-repetitive \
        --decontaminate-pairs strict &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, human_db, silva_db, trimmomatic, sequencer_source), shell=True).wait()

elif snakemake.config["human_db"] == "none" and snakemake.config["silva_db"] != "none":
    subprocess.Popen(
        """
        kneaddata \
        -i %s \
        -i %s \
        -o %s \
        -t %s \
        -p %s \
        --max-memory %s \
        -db %s \
        --trimmomatic-options %s \
        --sequencer-source %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --cat-final-output \
        --run-trim-repetitive \
        --decontaminate-pairs strict &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, human_db, trimmomatic, sequencer_source), shell=True).wait()

elif snakemake.config["human_db"] != "none" and snakemake.config["silva_db"] == "none":
    subprocess.Popen(
        """
        kneaddata \
        -i %s \
        -i %s \
        -o %s \
        -t %s \
        -p %s \
        --max-memory %s \
        -db %s \
        --trimmomatic-options %s \
        --sequencer-source %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --cat-final-output \
        --run-trim-repetitive \
        --decontaminate-pairs strict &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, silva_db, trimmomatic, sequencer_source), shell=True).wait()

elif snakemake.config["human_db"] == "none" and snakemake.config["silva_db"] == "none":
    subprocess.Popen(
        """
        kneaddata \
        -i %s \
        -i %s \
        -o %s \
        -t %s \
        -p %s \
        --max-memory %s \
        --trimmomatic-options %s \
        --sequencer-source %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --cat-final-output \
        --run-trim-repetitive \
        --decontaminate-pairs strict &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, trimmomatic, sequencer_source), shell=True).wait()
