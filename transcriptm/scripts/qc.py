import subprocess
import os
import shutil

short_reads_1 = snakemake.config["short_reads_1"]
short_reads_2 = snakemake.config["short_reads_2"]
output_dir_path = os.path.join(snakemake.config["output"], "qc")
max_memory = str(snakemake.config["max_memory"]) + 'g' # kneaddata takes bytes as input. change to GBs
human_db = snakemake.config["human_db"]
silva_db = snakemake.config["silva_db"]
other_db = snakemake.config["other_db"]
trimmomatic = snakemake.config["trimmomatic"]
sequencer_source = snakemake.config["sequencer_source"]
skip_qc = snakemake.config["skip_qc"]
trimmomatic_exec = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../bin/Trimmomatic-0.39/")


print('R1 = ' + short_reads_1)
print('R2 = ' + short_reads_2)
shutil.rmtree("qc/")


if skip_qc == True:
    subprocess.Popen(
        """
        echo '--skip-qc specified, skipping read QC steps...' &&
        mkdir qc
        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1 &&
        mkdir qc/clean_reads/R2 &&
        cp %s qc/clean_reads/R1 &&
        cp %s qc/clean_reads/R2 &&
        touch qc/done
        """ %
        (short_reads_1, short_reads_2), shell=True).wait()

elif human_db != "none" and silva_db != "none" and other_db != "none" and skip_qc != True:
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
        -db %s \
        --trimmomatic-options %s \
        --sequencer-source %s \
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*paired_1.fastq qc/clean_reads/R1 &&
        mv qc/*paired_2.fastq qc/clean_reads/R2 &&
        
        gzip qc/clean_reads/R1/* &&
        gzip qc/clean_reads/R2/* &&
        
        rm qc/*.fastq &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, human_db, silva_db, other_db, trimmomatic, sequencer_source, trimmomatic_exec), shell=True).wait()

elif human_db != "none" and silva_db != "none" and other_db == "none" and skip_qc != True:
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
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*paired_1.fastq qc/clean_reads/R1 &&
        mv qc/*paired_2.fastq qc/clean_reads/R2 &&
        
        gzip qc/clean_reads/R1/* &&
        gzip qc/clean_reads/R2/* &&
        
        rm qc/*.fastq &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, human_db, silva_db, trimmomatic, sequencer_source, trimmomatic_exec), shell=True).wait()

elif human_db == "none" and silva_db != "none" and other_db != "none" and skip_qc != True:
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
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*paired_1.fastq qc/clean_reads/R1 &&
        mv qc/*paired_2.fastq qc/clean_reads/R2 &&
               
        gzip qc/clean_reads/R1/* &&
        gzip qc/clean_reads/R2/* &&
        
        rm qc/*.fastq &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, silva_db, other_db, trimmomatic, sequencer_source, trimmomatic_exec), shell=True).wait()

elif human_db != "none" and silva_db == "none" and other_db != "none" and skip_qc != True:
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
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*paired_1.fastq qc/clean_reads/R1 &&
        mv qc/*paired_2.fastq qc/clean_reads/R2 &&
                
        gzip qc/clean_reads/R1/* &&
        gzip qc/clean_reads/R2/* &&
        
        rm qc/*.fastq &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, human_db, other_db, trimmomatic, sequencer_source, trimmomatic_exec), shell=True).wait()


elif human_db == "none" and silva_db != "none" and other_db == "none" and skip_qc != True:
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
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*paired_1.fastq qc/clean_reads/R1 &&
        mv qc/*paired_2.fastq qc/clean_reads/R2 &&
               
        gzip qc/clean_reads/R1/* &&
        gzip qc/clean_reads/R2/* &&
        
        rm qc/*.fastq &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, silva_db, trimmomatic, sequencer_source, trimmomatic_exec), shell=True).wait()


elif human_db != "none" and silva_db == "none" and other_db == "none" and skip_qc != True:
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
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*paired_1.fastq qc/clean_reads/R1 &&
        mv qc/*paired_2.fastq qc/clean_reads/R2 &&
                
        gzip qc/clean_reads/R1/* &&
        gzip qc/clean_reads/R2/* &&
        
        rm qc/*.fastq &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, human_db, trimmomatic, sequencer_source, trimmomatic_exec), shell=True).wait()


elif human_db == "none" and silva_db == "none" and other_db != "none" and skip_qc != True:
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
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --remove-intermediate-output \
        --bypass-trf \
		--run-trim-repetitive \
        --decontaminate-pairs strict &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*paired_1.fastq qc/clean_reads/R1 &&
        mv qc/*paired_2.fastq qc/clean_reads/R2 &&
                
        gzip qc/clean_reads/R1/* &&
        gzip qc/clean_reads/R2/* &&
        
        rm qc/*.fastq &&

        touch qc/done
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, other_db, trimmomatic, sequencer_source, trimmomatic_exec), shell=True).wait()


elif human_db == "none" and silva_db == "none" and other_db == "none" and skip_qc != True:
    #kneaddata doesn't reorder if not using bowtie2 decontamination
    subprocess.Popen(
        """
        kneaddata \
        -i %s \
        -i %s \
        -o %s \
        -p %s \
        -t %s \
        --max-memory %s \
        --sequencer-source %s \
        --trimmomatic %s \
        --run-fastqc-end \
        --reorder \
        --bypass-trf \
		--run-trim-repetitive \
        --remove-intermediate-output &&

        mkdir qc/clean_reads &&
        mkdir qc/clean_reads/R1/ &&
        mkdir qc/clean_reads/R2/ &&

        mv qc/*trimmed.1.fastq qc/clean_reads/R1 &&
        mv qc/*trimmed.2.fastq qc/clean_reads/R2 &&

        rm qc/*.fastq
        """ %
        (short_reads_1, short_reads_2, output_dir_path, snakemake.threads, snakemake.threads, max_memory, sequencer_source, trimmomatic_exec), shell=True).wait()

    subprocess.Popen(
        "echo 'Correcting read order...' && "

        "mv qc/clean_reads/R1/*.fastq qc/clean_reads/R1/unordered.fastq && "
        "mv qc/clean_reads/R2/*.fastq qc/clean_reads/R2/unordered.fastq && "

        "repair.sh in1=qc/clean_reads/R1/unordered.fastq in2=qc/clean_reads/R2/unordered.fastq out1=qc/clean_reads/R1/reordered_R1.fastq out2=qc/clean_reads/R2/reordered_R2.fastq && "
        
        "gzip qc/clean_reads/R1/reordered_R1.fastq && "
        "gzip qc/clean_reads/R2/reordered_R2.fastq && "

        "rm qc/clean_reads/R1/unordered.fastq && "
        "rm qc/clean_reads/R2/unordered.fastq && "

        "touch qc/done", shell=True).wait()
