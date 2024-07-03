#ls *R1.fastq.gz  | awk '{print $1"\t"$1"\t"$1}' | sed -e 's/.R1.fastq.gz//1' | sed -e 's/.R1.fastq.gz/.R2.fastq.gz/2' | awk 'BEGIN{print "Sample\tfwfq\trvfq"}{print $1"\t/media/3T_sdb/16S_datasets/cramtest/fastqfiles"$2"\t/media/3T_sdb/16S_datasets/cramtest/fastqfiles/"$3}' > ../mapping

snakemake -s 16Scram.smk --use-conda --cores 16
