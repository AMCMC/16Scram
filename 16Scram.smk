import pandas

# Generate a sample mapping file
# Example:
# ls BARIA.Metagenome.DIABAR*1.fq | awk '{print $1"\t"$1"\t"$1}' | sed -e 's/.1.fq//1' | sed -e 's/.1.fq/.2.fq/2' | awk 'BEGIN{print "Sample\tfwfq\trvfq"}{print $0}' > mapping

samplemapping = pandas.read_csv("mapping", sep="\t")

g_slurm_partition="rome"

# get sample info from mapping file

SAMPLES = list(samplemapping["Sample"])

dictionary1 = {}
for i in range(len(SAMPLES)):
    dictionary1[SAMPLES[i]] = list(samplemapping["fwfq"])[i]

dictionary2 = {}
for i in range(len(SAMPLES)):
    dictionary2[SAMPLES[i]] = list(samplemapping["rvfq"])[i]

def fq1_from_sample(wildcards):
  return dictionary1[wildcards.sample]

def fq2_from_sample(wildcards):
  return dictionary2[wildcards.sample]

rule all:
    input:
        expand("{sample}_merged.fq", sample=SAMPLES),
        expand("{sample}_ASV.fa", sample=SAMPLES),
        expand("{sample}.cram", sample=SAMPLES)

# AVS_call 
rule merge:
    input:
        fw_reads=fq1_from_sample,
        rv_reads=fq2_from_sample
    output:
        merge="{sample}_merged.fq"
    log:
    conda: "16Scram"
    threads: 4
    resources:
        mem_mb=12000,
        runtime=60,
        slurm_partition=g_slurm_partition
    shell:
        """
        vsearch --fastq_mergepairs {input.fw_reads} \
          --reverse {input.rv_reads} \
          --fastq_maxdiffs 100 \
          --fastq_allowmergestagger \
          --fastqout {output.merge}
        """

rule filter:
    input:
        merge="{sample}_merged.fq"
    output:
        filtered="{sample}_filtered.fa"
    log:
    conda: "16Scram"
    threads: 4
    resources:
        mem_mb=12000,
        runtime=60,
        slurm_partition=g_slurm_partition
    shell:
        """
        vsearch --fastq_filter {input.merge} \
           -fastq_maxee 30 -fastaout {output.filtered} \
           --quiet --fasta_width 0
        """

rule callASV:
    input:
        filtered="{sample}_filtered.fa"
    output:
        ASV="{sample}_ASV.fa"
    log:
    conda: "16Scram"
    threads: 4
    resources:
        mem_mb=12000,
        runtime=60,
        slurm_partition=g_slurm_partition
    shell:
        """
        vsearch --derep_fulllength \
          {input.filtered} --output - --sizeout | vsearch --cluster_unoise - \
          --centroids - --minsize 4 | vsearch --uchime3_denovo - \
          --nonchimeras - --fasta_width 0 > {output.ASV}
        """

rule Cram16:
    input:
        ASV="{sample}_ASV.fa",
        fw_reads=fq1_from_sample,
        rv_reads=fq2_from_sample
    params:
        btindex=temp("{sample}.btindex")
    output:
        cram="{sample}.cram"
    log:
    conda: "16Scram"
    threads: 4
    resources:
        mem_mb=12000,
        runtime=60,
        slurm_partition=g_slurm_partition
    shell:
        """
	bowtie2-build {input.ASV} {params.btindex}
        bowtie2 -1 {input.fw_reads} -2 {input.rv_reads} -x {params.btindex} | samtools view -C -T {input.ASV} -o {output.cram}
        """

