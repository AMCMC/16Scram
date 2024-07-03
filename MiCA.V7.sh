#!/usr/bin/bash

#### input ####

mappingfile=$1 #format sample.ID fq1 fq2
outputdir=$2 # provide full output path
subsample=$3 # maximum number of reads to use for analyses
ee=$4 # expected error filter

#### pipeline parameters ####

filterpurity=1 # DEFAULT:1 is yes 0 is false Not implemented
minprev=0.002  # Used for global reference ASV
minmaxabun=0.001  # Used for global reference ASV
nsamples=`(cat $mappingfile | awk '{print $1}' | sort | uniq -c | wc -l)` # check duplicates?

#### create working environment & directory ####

source /home/mdavids/miniconda3/etc/profile.d/conda.sh

conda init bash
conda activate mica

mkdir $outputdir
mkdir $outputdir/merged
mkdir $outputdir/filtered
mkdir $outputdir/ASV
mkdir $outputdir/abundance
mkdir $outputdir/local_reference
mkdir $outputdir/Data/

#### Functions wrappers ####

merge_fq () {
vsearch --fastq_mergepairs $2 --reverse $3 --fastq_maxdiffs 100 --fastq_allowmergestagger --fastqout $4/merged/$1.merged.fq
}
export -f merge_fq

filter_and_sub () {
subsample=`(echo $2 | awk '{print $0*2}')` # assumes 2 line fasta
sample=$1
outdir=$4
maxee=$5
threads=10
vsearch --fastq_filter $outdir/merged/$sample.merged.fq -fastq_maxee $maxee -fastaout - --quiet --threads 1 --fasta_width 0 | awk '{if (substr($1,1,1)==">"){x=1}}{if (substr($2,1,1)=="N"){x=0}}{if (x==1){print $0}}'  | head -n $subsample > $outdir/filtered/$sample.filtered.fa
}
export -f filter_and_sub

call_ASV () {
sample=$1
outdir=$4
minsize=4
cat $outdir/filtered/$sample.filtered.fa | vsearch --derep_fulllength - --output - --sizeout | vsearch --cluster_unoise - --centroids - --minsize $minsize | vsearch --uchime3_denovo - --nonchimeras - --fasta_width 0 | awk -F"size=" '{if (NF==2){size=$2}else{print ">"$0"\t"size"\n"$0}}' > $outdir/ASV/$sample.ASV.fasta
}
export -f call_ASV

local_abundance () {
sample=$1
outdir=$4
cat $outdir/filtered/$sample.filtered.fa  | vsearch --usearch_global - --db $outdir/ASV/$sample.ASV.fasta --otutabout $outdir/abundance_local/$sample.abundance --id 0.97
}
export -f local_abundance

global_abundance () {
sample=$1
outdir=$4
cat $outdir/filtered/$sample.filtered.fa  | vsearch --usearch_global - --db $outdir/reference.ASVs.fasta --otutabout $outdir/abundance_global/$sample.abundance --id 0.97
}
export -f global_abundance

priors_abundance () {
sample=$1
outdir=$4
awk '{if (ARGIND==1 && NR%2==0){arr[$1]++}}{if (ARGIND==2 && NR%2==1){header=$0}}{if (ARGIND==2 && NR%2==0){if ($1 in arr){print arr[$1]"\t"header"\t"$1}}}' $outdir/filtered/$sample.filtered.fa $outdir/Data/reference.ASVs.fasta | sort -nr | awk '{print $2"\n"$3}' > $outdir/local_reference/$sample.ref.fa
cat $outdir/filtered/$sample.filtered.fa  | vsearch --usearch_global - --db $outdir/local_reference/$sample.ref.fa --otutabout $outdir/abundance/$sample.abundance --id 0.97
}
export -f priors_abundance

#### Pipeline ####

#1 Merge paired-end reads
cat $mappingfile | awk '{print $0"\t""'"$outputdir"'"}' | xargs -l bash -c 'merge_fq $0 $1 $2 $3'

#2 Filter and subsample 
cat $mappingfile | awk '{print $1"\t""'"$subsample"'""\t""'"$filterpurity"'""\t""'"$outputdir"'""\t""'"$ee"'"}' | xargs -l bash -c 'filter_and_sub $0 $1 $2 $3 $4'

#3 Call ASVs
cat $mappingfile | awk '{print $0"\t""'"$outputdir"'"}' | xargs -l bash -c 'call_ASV $0 $1 $2 $3'

#4 Determine Global ASV reference set
cat $mappingfile | awk '{print "'"$outputdir"'""/ASV/"$1".ASV.fasta"}' | xargs cat | grep ">" | sed -e 's/>//g' | awk '{arr[$1]++}{if ($2>arr2[$1]){arr2[$1]=$2}}END{for (i in arr){print i"\t"arr[i]"\t"arr2[i]}}' | sort -nr -k2,2 -k3,3  | awk ' $3/"'"$subsample"'" > "'"$minmaxabun"'" || $2/"'"$nsamples"'" > "'"$minprev"'" ' | awk '{print ">ASVX_"NR"\n"$1}' > $outputdir/Data/reference.ASVs.fasta

#5 Priors ASV abundance
cat $mappingfile | awk '{print $0"\t""'"$outputdir"'"}' | xargs -l bash -c 'priors_abundance $0 $1 $2 $3'

#6 Generate tree

mafft --auto $outputdir/Data/reference.ASVs.fasta > $outputdir/Data/reference.ASVs.aln.fasta
FastTree -nt -gtr -gamma $outputdir/Data/reference.ASVs.aln.fasta > $outputdir/Data/reference.ASVs.tree
