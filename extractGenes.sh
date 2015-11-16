#!/bin/bash
# Author: Ioannis Moustakas, i.moustakas@uva.nl
# Title: Build sorted bam files containing only the genes of interest aiming to visualize on a genome browser. 
#	 Also build the corresponding references fasta file(s)

# Path where the sam files are stored
samPath=/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Scratch/pipeline

# output path
OutputPath=/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Results/visualizeGenes/dspAvr2
if [ ! -d "$OutputPath" ]; then
  mkdir $OutputPath
fi

# Transcriptome (reference) file
transcriptomeFile=/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1208-Frank_Takken/MAD1208-P001-DTL_Hotel/MAD1208-P001-E001_2014_RNASeq_Tomato_svleeuw1/Scratch/pipeline/reference/Solanum_lycopersicum_I2_Avr2.fa

# Genes to extract, as mentioned in reference file
genes=(dspAvr2)

# List of sam file sample names to grep the genes. Only the unique part of the file name
samples=(S01 S02 S03 S16 S17 S18 S29 S30 S31 S07 S08)


for gene in ${genes[@]}; 
do  
  # Get the gene sequence form the trancriptome file
  grep -A 1 $gene $transcriptomeFile | fold -w 120 > $OutputPath/$gene"_Sequence.fa"

  # go through all sam files and get all alignments of the gene of interest
  # then sam => indexed bam file
  for sample in ${samples[@]}; 
  do
    echo -e "Working on sample $sample \n"
    grep $gene $samPath/$sample"_tmap_transc.sam" > $OutputPath/$sample"_"$gene".sam"
    samtools view -Sb $OutputPath/$sample"_"$gene".sam" > $OutputPath/$sample"_"$gene".bam"
    samtools sort $OutputPath/$sample"_"$gene".bam" $OutputPath/$sample"_"$gene".sorted"
    samtools index $OutputPath/$sample"_"$gene".sorted.bam"
    rm $OutputPath/$sample"_"$gene".bam"
    numAlign=$(grep -c $gene $OutputPath/$sample"_"$gene".sam")
    echo -e "Reads in sample $sample aligned on gene $gene: $(($numAlign-1))"
  done
done
