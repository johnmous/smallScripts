#!/bin/bash
# Title: A script to calculate the avarage coverage per bam file
# Author: Ioannis Moustakas, i.moustakas@uva.nl

# sufix of bam files
sufix=sorted.bam

# list files in current directory
bamFiles=$(ls *$sufix)

echo -e "Sample\tAverageCoverage"

for filename in ${bamFiles[@]};
do
  prefix=$(echo $filename | tr "$sufix" "\n")
  prefix=$prefix
  averageCov=$(genomeCoverageBed -ibam $filename -d | awk '{sum+=$3}END{print sum/NR}')
  echo -e "$prefix\t$averageCov"
done