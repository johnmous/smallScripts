#!/bin/bash
##################################################################
### 2015-05-08 script to generate processed data for Vipergen
### Author: Oskar Bruning: o.bruning@uva.nl
### 2015-05-22 made a full pipeline script 
### 2015-11 Optimized by: Ioannis Moustakas, i.moustakas@uva.nl
### usage: include the sff file locations and a working directory
##################################################################


#### start to change before each run #######
# select the runs to use and output file prefix (Vippergen sample ID) of the final zip archive
declare -a RUNID=("RID0341" "RID0342")
prefix=VIPS135b

# fixate working folder
base_dir=/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1211-Vipergen/MAD1211-P001-Leif_Kongskov_Larsen/MAD1211-P001-E034_2016_VIPS135b_svleeuw1/
#### end to change before each run #######

# go to there
cd $base_dir

# create temp working folder
mkdir $base_dir/Scratch/temp

# Run location folder
Runs_dir=/zfs/datastore0/group_root/MAD-RBAB/04_MAD-RBAB-runs/data

# Download folder 
dl_dir=/zfs/datastore0/system/web/download/MAD1211

# loop for performing the processing steps
for i in "${RUNID[@]}"
do
  # create fna files from BAM files
  samtools view $Runs_dir/$i/rawlib.basecaller.bam \
  | awk '{if (length($10)>175) printf(">%s\n%s\n", $1, $10) }' > $base_dir/Scratch/temp/$i"_over175.fna"
done

# run in R
Rscript $base_dir/Scripts/fna_split.R

# make zipfile and cp to download location
# get unique id based on the runids
ZIPID=`echo ${RUNID[@]} | tr ' ' '_'`

cd $base_dir/Scratch/temp

zip -9 --junk-paths $base_dir/Scratch/$prefix"_"$ZIPID".zip" *start*
cp $base_dir/Scratch/*.zip $dl_dir

filesInZip=$(unzip -l $dl_dir/$prefix"_"$ZIPID".zip")

mailBody="Vipergen Script is finished and the output is saved in the follwing zip file that contains: \n $filesInZip"
# clean-up
#rm -r $base_dir/Scratch/temp/

# Script mails user after finishing up
echo -e "$mailBody" | mailx -s VipergenDone -c S.M.vanLeeuwen@uva.nl i.moustakas@uva.nl 

