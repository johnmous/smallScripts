# setup workspace
X11.options(type="Xlib")
options(stringsAsFactors = FALSE)
library(SeqLibR)

# set number of reads per file
nr_reads <- 1000000

setwd("./Scratch/temp")
# read in fasta files
for(R in grep("over175.fna", dir(), value="T")){
 temp_run <- ns.read.fasta(R)

 # create fasta files with numbers of reads as specified above
 for (i in seq(1, length(temp_run), nr_reads)){
  if(length(i:length(temp_run)) > nr_reads){
   temp_sub <- temp_run[i:((i+nr_reads)-1)]
   ns.write.fasta(paste(sub(".fna", "",  R),"_start", i, "_stop", ((i+nr_reads)-1), ".fna", sep=""), temp_sub)
  }
  if(length(i:length(temp_run)) < nr_reads){
   temp_sub <- temp_run[i:length(temp_run)]
   ns.write.fasta(paste(sub(".fna", "",  R),"_start", i, "_stop", length(temp_run), ".fna", sep=""), temp_sub)
  }
 }
}

