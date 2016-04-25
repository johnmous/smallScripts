# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:55:47 2015

@author: Ioannis Moustakas, i.moustakas@uva.nl
Title: From a bam file aligned with Tmap, Softclipping ON, extract 
aligned and softclipped fractions of the read sequence and put them seperately in a fasta file. 
Read ID in the fasta header
"""

import pysam

BAMFile = "unmapped_ArabidopsisPolyA_On_ATtransriptNorRNA_SoftclON.bam"
#"unmapped_tomatoPolyA_OnTranscriptomeRNASoftClOff_ON_TranscriptomeSoftClOn.bam"
#"/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1207-Puratos/MAD1207-P001-Bacillus_Subtilis/MAD1207-P001-E001_2014_RNASeq_25x_svleeuw1/Scratch/alignOnGenome/unmapedByBt2TmapSoftOn.sorted.bam"

samfile = pysam.AlignmentFile(BAMFile, "rb")
# open output fasta file
fastaFile = open("softClippedSeq.fa", 'w')

# Count total reads, number of softclipped reads, number of sequences extracted
readCount = 0
mappedReadsCount = 0 
softclReadsCount = 0
seqExtracted = 0
singleSoftclReads = 0 
twosidedSoftclReads = 0 
alignedSeqExceeding25 = 0 

# for each read in the Sam/Bam file
for read in samfile:
    readCount += 1
    cigarTuples = read.cigartuples
    readSequence = read.query_sequence
    readLength = read.query_length
    # Check if cigar string is emtpy (non aligned read)
    if cigarTuples is not None:
        countSoftclSides = 0
        mappedReadsCount += 1
        fastaHeader = '>'+read.query_name+'\n'
        #if read softclipped on the left, extact left softcliped part
        # Keep only the sequences that are larger than 25 nt long
        if (cigarTuples[0][0]==4) and (cigarTuples[0][1] > 25):
            countSoftclSides += 1 
            seqExtracted += 1 
            fastaFile.write(fastaHeader)
            fastaFile.write(readSequence[ :cigarTuples[0][1]]+'\n')
        
        #if read softclipped on the rigth, extact right softcliped part
        if (cigarTuples[-1][0]==4) and (cigarTuples[-1][1] > 25):
            countSoftclSides += 1 
            seqExtracted += 1 
            fastaFile.write(fastaHeader)
            fastaFile.write(readSequence[readLength - cigarTuples[-1][1]: ]+"\n")
                
        # if softclipped, right or left, extract the aligned sequence as well
        if (cigarTuples[0][0]==4) or (cigarTuples[-1][0]==4):
            alignedReadSeq = read.query_alignment_sequence
            softclReadsCount += 1            
            if len(alignedReadSeq) > 25:
                seqExtracted += 1
                alignedSeqExceeding25 += 1
                fastaFile.write(fastaHeader)
                fastaFile.write(alignedReadSeq+"\n")
        if (countSoftclSides == 1):
            singleSoftclReads += 1
        elif (countSoftclSides == 2):
            twosidedSoftclReads += 1
            
                
print("Number of reads: ", readCount)
print("Number of reads aligned on reference: ", mappedReadsCount)
print("Number of softclipped reads: ", softclReadsCount)
print("Number of softclipped sequences extracted: ", seqExtracted)
print("Number of single side softclipped reads (softcl sequence > 25): ", singleSoftclReads)
print("Number of two sided softclipped reads (softcl sequence > 25): ", twosidedSoftclReads)
print("Number of aligned sequences exceeding length threshold (25): ", alignedSeqExceeding25)
                
fastaFile.close()
samfile.close()


#==============================================================================
#             print("right")        
#             print(read.cigarstring)
#             print(read.query_alignment_sequence)
#             print(read.query_length)
#             print(readSequence)
#             print(readSequence[ -cigarTuples[-1][1]: ])
#             
#             print("left")
#             print(read.cigarstring)
#             print(read.query_alignment_sequence)
#             print(readSequence[cigarTuples[0][1]: ])
#             print(readSequence[ :cigarTuples[0][1]])
#             print(read.query_length)
#             print(">" + read.query_name)
#==============================================================================
