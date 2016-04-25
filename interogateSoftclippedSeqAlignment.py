# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:44:44 2015

@author: Ioannis Moustakas, i.moustakas@uva.nl

@Title: A script to trace the genes sequence fragments map on the arabidopsis tracriptome. 
Sequence fragments derived from the same read have the same ID
"""

import pysam

# open the bam file to interogate
BAMFile = "softClippedSeqOnATtranscriptTAIR10_No_rRNA.SoftclOFF.sorted.bam"
samfile = pysam.AlignmentFile(BAMFile, "rb")

# Make a dictionary readID => GeneID out of the BAM file and save it
queryName2GeneName = {}

# Go through all entries in the BAM file
for read in samfile:
    # Skip unmapped reads    
    if read.is_unmapped == False:
        readID = read.query_name
        referenceID = samfile.getrname(read.tid)
        # Skip reads that map on nothing (redundant, but I keep it in)
        if referenceID != None:
            # if readID does not exist as key, create readID => [referenceID]
            if readID not in queryName2GeneName: 
                queryName2GeneName[readID] = [referenceID]
            # else append to the already exitsting value [list] of that key                 
            else:
                queryName2GeneName[readID].append(referenceID)
            
readInSameGene = 0            
readInSameChrom= 0
readInDifferentChrom = 0

# iterate thtough the dictionary and calculate some statistics
for x in queryName2GeneName:
    if len(queryName2GeneName[x])==2:
        gene1 = queryName2GeneName[x][0][:-2]
        gene2 = queryName2GeneName[x][1][:-2]
        if gene1 == gene2:
            readInSameGene += 1
        elif gene1[2] == gene2[2]:
            readInSameChrom += 1
        else:
            readInDifferentChrom += 1
    elif len(queryName2GeneName[x]) == 3:
        gene1 = queryName2GeneName[x][0]
        gene2 = queryName2GeneName[x][1]
        gene3 = queryName2GeneName[x][2]
        if gene1 == gene2 == gene3:
            readInSameGene += 1
        elif gene1[2] == gene2[2] == gene3[2]:
            readInSameChrom += 1
        else:
            readInDifferentChrom += 1
        
print("Number of reads with fragments in same gene: ", readInSameGene)
print("Number of reads with fragments in same chromosome: ", readInSameChrom)
print("Number of reads with fragments in different chromosome: ", readInDifferentChrom)
