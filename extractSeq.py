# Author: Ioannis Moustakas
# Title: extract sequences from fasta using GFF annotation
# 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import re
import click
 

# read the fasta file and save it in a dictionary
# fasta header => sequence
@click.group()
def cli():
    pass 

@cli.command()
@click.option('--gff_file',
              type=str,
              help='GFF annotation file')
@click.option('--fasta_file',
              type=str,
              help='Fasta file')
@click.option('--output_fasta_file',
              type=str,
              help='Output fasta file where the sequences are saved')
def extract_seq(fasta_file, gff_file, output_fasta_file):          
    fastaHeadToSeq = SeqIO.to_dict(SeqIO.parse(open(fasta_file),'fasta'))
    
    listOfSeq = []
    outHandle = output_fasta_file
    
    in_handle = open(gff_file)
    # skip first line
    next(in_handle)
    reader = csv.reader(in_handle, delimiter='\t')
    for line in reader:
        seqID = line[0]
        typeOfSeq = line[2]
        start = int(line[3])
        end = int(line[4])
        strand = line[6]
        attributes = line[8]
        listOfAttr=attributes.split(sep=';')   
        dict={}
        for element in listOfAttr:
            elements = element.split(sep='=')
            if len(elements)>1:
                dict[elements[0]]=elements[1]
        # key used here is from the attributes field of GFF file
        uniqueID=dict["ID"]
        annotation=dict["Annotation"]
        fastaHeader=uniqueID+"_"+annotation
        fastaHeader=re.sub(' ', '_', fastaHeader)    
    
        # get the sequence using the seqID (fasta header in the original file)
        sequence = str(fastaHeadToSeq[seqID].seq)
        subSequence = sequence[start:end]
        record = SeqRecord(Seq(subSequence),
                       id=uniqueID, description=annotation) 
        # if the sequence is on the minus strand then reverse complement
        if strand=="-":
            record.reverse_complement
            
        listOfSeq.append(record)
     
    SeqIO.write(listOfSeq, outHandle, "fasta")

if __name__ == '__main__':
    cli()