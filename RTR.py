#!/usr/bin/env python 
import re # to use regular expressions
import argparse # To allow users setting flags in order to define what to filter for.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Selecting MethylRAD tags with less degenerate adaptor ends")
parser.add_argument("-f","--fasta", help="fasta file with MethylRAD tags (output from the InSilicoTypeIIbDigestion.py script). The tag endings must be the complementary bases to adaptor 1 or 2.")
parser.add_argument("-b1","--base1", help="specific base of adaptor 1 that is complementary to the cutting site overhang. The script converts it to the complementary base on the 5 prime end of the MethylRAD tag")
parser.add_argument("-b2","--base2", help="specific base of adaptor 2 that is complementary to the cutting site overhang. The script converts it to the complementary base on the 5 prime end of the MethylRAD tag")
parser.add_argument("-l1","--location1", help="location of specific base 1 in the cutting site overhang. Location=1 refers to the first outermost base on either site of the MethylRAD tag")
parser.add_argument("-l2","--location2", help="location (in bp) of specific base 2 in the cutting site overhang. Location=1 refers to the first outermost base on either site of the MethylRAD tag")

args=parser.parse_args()


def RTR(fasta,base1,base2,location1,location2):
        '''
        Selecting MethylRAD tags with degenerate adaptor ends
        '''
        complementbases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        
        b1 = str(base1) 
        b2 = str(base2) 
        c1 = complementbases[str(base1)]
        c2 = complementbases[str(base2)]
        l1 = int(location1)
        l2 = int(location2)

#        iterator = SeqIO.parse(fasta, "fasta") # Open the fasta file 
#        fastadictionary = SeqIO.to_dict(iterator) # Create a dictionary from the fasta file

        tagcount = {} # Will count the number of tags that fit to these adaptor endings

        recordnumber=0
        hitnumber=0

        fastafilemarked=open(str(fasta.replace(".fasta","_RTRMarked.fasta")),"w")
        fastafileextracted=open(str(fasta.replace(".fasta","_RTRExtracted.fasta")),"w")
        
        for record in SeqIO.parse(fasta, "fasta"):
            recordnumber += 1
            # I need to take here the reverse complement
            if c1 == str(record.seq)[l1-1]:
                if b2 == str(record.seq)[len(str(record.seq))-l2]:
                    #print(str(record.id))
                    #print(str(record.seq))
                    hitnumber +=1
                    currentsequence=SeqRecord(record.seq,id=record.id,description="RTR tag")
                    SeqIO.write(currentsequence,fastafilemarked,"fasta")
                    SeqIO.write(currentsequence,fastafileextracted,"fasta")
            elif c2 == str(record.seq)[l1-1]:
                if b1 == str(record.seq)[len(str(record.seq))-l2]:
                    #print(str(record.id))
                    #print(str(record.seq))
                    hitnumber +=1
                    currentsequence=SeqRecord(record.seq,id=record.id,description="RTR tag")
                    SeqIO.write(currentsequence,fastafilemarked,"fasta")
                    SeqIO.write(currentsequence,fastafileextracted,"fasta")
            else:
                currentsequence=SeqRecord(record.seq,id=record.id,description="Not included")
                SeqIO.write(currentsequence,fastafilemarked,"fasta")


        percent = 100*float(hitnumber)/float(recordnumber)

        
        fastafilemarked.close()
        fastafileextracted.close()

        print ( str(hitnumber) + " of " + str(recordnumber) + " (" + str(percent) + "%) tags fit to these adaptor endings")

                
RTR(str(args.fasta),str(args.base1),str(args.base2),int(args.location1),int(args.location1))
