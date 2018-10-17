#!/usr/bin/env python 
import re # to use regular expressions
import argparse # To allow users setting flags in order to define what to filter for.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser(description="Extracts sequence windows from a fasta file based on sequence IDs and positions and writes them to a new fasta file")
parser.add_argument("-f","--fasta", help="fasta file to extract from")
parser.add_argument("-t","--text", help="Tab-delimited text file containing the sequence-id in the first column and the
                        position of the targeted window-centers in the second column")
parser.add_argument("-w","--windowsize", help="sum of the number of basepairs preceding and following the position that is specified in the text file. A window of 10 kb means that 5 kb shall precede and 5 kb shall follow the speciied position.")
parser.add_argument("-p","--positionspecifier", help="string that specified the position which is provided in the second column of the text file. This could be for example 'SNP' or 'windowcenter' and will appear as comment following the sequence IDs in the output fasta files")
args=parser.parse_args()


def WindowsFromFasta(fasta,text,windowsize,positionspecifier):
        '''
        Extracts sequence windows from a fasta file based on sequence IDs and positions and writes them to a new fasta file
        '''

        windowsize = int(windowsize)
        
        iterator = SeqIO.parse(fasta, "fasta") # Open the fasta file 
        fastadictionary = SeqIO.to_dict(iterator) # Create a dictionary from the fasta file
        
        # open the text file with the sequence IDs and the start and end positions
        textfile = open(text,'U')
        unique = {} # serves to connect Sequence IDs to fasta comments 
    
        for line in textfile:
            line= line.strip().split('\t') # read from a tab-delimited text file
            # a combination of the identifier, and the sequence position (the center of the window that shall be extracted).
        
            p = int(line[1])
            s = int(line[1])-windowsize/2# If this is negative, I want to set it to 0
            if s<1: s=1
            e = int(line[1])+windowsize/2 # An index that is too large doesn't give an error if it is used in slicing. So, if the string has no base at position 20000, it will print out the last base before that. However, I want to inform then what the last base is
            if e > int(len(fastadictionary[line[0]].seq)): e=int(len(fastadictionary[line[0]].seq))

            outname = line[0]

            # Remove duplicate windows but keep all positionspecifiers
      
            u = str(outname + str(s) + '-' + str(e))
            if u not in unique:
                unique[u] =  'window: ' + str(s) + '-' + str(e) + ', ' + positionspecifier + ' at ' + str(line[1]) # Create a comment for this record      
            else:
                existingcomment = unique[u]
                unique[u] =  existingcomment +  ', ' +positionspecifier + ' at ' + str(line[1]) # update comment for this record

        textfile.close()
        
        textfile = open(text,'U')
        unique2=[] # serves to check for duplicate windows

        for l in textfile:
            l= l.strip().split('\t') # read from a tab-delimited text file
            # a combination of the identifier, and the sequence position (the center of the window that shall be extracted).
      
            p = int(l[1])
            s = int(l[1])-windowsize/2# If this is negative, I want to set it to 1
            if s<1: s=1
            e = int(l[1])+windowsize/2 # An index that is too large doesn't give an error if it is used in slicing. So, if the string has no base at position 20000, it will print out the last base before that. However, I want to inform then what the last base is
            if e > int(len(fastadictionary[l[0]].seq)): e=int(len(fastadictionary[l[0]].seq))

            outname = l[0]

            # Print windows only once but with all positionspecifiers
            u = str(outname + str(s) + '-' + str(e))
            if u not in unique2:
                unique2.append(u)
                sequence_record=SeqRecord(fastadictionary[l[0]].seq[s-1:e],id = outname,description = unique[u])
                print sequence_record.format("fasta")

        

            
WindowsFromFasta(str(args.fasta),str(args.text),int(args.windowsize),str(args.positionspecifier))   