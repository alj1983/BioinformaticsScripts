#!/usr/bin/env python 


import re # to use regular expressions
import argparse # To allow users setting flags in order to define what to filter for.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

WindowIdentificationPattern = re.compile('.*window: (\d+)-(\d+).*') # Extract the window start and end positions from the fasta descriptions
parser = argparse.ArgumentParser(description="Merges windows that were extracted from the same contig, removes duplicate sequences, and merges the sequence-ids of duplicate sequences")
parser.add_argument("-f","--fasta", help="fasta files to merge",nargs='+')
parser.add_argument("-r","--referencefasta", help="fasta file of the genome/transcriptome that was used to extract the windows in the single fasta files")
args=parser.parse_args()


def MergeSequenceWindows(fasta,referencefasta):
        '''
        Merges windows that were extracted from the same contig, removes duplicate sequences, and merges the sequence-ids of duplicate sequences
        '''

        openreferencefasta = SeqIO.parse(referencefasta, "fasta") # Open the fasta file 
        referencefastadictionary = SeqIO.to_dict(openreferencefasta) # Create a dictionary from the fasta file
        
        uniquewindows = {} # serves to connect contigs to different extracted windows
        uniquedescriptions = {} # serves to connect contigs to different descriptions (location of SNPs or windowcenters)
        for f in fasta:
            iterator = SeqIO.parse(f, "fasta") # Open the fasta file 
            #fastadictionary = SeqIO.to_dict(iterator) # Create a dictionary from the fasta file
            records = list(iterator)
            for sequence in records:
                #print sequence.description
                windowidentification = WindowIdentificationPattern.search(str(sequence.description))
                actualwindowstart = windowidentification.group(1)
                actualwindowend = windowidentification.group(2)

                # Get the single SNPs, or windowcenters in a list (will be added to the sequence description)
                l= str(sequence.description).strip().split(',')
                actualdescriptions = list(l[1:len(l)])
                if sequence.id not in uniquewindows:
                    uniquewindows[sequence.id] = [actualwindowstart,actualwindowend] # This saves the start and end of the window
                    uniquedescriptions[sequence.id] = actualdescriptions
                else:
                    existingwindowstart = uniquewindows[sequence.id][0]
                    existingwindowend = uniquewindows[sequence.id][1]
                    minimumwindowstart = min([int(existingwindowstart),int(actualwindowstart)])
                    maximumwindowend = max([int(existingwindowend),int(actualwindowend)])
                    uniquewindows[sequence.id] = [int(minimumwindowstart),int(maximumwindowend)]
                    
                    existingdescriptions = uniquedescriptions[sequence.id]
                    existingdescriptions.extend(actualdescriptions)
                    uniquedescriptions[sequence.id] = existingdescriptions


        for uniqueid in uniquewindows:
            actualsequence = referencefastadictionary[uniqueid].seq
            actualwindowstart = int(uniquewindows[uniqueid][0])
            actualwindowend = int(uniquewindows[uniqueid][1])
            actualsequencewindow=actualsequence[actualwindowstart-1:actualwindowend-1]
            windowlength = actualwindowend-actualwindowstart+1
            actualdescriptions = str('window with ' + str(windowlength) + ' bp: ' + str(actualwindowstart) + '-' + str(actualwindowend) + ', ' + ','.join(uniquedescriptions[uniqueid]))
            sequence_record=SeqRecord(actualsequencewindow,id = str(uniqueid),description = actualdescriptions)                
            print sequence_record.format("fasta")               
        
MergeSequenceWindows(args.fasta,str(args.referencefasta))   