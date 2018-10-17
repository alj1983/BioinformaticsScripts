#!/usr/bin/env python 
import re # to use regular expressions
import argparse # To allow users setting flags in order to define what to filter for.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Counts recognition sites of restriction endonucleases and extracts the fragments resulting from an in silico digest with type IIb enzymes")
parser.add_argument("-f","--fasta", help="fasta file to extract from")
parser.add_argument("-r","--recognition", help="Comma separated list of recognition sites")
parser.add_argument("-d5","--distance5prime", help="Comma separated list of distances (in bp, for each recognition site) from the first base of the recognition site to the cutting site in 5 prime direction")
parser.add_argument("-l","--length", help="Comma separated list of total length (in bp) to cut out for each of the recognition sites")
args=parser.parse_args()


def InSilicoDigest(fasta,recognition,distance5prime,length):
        '''
        Counts recognition sites of IIb-restriction enzymes and extracts the fragments resulting from an in silico digest
        '''

        
        distance5prime = str(distance5prime).split(',')
        length = str(length).split(',')
        recognitionlist=str(recognition).split(',') # Split the lines by comma separators
        
        iterator = SeqIO.parse(fasta, "fasta") # Open the fasta file 
        fastadictionary = SeqIO.to_dict(iterator) # Create a dictionary from the fasta file
        
        recognitionsitenumbers = {} # Will contain the counts of occurrence for each recognition site

        basesinfasta = int(0) # count total number of bases in fasta file
        for record in SeqIO.parse(fasta, "fasta"):
                basesinfasta=int(basesinfasta)+int(len(record.seq))

        print(str("Total bases in fasta file: " + str(basesinfasta)))
                
        

        recs = [] # Saves the checked recognition sites for each hit
        starts = [] #saves the start positions of recognition site hits
        ends = [] #saves the end positions of recognition site hits
        allids = [] # saves ids for all recognition sites

        n = 0
        for r in recognitionlist: # Going through all the recognition sites
            regexp = re.compile(str(r),re.IGNORECASE)
            ids = [] # saves the ids with recognition site hits
            for record in SeqIO.parse(fasta, "fasta"):
                rhit = regexp.finditer(str(record.seq))
                for s in rhit:
                    ids.append(record.id)
                    allids.append(record.id)
                    starts.append(s.start()-int(distance5prime[n]))
                    ends.append(s.start()-int(distance5prime[n])+int(length[n]))
                    recs.append(str(r))

            recognitionsitenumbers[str(r)]=len(ids)
            n += 1

        for r in recognitionsitenumbers:
                print str(str(r) + " recognition sites: " + str(recognitionsitenumbers[r]) + ", density in fasta file: " + str(int(basesinfasta)/int(recognitionsitenumbers[r])))

        
        fastafile=open(str(fasta.replace(".fasta","_digested.fasta")),"w")
        # Now extract the recognition sites


        
        for x in range(0,len(allids)):

            startofrange = int(starts[x])
            endofrange = int(ends[x])
            
            # if the start is negative, I want to set it to 1
            if startofrange<1: startofrange=1
                
            # if the end is larger than the actual contig length, set it to the end of the actual contig length
            if endofrange > int(len(fastadictionary[allids[x]].seq)): endofrange=int(len(fastadictionary[allids[x]].seq))
            
            currentidentifier = str(str(allids[x]) + ',' + str(recs[x]) + ',' + str(startofrange) + '-' + str(endofrange))
            currentsequence=SeqRecord(fastadictionary[allids[x]].seq[startofrange:endofrange],id=currentidentifier)
            SeqIO.write(currentsequence,fastafile,"fasta")

        fastafile.close()
        
InSilicoDigest(str(args.fasta),str(args.recognition),str(args.distance5prime),str(args.length))   
