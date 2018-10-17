#!/usr/bin/python

import re # to use regular expressions
#import sys # to read user arguments
import argparse # To allow users setting flags in order to define what to filter for.

parser = argparse.ArgumentParser(description="Filters bowtie2 alignments (.sam files)")

parser.add_argument("-mq", "--MappingQuality", help="Filter out reads with a mapping quality below 20",
                    action="store_true")
parser.add_argument("-u", "--unmapped", help="Filter out unmapped reads",
                    action="store_true")
parser.add_argument("-d", "--discordant", help="Filter out discordantly mapped reads",
                    action="store_true")
parser.add_argument("-a", "--ambiguous", help="Filter out ambiguous (non-unique) mappings",
                    action="store_true")
parser.add_argument("-s","--sam", help="sam file name")
args=parser.parse_args()
# Define the search patterns:

HeaderPattern= re.compile('^@')
DiscordantPattern= re.compile('YT:Z:DP') # This field is specific to mappings that were done with Bowtie2
NonUniquePattern= re.compile('XS:i') # This field is specific for mappings that mapped more than once
SamFileNamePattern = re.compile('(.*)\.sam')


def MappingFilter(samfile):
     '''
     Filters
     the percentage of reads with a mapping quality below 20;
     the percentage of unmapped reads;
     the percentage of discordantly mapped reads;
     and the percentage of non-unique mappings.
     '''
     sam=open(samfile)
     OutfilePrefix1 =SamFileNamePattern.match(samfile) # define outfile 
     OutfilePrefix=OutfilePrefix1.group(1)      
     outfile=open(OutfilePrefix + 'filtered.sam','w')
     total=0
     filtered=0
     for line in sam:
         line=line.strip()
         if HeaderPattern.search(line): # If this is a header line just write it to the outfile
              outfile.write(line+'\n')
         else: # if there is no @ at the beginning of the line
              tag='keep' # Keep the line by default, and only throw it out if any of the conditions below are not met
              total=total+1 # counting the total number of reads
	      splitted=line.split('\t') # Split the lines by tab separators
              if args.MappingQuality: # if the user wants to filter out reads below a mapping quality of 20
                   if int(splitted[4])>19: # The fifth column -- which is the mapping quality -- shall be over 19
	                tag='keep'
                   else: 
                        tag='discard'
              if tag is 'keep': # Continue only if the alignment is not yet filtered out
                   if args.unmapped: # if the unmapped reads shall be filtered out
                        if splitted[2] is not '*': # column 3 in the SAM file shows "Name of reference sequence where alignment occurs" or a * for unmapped reads.
	                     tag='keep'
                        else: 
                             tag='discard'
                   if tag is 'keep': # Continue only if the alignment is not yet filtered out
                        if args.discordant: # if the discordant reads shall be filtered out
                             if DiscordantPattern.search(line): # Look for the field indicating discordant mappings
                                  tag='discard'
                             else: 
                                  tag='keep'
                        if tag is 'keep': # Continue only if the alignment is not yet filtered out
                             if args.ambiguous: # if the ambiguous reads shall be filtered out
                                  if NonUniquePattern.search(line): # Look for the field indicating discordant mappings
                                       tag='discard'
                                  else: 
                                       tag='keep'
                             if tag is 'keep': # Write the alignment only if it was not filtered out before
                                  filtered = filtered+1 # count the number of reads remaining after filtering
                                  outfile.write(line+'\n')
     print samfile, ' - unfiltered reads:', '%13i' % total
     print samfile, ' -   filtered reads:', '%13i' % filtered
     

MappingFilter(str(args.sam))