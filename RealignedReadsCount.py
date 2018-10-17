import re # to use regular expressions
import argparse # To allow users setting flags in order to define what to filter for.


parser = argparse.ArgumentParser(description="Counts the total number of aligned reads and the number of reads that were realigned with GATK (based on the OC tag)")
parser.add_argument("-s","--sam", help="sam file name")
args=parser.parse_args()

# Define the search patterns:
HeaderPattern= re.compile('^@')
RealignedPattern= re.compile('OC') # This field is a specific tag storing the CIGAR string before realignment was performed

def RealignedReadsCount(samfile):
     '''
     Counts the total number of aligned reads;
     Count the number of realigned reads;
     '''
     sam=open(samfile)
     total=0
     realigned=0
     for line in sam:
         line=line.strip()
         if not HeaderPattern.search(line): # if there is no @ at the beginning of the line
             total=total+1
             if RealignedPattern.search(line): # Look for the field indicating realigned reads
                  realigned=realigned+1

     print samfile, '    total number of reads:', '%10i' % total
     print samfile, 'number of realigned reads:', '%10i' % realigned
     
RealignedReadsCount(str(args.sam))