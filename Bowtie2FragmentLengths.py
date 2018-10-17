import re # to use regular expressions
import argparse # To allow users setting flags in order to define what to filter for.
import math # to calculate sqrt
from pylab import figure, title, xlabel, ylabel, hist, axis, grid, savefig # to draw a historgram

parser = argparse.ArgumentParser(description="Creates a histogram on fragment lengths covered by paired end reads")
parser.add_argument("-s","--sam", help="sam file name")
args=parser.parse_args()
# Define the search patterns:

HeaderPattern= re.compile('^@')
SamFileNamePattern = re.compile('(.*)\.sam')



def Fraglength(samfile):
     '''
     Creates a histogram on fragment lengths covered by paired end reads
     '''
     sam=open(samfile)
     OutfilePrefix1 =SamFileNamePattern.match(samfile) # define outfile 
     OutfilePrefix=OutfilePrefix1.group(1)      
     Lengths=[] #Create an empty list
     for line in sam:
         line=line.strip()
         if not HeaderPattern.search(line): # If this is a header line then skip it
              splitted=line.split('\t') # Split the lines by tab separators
              Lengths.append(abs(int(splitted[8]))) # Pick out the 9th column, which is the inferred fragment length (output by the bowtie2 mapper) and get the absolute integer value from it
     mean=float(sum(Lengths)/len(Lengths))
     total = 0.0
     for value in Lengths:
          total += (value - mean) ** 2
     stddev = math.sqrt(total / len(Lengths))
     stderr = stddev / math.sqrt(len(Lengths))

     print samfile, 'Fragment length -           mean:', '%10.2f' % mean
     print samfile, 'Fragment length - standard error:', '%10.2f' % stderr

     # Drawing the historgram
     n_bins = 10
     figure()
     num, bins, patches = hist(Lengths, n_bins, histtype='bar', facecolor='green')
     title('Histogram')
     xlabel('Fragment length')
     ylabel('Frequency')
     axis()
     grid(True)
     savefig(OutfilePrefix + 'FragmentLength_Histogram.png')

     

Fraglength(str(args.sam))