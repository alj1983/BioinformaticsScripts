import re # to use regular expressions
import math
import argparse # To allow users setting flags in order to define what to filter for.
from pylab import figure, title, xlabel, ylabel, hist, axis, grid, savefig # to draw a historgram

parser = argparse.ArgumentParser(description="Creates a cumulative histogram on read coverage (excluding coverage of 0) and prints out upper and lower coverage thresholds")
parser.add_argument("-p","--pileup", help="pileup file name")
parser.add_argument("-ma","--max", help="maximum coverage to be shown")
parser.add_argument("-mi","--min", help="minimum coverage to be shown")
args=parser.parse_args()
# Define the search patterns:

PileupFileNamePattern = re.compile('(.*)\.pileup')



def Coverages(pileupfile,maximum,minimum):
     '''
     Creates a cumulative histogram on read coverage (excluding coverage of 0) and prints out upper and lower coverage thresholds
     '''
     pileup=open(pileupfile)
     OutfilePrefix1 =PileupFileNamePattern.match(pileupfile) # define outfile 
     OutfilePrefix=OutfilePrefix1.group(1)      
     Coverage=[] #Create an empty list
     for line in pileup:
         line=line.strip()
         splitted=line.split('\t') # Split the lines by tab separators
         if abs(int(splitted[3])) > 0: # excluding all coverages of 0
              Coverage.append(abs(int(splitted[3]))) # Pick out the 4th column, which is the number of reads that covered this position
     # Remove coverages of 0
     CoverageAboveZero=[]
     for c in Coverage:
         if c>0:              
              CoverageAboveZero.append(c)
     # sort the coverages
     CoverageAboveZeroSorted = CoverageAboveZero.sort()
     TotalBasesCovered = len(CoverageAboveZeroSorted)
     # Identify the lower 2% limit
     LowerLimit=CoverageAboveZeroSorted[int((0.02*TotalBasesCovered)-1)]
     # Identify the upper 2% limit
     UpperLimit=CoverageAboveZeroSorted[int((0.98*TotalBasesCovered)-1)]
     # Calculating the standard deviation
     mean=float(sum(CoverageAboveZeroSorted)/len(CoverageAboveZeroSorted))
     total = 0.0
     for value in CoverageAboveZeroSorted:
          total += (value - mean) ** 2
     stddev = math.sqrt(total / (len(CoverageAboveZeroSorted)-1))
     print pileupfile, 'Coverage -           mean:', '%10.2f' % mean
     print pileupfile, 'Coverage -            max:', '%10i' % max(CoverageAboveZeroSorted)
     print pileupfile, 'Coverage -            min:', '%10i' % min(CoverageAboveZeroSorted) 
     print pileupfile, 'Coverage -  mean-2*stddev:', '%10i' % mean-2*stddev
     print pileupfile, 'Coverage -  mean+2*stddev:', '%10i' % mean+2*stddev
     print pileupfile, 'Coverage -       upper 2%:', '%10i' % UpperLimit
     print pileupfile, 'Coverage -       lower 2%:', '%10i' % LowerLimit

     # Drawing the historgram
     figure()
     hist(CoverageAboveZeroSorted, range=(minimum,maximum), bins=50, histtype='bar', facecolor='green', normed=1,cumulative=1)
     title('Cumulative Histogram')
     xlabel('Coverage')
     ylabel('Frequency')
     axis()
     grid(True)
     savefig(OutfilePrefix + 'Coverage_Histogram.png')

     

Coverages(str(args.pileup),int(args.max),int(args.min))
