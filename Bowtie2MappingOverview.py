import re # to use regular expressions
import sys # to read user arguments

# Define the search patterns:

HeaderPattern= re.compile('^@')
DiscordantPattern= re.compile('YT:Z:DP') # This field is specific to mappings that were done with Bowtie2
NonUniquePattern= re.compile('XS:i') # This field is specific for mappings that mapped more than once
outfile=open('out.txt','w')
headline= ['sequence','total','mapq20','unmapped','discordant','nonunique']
headlineout='\t'.join(headline)
outfile.write(headlineout)
outfile.write('\n')

def MappingOverview(samfile):
     '''
     Counts the total number of aligned reads;
     the percentage of reads with a mapping quality below 20;
     the percentage of unmapped reads;
     the percentage of discordantly mapped reads;
     and the percentage of non-unique mappings.
     '''
     sam=open(samfile)
     total=0
     mapq20=0
     unmapped=0
     discordant=0
     nonunique=0
     for line in sam:
         line=line.strip()
         if not HeaderPattern.search(line): # if there is no @ at the beginning of the line
             total=total+1
	     splitted=line.split('\t') # Split the lines by tab separators
             if int(splitted[4])>19: # The fifth column -- which is the mapping quality -- shall be over 19
	          mapq20=mapq20+1
             if splitted[2] is '*': # column 3 in the SAM file shows "Name of reference sequence where alignment occurs" or a * for unmapped reads.
	          unmapped=unmapped+1
             if DiscordantPattern.search(line): # Look for the field indicating discordant mappings
                  discordant=discordant+1
             if NonUniquePattern.search(line): # Look for the field indicating concordant mappings
                  nonunique=nonunique+1
     outfile.write('%s\t%12i\t%7.6f\t%7.6f\t%7.6f\t%7.6f\n' % (samfile, total, float(mapq20)/float(total), float(unmapped)/float(total), float(discordant)/float(total), float(nonunique)/float(total)))
     

for s in sys.argv[1:]:
     MappingOverview(str(s))         