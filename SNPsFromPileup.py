import re # to use regular expressions
import argparse # To allow users setting flags in order to define what to filter for.

parser = argparse.ArgumentParser(description="Reduces pileup files to specified SNP locations")
parser.add_argument("-s","--snps", help="text file containing SNP locations. Each line contains the chromosome specification and the TAB-separated position of a SNP")
parser.add_argument("-p","--pileup", help="pileup file, created with samtools mpileup")
args=parser.parse_args()


def ExtractingSNPs(snpfile):
    '''
    Stores the chromosome specification and the position of the SNPs into a table;
    '''
     
    snplocations=open(snpfile)
    #chrompos_table = []
    chrompos_list = []
    for line in snplocations:
        line=line.strip()
        columns = line.split('\t')
        #chrompos_table.append([columns[0],columns[1]])
        chrompos_list.append(str(columns[0]+columns[1]))
        
    return chrompos_list

SNPlist = ExtractingSNPs(str(args.snps))

# Extract the SNP locations from the pileup file
outfile=open(str(args.pileup)+'_ReducedToSNPs.pileup','w')

def SNPsInPileup(pileupfile):
      '''
      Reads the pileup file and extracts only those lines that contain SNPs
      '''
      pileups=open(pileupfile)
      for line in pileups:
           line=line.strip()
           columns = line.split('\t')
           chrompos= str(columns[0]+columns[1])
           if chrompos in SNPlist:
                lines = '\t'.join(columns)
                outfile.write(lines+'\n')

SNPsInPileup(str(args.pileup))

outfile.close()
