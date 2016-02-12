#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Calculates the mean quality, as a function of position
within the read (since quality typically drops at the edges), of all reads found
in all fastq.gz files passed as arguments.'''

import os
import sys
import argparse
import itertools as it
#import matplotlib.pyplot as plt
from Bio import SeqIO
import gzip

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
#parser.add_argument('OutputPDFname')
parser.add_argument('FastqFile', type=File, nargs='+')
args = parser.parse_args()

#for FastqFile in args.FastqFile:
#  if not ((len(FastqFile) > 3 and FastqFile[-3:] == '.fq') or \
#  (len(FastqFile) > 6 and FastqFile[-6:] == '.fastq')):
#    print('Warning: at least one input file ends in something other than .fq', \
#    ' or .fastq. If a file is not in fastq format, the python Bio.SeqIO',\
#    'module produces an incomprehensible error message. Note that zipped',\
#    'fastqs must be unzipped first.', file=sys.stderr)
#    break

SumsOfAllQs = []
PositionCounts = []
print('Of', len(args.FastqFile), 'files, number processed: ', end='')
for NumFilesDone, FastqFile in enumerate(args.FastqFile):
  if NumFilesDone % 10 == 0:
    print(NumFilesDone, end=' ')
    sys.stdout.flush()
  with gzip.open(FastqFile, 'r') as f:
    for read in SeqIO.parse(f, 'fastq'):
      qualities = read.letter_annotations["phred_quality"]
      SumsOfAllQs = \
      [x+y for x,y in it.izip_longest(SumsOfAllQs, qualities, fillvalue=0)]
      PositionCounts = [x+y for x,y in it.izip_longest(PositionCounts, \
      [1]*len(qualities), fillvalue=0)]
print()

if len(SumsOfAllQs) != len(PositionCounts):
  print('Malfunction of the code: len(SumsOfAllQs) != len(PositionCounts)' +\
  '. Quitting.', file=sys.stderr)
  exit(1)

MeanQs = [float(total)/count for total,count in \
it.izip(SumsOfAllQs, PositionCounts)]

for PosMin1, q in enumerate(MeanQs):
  print(PosMin1+1, q)

#plt.plot(range(1,len(MeanQs)+1), MeanQs)
#plt.xlabel('position of base in read')
#plt.ylabel('mean quality')
#plt.savefig(args.OutputPDFname)

