#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script finds those read pairs in a bam file for
 which the insert size is a user-specified value. Amongst these read pairs, the
 number of times each unique read occurs is counted. The most common reads (call
 with the -h flag for more detail) are printed to stdout in fasta format. The
intended purpose is to find cases where one read pair has been exactly
duplicated a large number of times, leading to a spike in an otherwise smooth
insert-size distribution.
'''

import pysam
import argparse
import os
from subprocess import call

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BamFile')
parser.add_argument('InsertSize', type=int)
parser.add_argument('RelativeCountThreshold', type=float, \
help='The desired reads will be sorted by how many times they occur; only '+\
'those whose count is at least the maximum count times this fraction will be '\
+'kept. e.g. If the most common read occurs 1000 times and this fraction is'+\
' set to be 0.15, only those reads appearing at least 150 times will be kept.')
args = parser.parse_args()

# Check that the BamFileName exists and is a file
if not os.path.isfile(args.BamFile):
  print(args.BamFile, 'does not exist or is not a file.', file=sys.stderr)
  exit(1)

# If the BamFileName does not have an index file, make one.
if not os.path.isfile(args.BamFile+'.bai'):
  call(['samtools', 'index', args.BamFile])

basename = os.path.basename(args.BamFile).rsplit('.',1)[0]
BamFile = pysam.AlignmentFile(args.BamFile, "rb")

# Check that the RelativeCountThreshold is between 0 and 1 inclusive.
if not 0 <= args.RelativeCountThreshold <= 1:
  print('The RelativeCountThreshold should be between 0 and 1 inclusive.'+\
  '\nQuitting', file=sys.stderr)
  exit(1)

# For inserts of the desired size, count how many times each read occurs
UniqueReadCounts = {}
for read in BamFile.fetch():
  if read.is_paired and abs(read.tlen) == args.InsertSize:
    if read.query_sequence in UniqueReadCounts:
      UniqueReadCounts[read.query_sequence] += 1
    else:
      UniqueReadCounts[read.query_sequence] = 1

# Print the most common reads
LastCount = None
for read, count in \
sorted(UniqueReadCounts.items(), key=lambda x: x[1], reverse=True):
  if LastCount != None and float(count)/LastCount < args.RelativeCountThreshold:
    break
  SeqHeader = '>'+basename+'_'+str(count)
  print(SeqHeader+'\n'+read)
  LastCount = count

