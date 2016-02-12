#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script takes one fasta file as an argument and
prints the length of each sequence found therein, optionally ignoring gaps.'''

import argparse
import os
import sys
from Bio import SeqIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
parser.add_argument('-g', '--include-gaps', action='store_true', \
help='include gap characters, "-" and "?", ignored by default')
args = parser.parse_args()


SeqLengths = []
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  if not args.include_gaps:
    seq.seq = seq.seq.ungap("-").ungap("?")
  SeqLengths.append([seq.id, len(seq.seq)])

for [SeqName,SeqLength] in SeqLengths:
  print(SeqName, SeqLength)


