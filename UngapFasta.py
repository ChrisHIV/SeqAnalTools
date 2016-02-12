#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script removes the gap character "-" from sequences
in a fasta file. Output is printed to stdout.'''

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
args = parser.parse_args()

UngappedSeqs = []
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  seq.seq = seq.seq.ungap("-")
  UngappedSeqs.append(seq)

SeqIO.write(UngappedSeqs, sys.stdout, "fasta")

