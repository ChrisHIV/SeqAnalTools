#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script replaces any "-" character (taken to mean a deletion at this position) that neighbours a "?" character (taken to mean unknown / missing data) by a "?" character, iteratively until no more replacements need to be made. The idea is that a deletion should not be called if it is not known what's next to it. Output is printed to stdout, suitable for redirection to a fasta file.'''

NoCoverageChar='?'
GapChar='-'

import argparse
import os
import sys
from Bio import SeqIO
from Bio import Seq

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

# Where NoCoverageChars neighbour GapChars, propagate the former outwards until
# they touch bases on both sides (because insertions should only be called when
# the bases on either side are known). e.g.
# ACTG---?---ACTG
# becomes
# ACTG???????ACTG
def PropagateNoCoverageChar(seq, LeftToRightDone=False):
  '''Replaces gaps that border "no coverage" by "no coverage".'''
  
  if LeftToRightDone:
    seq = seq[::-1]
  BaseToLeftIsNoCoverage = False
  ResultingSeq = ''
  for base in seq:
    if base == NoCoverageChar:
      BaseToLeftIsNoCoverage = True
      ResultingSeq += NoCoverageChar
    elif base == GapChar:
      if BaseToLeftIsNoCoverage:
        ResultingSeq += NoCoverageChar
      else:
        ResultingSeq += GapChar
    else:
      BaseToLeftIsNoCoverage = False
      ResultingSeq += base
  if LeftToRightDone:
    ResultingSeq = ResultingSeq[::-1]
  else:
    ResultingSeq = PropagateNoCoverageChar(ResultingSeq, True)
  return ResultingSeq

seqs = []
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  SeqAsString = str(seq.seq)
  SeqAsString = PropagateNoCoverageChar(SeqAsString)
  seq.seq = Seq.Seq(SeqAsString)
  seqs.append(seq)

SeqIO.write(seqs, sys.stdout, "fasta")

