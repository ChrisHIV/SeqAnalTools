#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script trims gap characters "-" and "?" from the
ends of all sequences in a fasta file. (Note that if the fasta file is an
alignment, this operation will, in general, unalign the sequences.)'''

GapChars = '-?'

import argparse
import os
import sys
from Bio import SeqIO
from collections import OrderedDict

# The following stops an error message from a python compiler bug when the
# output of this script is piped to something that does not read all the output,
# e.g.  awk '/^>/{if(N)exit;++N;} {print;}'  which prints only the first
# sequence.
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

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

def FindSeqStartAndEnd(SeqName,seq):
  '''Find the 0-based positions of the start and end of the sequence.'''
  StartOfSeq = 0
  try:
    while seq[StartOfSeq] in GapChars:
      StartOfSeq += 1
  except IndexError:
    print(SeqName, "has no bases - it's just one big gap. Quitting.", \
    file=sys.stderr)
    exit(1)
  EndOfSeq = len(seq)-1
  while seq[EndOfSeq] in GapChars:
    EndOfSeq -= 1
  return StartOfSeq,EndOfSeq

# Thanks Stackoverflow:
def insert_newlines(string, every=50):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

SeqDict = OrderedDict()
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  if seq.id in SeqDict:
    print('Two (or more) sequences in', args.FastaFile, 'are called', seq.id+\
    '. Sequence names should be unique. Quitting.', file=sys.stderr)
    exit(1)
  SeqAsString = str(seq.seq)
  StartOfSeq,EndOfSeq = FindSeqStartAndEnd(seq.id, SeqAsString)
  SeqDict[seq.id] = SeqAsString[StartOfSeq:EndOfSeq+1]

for SeqName,seq in SeqDict.items():
  print('>'+SeqName+'\n'+insert_newlines(seq))




