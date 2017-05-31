#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Given an alignment of sequences, this script prints the
(one-based) coordinates of positions (if any exist) at which a named sequence
has unique insertions, i.e. positions where only that sequence has bases, all
other sequences have gaps. '''

import argparse
import os
import sys
from Bio import SeqIO

GapChar = '-'

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
parser.add_argument('NameOfSeqOfInterest')
args = parser.parse_args()

# Read in the alignment
AlignmentLength = None
ChosenRef = None
OtherRefs = []
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  SeqAsString = str(seq.seq)
  if AlignmentLength == None:
    AlignmentLength = len(SeqAsString)
  elif len(SeqAsString) != AlignmentLength:
    print(args.FastaFile, 'is not an alignment. Quitting.', file=sys.stderr)
    exit(1)
  if seq.id == args.NameOfSeqOfInterest:
    ChosenRef = SeqAsString
  else:
    OtherRefs.append(SeqAsString)

# Check we found the chosen ref and at least one other
if ChosenRef == None:
  print(args.NameOfSeqOfInterest, 'not found in', args.FastaFile +'. Quitting.',
  file=sys.stderr)
  exit(1)
if OtherRefs == []:
  print('No sequences other than', args.NameOfSeqOfInterest, 'found in', \
  args.FastaFile +'. Quitting.', file=sys.stderr)
  exit(1)

UniqueInsertionPositions = []
for position,ChosenRefBase in enumerate(ChosenRef):
  if ChosenRefBase != GapChar:
    NoOtherRefHasBaseHere = True
    for OtherRefSeq in OtherRefs:
      if OtherRefSeq[position] != GapChar:
        NoOtherRefHasBaseHere = False
        break
    if NoOtherRefHasBaseHere:
      UniqueInsertionPositions.append(position+1)

print(' '.join(map(str, UniqueInsertionPositions)))

