#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script removes positions, specified with respect to
one sequence therein, from an alignment.'''

import argparse
import os
import sys
import collections
from Bio import SeqIO
from Bio import Seq

GapChar = '-'

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
parser.add_argument('ReferenceName', help='The name of the sequence with '+\
'respect to which the positions to remove are specified.')
parser.add_argument('PositionToRemove', type=int, nargs='+')
args = parser.parse_args()

# Order the coordinates. Check that each was specified only once. The smallest 
# should be greater than 0.
PositionsToRemove = args.PositionToRemove
if min(PositionsToRemove) < 1:
  print('position args should be greater than zero. Quitting.', file=sys.stderr)
  exit(1)
CounterObject = collections.Counter(PositionsToRemove)
DuplicatedPositionsToRemove = [pos for pos in CounterObject if \
CounterObject[pos] > 1]
if len(DuplicatedPositionsToRemove) != 0:
  print('Warning: the following coordinates were multiply specified and have', 
  'only been counted once:', ', '.join(map(str,DuplicatedPositionsToRemove)),
  file=sys.stderr)
  PositionsToRemove = list(set(PositionsToRemove))

# Read in the alignment
AlignmentLength = None
seqs = []
ref = None
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  SeqLength = len(seq.seq)
  if AlignmentLength == None:
    AlignmentLength = SeqLength
  elif SeqLength != AlignmentLength:
    print('Sequence', seq.id, 'has length', str(SeqLength) +', c.f. the first',\
    'sequence in', args.FastaFile, 'which has length', str(AlignmentLength)+\
    '. Aligned sequences are required. Quitting.', file=sys.stderr)
    exit(1)
  if seq.id == args.ReferenceName:
    if ref != None:
      print('The sequence name', seq.id, 'is not unique in', args.FastaFile + \
      '. Quitting.', file=sys.stderr)
      exit(1)
    ref = seq
    RefLength = len(seq.seq.ungap(GapChar))
  seqs.append(seq)
if ref == None:
  print('Sequence', args.ReferenceName, 'was not found in', args.FastaFile +\
  '. Quitting.', file=sys.stderr)
  exit(1)

# Check no position to remove is after the end of the desired ref.
if max(PositionsToRemove) > RefLength:
  print('The largest position arg,', str(max(PositionsToRemove)) +\
  ', is greater than the length of', args.ReferenceName +'. Quitting.', \
  file=sys.stderr)
  exit(1)

# Find the (zero-based) coordinates with respect to the alignment
PositionsToRemove_AlignmentCoords = []
PositionInRef = 0
for AlignmentPostitionMin1,base in enumerate(str(ref.seq)):
  if base == GapChar:
    continue
  PositionInRef += 1
  if PositionInRef in PositionsToRemove:
    PositionsToRemove_AlignmentCoords.append(AlignmentPostitionMin1)

# Iterate backwards through the coordinates, remove them from each sequence.
PositionsToRemove_AlignmentCoords = \
sorted(PositionsToRemove_AlignmentCoords, reverse=True)
for seq in seqs:
  StrippedSeq = str(seq.seq)
  for position in PositionsToRemove_AlignmentCoords:
    if seq.id == args.ReferenceName:
      print('removing base', StrippedSeq[position], file=sys.stderr)
    StrippedSeq = StrippedSeq[:position] + StrippedSeq[position+1:]
  seq.seq = Seq.Seq(StrippedSeq)

SeqIO.write(seqs, sys.stdout, "fasta")
