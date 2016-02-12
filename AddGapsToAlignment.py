#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script adds gaps into an alignment. The coordinates
are moved through from smallest to largest; at each coordinate, that column of
the alignment is moved one position to the right and a column of gaps is
inserted. For example, specifying the coordinates 1, 3 and 4 would insert one 
column of gaps onto the left of the alignment, and two columns of gaps in
between the first and second columns.'''

import argparse
import os
import sys
import collections
from Bio import Seq
from Bio import SeqIO

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
parser.add_argument('GapCoord', type=int, nargs='+')
args = parser.parse_args()

# Order the coordinates. Check that each was specified only once. The smallest 
# should be greater than 0.
coords = sorted(args.GapCoord)
if min(coords) < 1:
  print('Coordinates should be greater than zero. Quitting.', file=sys.stderr)
  exit(1)
CounterObject = collections.Counter(coords)
DuplicatedCoords = [coord for coord in CounterObject if \
CounterObject[coord] > 1]
if len(DuplicatedCoords) != 0:
  print('The following coordinates were multiply specified:', \
  ' '.join(DuplicatedArgs) +'; all arguments should be unique. Quitting.', \
  file=sys.stderr)
  exit(1)

# Read in the alignment
AlignmentLength = None
seqs = []
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  SeqLength = len(seq.seq)
  if AlignmentLength == None:
    AlignmentLength = SeqLength
  elif SeqLength != AlignmentLength:
    print('Sequence', seq.id, 'has length', str(SeqLength) +', c.f. the first',\
    'sequence in', args.FastaFile, 'which has length', str(AlignmentLength)+\
    '. Aligned sequences are required. Quitting.', file=sys.stderr)
    exit(1)
  seqs.append(seq)
NumSeqs = len(seqs)
if NumSeqs == 0:
  print('No sequences found in', args.FastaFile +'. Quitting.', file=sys.stderr)
  exit(1)

# Make the coords zero-based
coords = [coord-1 for coord in coords]

# Add the gaps to each sequence!
for i in range(0,NumSeqs):
  SeqAsString = str(seqs[i].seq)
  for NumGapsAddedSoFar,coord in enumerate(coords):
    SeqLengthSoFar = AlignmentLength + NumGapsAddedSoFar
    if coord < SeqLengthSoFar:
      SeqAsString = SeqAsString[:coord] + GapChar + SeqAsString[coord:]
    elif coord == SeqLengthSoFar:
      SeqAsString += GapChar
    else:
      print('Coordinate', coord+1, 'does not make sense: the input alignment was'\
      ,str(AlignmentLength)+'bp long, and', NumGapsAddedSoFar, \
      'gaps have been added so far. Quitting.', file=sys.stderr)
      exit(1)
  seqs[i].seq = Seq.Seq(SeqAsString)

SeqIO.write(seqs, sys.stdout, "fasta")
