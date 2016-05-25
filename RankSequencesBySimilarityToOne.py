#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Calculates the similarity of every other sequence to
one named sequence in a fasta-format alignment. Comparing each other sequence in
turn we find the number of matching bases in the region of the alignment
covered by both sequences (ignoring missing coverage - "?" positions), the
length of that region, and the ratio of the former to the latter i.e. the
fractional identity in the overlap. Results are printed to stdout, in csv
format, sorted by this ratio. (For non-overlap this ratio is 0/0: we report 0.)
'''

import argparse
import os
import sys
import itertools
from Bio import AlignIO
from Bio import Seq  
from Bio import SeqIO  
import collections

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('ChosenSeq')
parser.add_argument('-C', '--case-sensitive', action='store_true', \
help='consider upper- and lower-case occurences of the same base to be '+\
'unequal, e.g. A != a. (By default, case is ignored.)')
parser.add_argument('-I', '--indels-only', action='store_true', \
help='Only count differences in insertions and deletions as differences, '+\
'not base changes.')
args = parser.parse_args()

# For brevity & speed
IgnoreCase = not args.case_sensitive

# Read in the alignment
try:
  alignment = AlignIO.read(args.alignment, "fasta")
except:
  print('Problem reading', args.alignment + ':', file=sys.stderr)
  raise
AlignmentLength = alignment.get_alignment_length()

# Find the chosen seq; make sure it's there and uniquely named.
NumSeqsWithMatchingName = 0
for seq in alignment:
  if seq.id == args.ChosenSeq:
    ChosenSeq = str(seq.seq)
    if IgnoreCase:
      ChosenSeq = ChosenSeq.upper()
    NumSeqsWithMatchingName += 1
if NumSeqsWithMatchingName != 1:
  print('Error: there are', NumSeqsWithMatchingName, 'sequences named', \
  args.ChosenSeq, 'in', args.alignment + '. There should be exactly 1.',
  'Quitting.', file=sys.stderr)
  exit(1)

def FindSeqStartAndEnd(SeqName,seq):
  '''Find the 0-based positions of the start and end of the sequence.'''
  StartOfSeq = 0
  try:
    while seq[StartOfSeq] in ['?', '-']:
      StartOfSeq += 1
  except IndexError:
    print(SeqName, "has no bases - it's just one big gap. Quitting.", \
    file=sys.stderr)
    exit(1)
  EndOfSeq = AlignmentLength-1
  while seq[EndOfSeq] in ['?', '-']:
    EndOfSeq -= 1
  return StartOfSeq, EndOfSeq

ChosenSeqStart, ChosenSeqEnd = FindSeqStartAndEnd(args.ChosenSeq, ChosenSeq)

# Iterate through all other sequences.
SimilaritiesDict = {}
for seq in alignment:
  if seq.id == args.ChosenSeq:
    continue

  # Check names are unique
  if seq.id in SimilaritiesDict:
    print('Error: encountered sequence', seq.id, 'twice in', args.alignment + \
     '. Sequence names should be unique. Quitting.', file=sys.stderr)
    exit(1)

  # Calculate similarity
  ThisSeq = str(seq.seq)
  if IgnoreCase:
    ThisSeq = ThisSeq.upper()
  ThisSeqStart, ThisSeqEnd = FindSeqStartAndEnd(seq.id, ThisSeq)
  OverlapStart = max(ThisSeqStart, ChosenSeqStart)
  OverlapEnd   = min(ThisSeqEnd,   ChosenSeqEnd)
  NumAgreeing = 0
  overlap = 0
  if OverlapEnd >= OverlapStart:
    for ChosenSeqBase, ThisSeqBase in itertools.izip( \
    ChosenSeq[OverlapStart:OverlapEnd+1], ThisSeq[OverlapStart:OverlapEnd+1]):
      if ChosenSeqBase == '?' or ThisSeqBase == '?':
        continue
      overlap += 1
      if ChosenSeqBase == ThisSeqBase:
        NumAgreeing += 1
      elif args.indels_only and ChosenSeqBase != '-' and ThisSeqBase != '-':
        NumAgreeing += 1
  if overlap == 0:
    FractionalSimilarity = 0
  else:
    FractionalSimilarity = float(NumAgreeing)/overlap
  SimilaritiesDict[seq.id] = [NumAgreeing, overlap, FractionalSimilarity]

# Sort by fractional identity
output = '"sequence","number of positions agreeing with ' + args.ChosenSeq + \
' in their overlap","overlap length","fraction of positions agreeing with ' + \
args.ChosenSeq + ' in their overlap (ratio of columns 2 to 3)"'
for ID, values in sorted(SimilaritiesDict.items(), key=lambda x: x[1][2], \
reverse=True):
  output += '\n' + ID +',' + ','.join(map(str,values))
print(output)
