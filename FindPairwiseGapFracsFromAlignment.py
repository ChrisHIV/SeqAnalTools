#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script, taking an sequence alignment in fasta
format as input, considers all possible pairs of sequences therein and
calculates the gap fraction of each sequence if those two sequences were aligned
on their own, i.e. not counting any position at which both sequences have a gap.
(Instead of actually removing such positions, the calculation is made using a
look-up table of hashed gap coordinates, easing the pain of the unavoidable
O(N^2) runtime for N sequences by improving speed 50-fold.)
'''

################################################################################
# USER INPUT
NumHistBins=30
XaxisLabel='gap fraction'
YaxisLabel='unit-normalised frequency'
title=''
xlabelfontsize = 17
ylabelfontsize = 17
titlefontsize = 17
tickfontsize = 13.5
legendfontsize = 17
colour = 'blue'
################################################################################

import argparse
import os
import sys
import itertools
from Bio import AlignIO
from Bio import Seq  
from Bio import SeqIO  
import collections
try:
  import matplotlib.pyplot as plt
except:
  GotMatPlotLib = False
else:
  GotMatPlotLib = True

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('OutputPDFname')
args = parser.parse_args()

# Read in the alignment
try:
  alignment = AlignIO.read(args.alignment, "fasta")
except:
  print('Problem reading', args.alignment + ':', file=sys.stderr)
  raise
AlignmentLength = alignment.get_alignment_length()

if len(alignment) < 2:
  print('Need at least two sequences. Quitting.', file=sys.stderr)
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
  EndOfSeq = len(seq)-1
  while seq[EndOfSeq] in ['?', '-']:
    EndOfSeq -= 1
  return StartOfSeq, EndOfSeq


# Find the start, end, and internal gap coordinates of all sequences.
SeqDict = collections.OrderedDict()
for seq in alignment:

  # Check names are unique
  if seq.id in SeqDict:
    print('Error: encountered sequence', seq.id, 'twice in', args.alignment + \
     '. Sequence names should be unique. Quitting.', file=sys.stderr)
    exit(1)

  ThisSeq = str(seq.seq)
  ThisSeqStart, ThisSeqEnd = FindSeqStartAndEnd(seq.id, ThisSeq)
  InternalGaps = []
  for position,base in enumerate(ThisSeq[ThisSeqStart:ThisSeqEnd+1]):
    if base == '-':
      InternalGaps.append(position + ThisSeqStart)
  SeqDict[seq.id] = [ThisSeqStart, ThisSeqEnd, InternalGaps, set(InternalGaps)]

def CountExternalGaps(SeqStart1, SeqEnd1, InternalGaps2):
  "Counts how many gaps are inside one seq and outside another."
  NumInternalGaps = len(InternalGaps2)
  i = 0
  while InternalGaps2[i] < SeqStart1:
    i += 1
    if i == NumInternalGaps:
      # Every gap from 2 is before the start of 1
      return NumInternalGaps
  j = 0
  while InternalGaps2[len(InternalGaps2) - 1 -j ] > SeqEnd1:
    j += 1
    if j == NumInternalGaps:
      # Every gap from 2 is after the end of 1
      return NumInternalGaps
  return i + j
  
count = 0
total = (len(SeqDict) * (len(SeqDict) - 1)) / 2
output = []
for i, (ID1, [SeqStart1, SeqEnd1, InternalGaps1, InternalGapSet1]) in \
enumerate(SeqDict.items()):
  for ID2, [SeqStart2, SeqEnd2, InternalGaps2, InternalGapSet2] in \
  SeqDict.items()[i+1:]:

    count +=1
    if count % 10000 == 0:
      print('Now working on pair', count, 'of', total, file=sys.stderr)

    # Common gaps inside both sequences:
    NumCommonInternalGaps = len(InternalGapSet1.intersection(InternalGapSet2))

    # Now find any gaps inside one sequence & outside the other:
    NumCommonGapsIn1 = NumCommonInternalGaps + \
    CountExternalGaps(SeqStart2, SeqEnd2, InternalGaps1)
    NumCommonGapsIn2 = NumCommonInternalGaps + \
    CountExternalGaps(SeqStart1, SeqEnd1, InternalGaps2)

    GapFrac1 = float(len(InternalGaps1) - NumCommonGapsIn1) / \
    (SeqEnd1 - SeqStart1 + 1 - NumCommonGapsIn1)
    GapFrac2 = float(len(InternalGaps2) - NumCommonGapsIn2) / \
    (SeqEnd2 - SeqStart2 + 1 - NumCommonGapsIn2)

    output.append([GapFrac1, GapFrac2, ID1, ID2])

AllGapFracs = list(itertools.chain.from_iterable(x[:2] for x in output))

print('\n'.join([','.join(map(str, data)) for data in output]))

if not GotMatPlotLib:
  print('Failed to import matplotlib. Quitting.', file=sys.stderr)  
  exit(0)

ax = plt.figure().add_subplot(111)
plt.hist(AllGapFracs, bins=NumHistBins, normed=True, color=colour)
plt.xlabel(XaxisLabel, fontsize=xlabelfontsize)
plt.ylabel(YaxisLabel, fontsize=ylabelfontsize)
ax.tick_params(axis='both', which='major', labelsize=tickfontsize)
plt.title(title, fontsize=titlefontsize)
plt.savefig(args.OutputPDFname)



# For testing: a method that's slower but easier to understand.
'''
output = []
for i,seq1 in enumerate(alignment):
  ID1 = seq1.id
  for j,seq2 in enumerate(alignment[i+1:]):
    PureGapPositions = []
    ID2 = seq2.id
    SeqAsString1 = str(seq1.seq)
    SeqAsString2 = str(seq2.seq)
    for pos,(base1,base2) in enumerate(itertools.izip(SeqAsString1, \
    SeqAsString2)):
      if base1 == '-' == base2:
        PureGapPositions.append(pos)
    PureGapPositions.reverse()
    for pos in PureGapPositions:
      SeqAsString1 = SeqAsString1[:pos] + SeqAsString1[pos+1:]
      SeqAsString2 = SeqAsString2[:pos] + SeqAsString2[pos+1:]
    start1, end1 = FindSeqStartAndEnd(ID1, SeqAsString1)
    start2, end2 = FindSeqStartAndEnd(ID2, SeqAsString2)
    GapFrac1 = float(SeqAsString1[start1:end1+1].count('-'))/(end1 - start1 + 1)
    GapFrac2 = float(SeqAsString2[start2:end2+1].count('-'))/(end2 - start2 + 1)
    output.append([GapFrac1, GapFrac2, ID1, ID2])

print('\n'.join([','.join(map(str, data)) for data in output]))
'''
