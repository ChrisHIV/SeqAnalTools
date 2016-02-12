#!/usr/bin/env python2

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview: call this script from the command line with a single fasta file as
## an argument, that file containing two sequences. This script moves the
## shorter one along the longer one to find the best match.
##
## Output is printed to stdout, suitable for redirection to a target fasta file.

import os.path, sys
from AuxiliaryFunctions import ReadSequencesFromFile

# Check this file is called from the command line with one argument
if len(sys.argv[1:]) != 1:
  print 'Incorrect number of arguments given.'
  print 'Usage:\n', sys.argv[0], 'NameOfYourFastaFile.fasta'
  exit(1)
DataFile = sys.argv[1]

# Check that the argument exists and is a file
if not os.path.isfile(DataFile):
  print DataFile, 'does not exist or is not a file.'
  exit(1)

# Read in the sequences as a dictionary. They are (nominally) not aligned.
Aligned=False
SeqDict, FirstSeqLength = ReadSequencesFromFile(DataFile, Aligned)

# We are expecting two sequences.
if len(SeqDict) != 2:
  print 'Expected 2 sequences;', DataFile, 'contains', str(len(SeqDict)) +\
  '.\nQuitting.'
  exit(1)

SeqNames = SeqDict.keys()
Seqs     = SeqDict.values()

# If the two sequences are the same length, there is no shuffling to do. Print
# them as they are.
if len(Seqs[0]) == len(Seqs[1]):
  for SeqName, seq in SeqDict.items():
    print '>'+SeqName
    print seq
  exit(0)

# Find which of the two sequences is the shorter one
if len(SeqDict[SeqNames[0]]) < len(SeqDict[SeqNames[1]]):
  ShorterSeqName, ShorterSeq, LongerSeqName, LongerSeq = \
  SeqNames[0], SeqDict[SeqNames[0]], SeqNames[1], SeqDict[SeqNames[1]]
else:
  ShorterSeqName, ShorterSeq, LongerSeqName, LongerSeq = \
  SeqNames[1], SeqDict[SeqNames[1]], SeqNames[0], SeqDict[SeqNames[0]]

LenLongerSeq = len(LongerSeq)
deficit = LenLongerSeq - len(ShorterSeq)

# When aligning the sequences, we want upper and lower case versions of the same
# letter to still count as a match. Create copies which are all upper case.
ShorterSeqUpper = ShorterSeq.upper()
LongerSeqUpper  = LongerSeq.upper()

# Try all ways of moving the shorter sequence along the longer one and keep
# track of how many matches there are in each case.
NumbersOfMatches = []
for LeftPadding in range(0,deficit+1):
  PaddedShorterSeq = '-'*LeftPadding + ShorterSeqUpper + \
  '-'*(deficit - LeftPadding)
  MatchingPositions = 0
  for position in range(0,LenLongerSeq):
    if LongerSeqUpper[position] == PaddedShorterSeq[position]:
      MatchingPositions += 1
  # NB the position in NumbersOfMatches we are about to fill (beginning with the
  # zeroth element) equals the number of gaps we padded the left hand side with.
  NumbersOfMatches.append(MatchingPositions)

# Imagine that NumbersOfMatches = [3,100,6]. That would mean that adding one gap
# to the left and one to the right gives the maximum number of matches (100) 
# between the two sequences. The index of the maximum element in 
# NumbersOfMatches is the number of gaps to pad the left hand side with.
BestLeftPadding = NumbersOfMatches.index(max(NumbersOfMatches))

print '>'+ShorterSeqName
print '-'*BestLeftPadding + ShorterSeq + '-'*(deficit - BestLeftPadding)
print  '>'+LongerSeqName
print LongerSeq
