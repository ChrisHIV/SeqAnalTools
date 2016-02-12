#!/usr/bin/env python2

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview: this script, called from the command line with a fasta file as an
## argument, calculates the entropy of the (nominally aligned) set of sequences
## contained therein. IUPAC ambiguity codes are interpreted as counting for half
## or one third of each of their associated letters. The characters '-', 'n' and
## 'N' are deemed to be equivalent and are treated as a fifth distinct letter.

import os.path, sys, collections, numpy
from AuxiliaryFunctions import ReadSequencesFromFile, IUPACdict, \
ReverseIUPACdict, InterpretIUPAC

GapChars = ['-','.','?']

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

# Read in all the sequences as a dictionary
AllSequences, SequenceLength = ReadSequencesFromFile(DataFile)

SequenceNames = AllSequences.keys()
sequences     = AllSequences.values()
NumSequences  = len(AllSequences)

DataFromAllPositions = []
SetOfAllBasesEncountered = []

# For brevity
acgt = ['A','C','G','T']


AllExpectedBases = acgt + IUPACdict.keys() + GapChars



TotalEntropy = 0    

# Consider each position (column) in the alignment in turn:
for i in range(0,SequenceLength):

  # A list of all bases (including '-' etc.) occurring at this position
  AllBasesHere = [sequence[i] for sequence in sequences]

  # A dictionary counting how many times each base occurs
  NumbersOfEachBaseHere = collections.Counter(AllBasesHere)

  # Quit if we encounter unexpected bases
  UnexpectedBases = [base for base in NumbersOfEachBaseHere.keys() if not \
  base in AllExpectedBases]
  if len(UnexpectedBases) > 0:
    print 'Encountered the following unexpected bases at position', i+1, \
    'in the alignment:', ' '.join(UnexpectedBases), '\nQuitting.'
    exit(1)

  # Add together all ways of expressing a gap
  NumGaps = 0
  for GapChar in GapChars:
    if GapChar in NumbersOfEachBaseHere:
      NumGaps += NumbersOfEachBaseHere[GapChar]

  # Divide the values for IUPAC ambiguity codes equally between their
  # corresponding 'normal' letters (a,c,g,t). This forgets about gaps.
  NumbersOfEachBaseHere = InterpretIUPAC(NumbersOfEachBaseHere)

  # Remember the gaps.
  if NumGaps != 0:
    NumbersOfEachBaseHere['-'] = NumGaps

  EntropyHere = 0
  for BaseCount in NumbersOfEachBaseHere.values():
    # NB BaseCount > 0 because it counts things that occured, so the log is fine
    p = float(BaseCount) / NumSequences 
    EntropyHere += -p * numpy.log(p)

  TotalEntropy += EntropyHere

print TotalEntropy
