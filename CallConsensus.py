#!/usr/bin/env python2
from __future__ import print_function
#
## A script to illustrate reading in data files with python,
## by Chris Wymant c.wymant@imperial.ac.uk
##
## Overview:
ExplanatoryMessage = '''
This script interprets a base frequency file of the format produced
by AnalysePileup.py and calls the consensus sequence.
Usage: call it from the command line with the first argument being a base
frequency file, the second argument being the minimum number of reads in
agreement with each other before a base is called, and the third argument
being the minimum number of reads in agreement with each other before upper
case is used for the base present.'''
##
################################################################################
# USER INPUT
# The minimum fraction of reads at a position before we call that base (or those
# bases, when one base alone does not reach that threshold fraction).
MinFraction = 0.6
# The character we call when coverage is below the specified minumum.
CharForLowCoverage = '?'
# Gap characters in the base frequency file.
GapChars = '-.?'
# On the line in the base frequency file containing the reference name, what
# string comes before the reference name itself?
RefNamePrefix = '>Reference_'
# Wrap the sequence in the output fasta file to this length per line
FastaSeqLineLength=50
PrintRefToo=True
################################################################################


# Import some modules we'll need.
import os.path
import sys
from AuxiliaryFunctions import ReverseIUPACdict2
import argparse

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BaseFreqFile', type=File)
parser.add_argument('MinCoverage', type=int)
parser.add_argument('MinCovForUpper', type=int)
parser.add_argument('-A', '--require-read-agreement', action='store_true', \
help='Interpret the minumum coverage as the minimum number of reads supporting'\
+"the most common base here (by default, it's just the minimum number of"+\
"reads).")
args = parser.parse_args()

BaseFreqFile = args.BaseFreqFile
MinCoverage = args.MinCoverage
MinCovForUpper = args.MinCovForUpper

# Check that MinCoverage and MinCovForUpper are positive integers, the 
# latter not smaller than the former.
if MinCoverage < 1:
  print('The specified MinumumCoverageToCallBase of', MinCoverage, \
  'is less than 1. Quitting.', file=sys.stderr)
  exit(1)
if MinCovForUpper < MinCoverage:
  print('The specified MinumumCoverageToUseUpperCase of', MinCoverage, \
  'is less than the specified MinumumCoverageToCallBase. Quitting.', \
  file=sys.stderr)
  exit(1)


# Read in the base frequency file
consensus = ''
RefSeq = ''
FoundRefName = False
RefNamePrefixLength = len(RefNamePrefix)
with open(BaseFreqFile, 'r') as f:

  # Loop through all lines in the file
  for line in f:

    if not FoundRefName and len(line) >= RefNamePrefixLength and \
    line[0:RefNamePrefixLength] == RefNamePrefix:
      FoundRefName = True
      RefName = line[RefNamePrefixLength:].rstrip()
      continue

    # Split up the line into fields separated by commas
    fields = line.split(',')

    # Try to get the columns we want
    try:
      RefBase = fields[2]
      coverage = fields[3]
      BasesHere = fields[5]
      freqs = fields[6:]
    except IndexError:
      # Not enough fields! Print the line and exit with an error.
      print('Line\n', line, 'in', BaseFreqFile, 'contains only', \
      len(fields), 'fields; expected at least 7.\nQuitting', file=sys.stderr)
      exit(1)

    # Try to convert coverage to an int. Complain if it's negative. If it's
    # below threshold, call the character reserved for that case.
    try:
      coverage = int(coverage)
    except ValueError:
      print('Could not understand the coverage in line\n', line, \
      'in', BaseFreqFile, 'as an integer.\nQuitting', file=sys.stderr)
      exit(1)
    if coverage < 0:
      print('Negative coverage on line', line, 'in', BaseFreqFile+\
      '.\nQuitting', file=sys.stderr)
      exit(1)

    # Append the reference base
    if len(RefBase) == 1:
      RefSeq += RefBase
    else:
      print('The reference base on line', line, 'in', BaseFreqFile, 'is', \
      RefBase+'. One character only was expected.\nQuitting.', file=sys.stderr)
      exit(1)

    # Call the appropriate character if coverage is below the threshold.
    # We first test whether coverage < MinCoverage; if it is, the majority
    # coverage is certainly < MinCoverage. We test this first in case the
    # coverage is zero - then we don't want to try to interpret the frequencies.
    if coverage < MinCoverage:
      consensus += CharForLowCoverage
      continue

    # Try to convert the frequencies to floats.
    NumBasesHere = len(freqs)
    try:
      for i in range(0,NumBasesHere):
        freqs[i] = float(freqs[i])
    except ValueError:
      print('Could not understand a frequency in line\n', line, \
      'in', BaseFreqFile, 'as a float.\nQuitting', file=sys.stderr)
      exit(1)

    # Sanity checks on the frequency values
    if NumBasesHere != len(BasesHere):
      print('The number of bases and frequencies differ on line', line, 'in', \
      BaseFreqFile+'.\nQuitting', file=sys.stderr)
      exit(1)
    if freqs != sorted(freqs, reverse=True):
      print('The frequencies on line', line, 'in', BaseFreqFile,\
      'do not decrease in size.\nQuitting', file=sys.stderr)
      exit(1)
    if not 0.999 < sum(freqs) < 1.001:
      print('The frequencies on line', line, 'in', BaseFreqFile,\
      'do not sum to 1.\nQuitting', file=sys.stderr)
      exit(1)

    # Call the appropriate character if coverage is below the threshold.
    MajorityCoverage = coverage * freqs[0]
    if args.require_read_agreement and MajorityCoverage < MinCoverage-0.5:
      consensus += CharForLowCoverage
      continue

    # Keep including bases until the threshold fraction is reached.
    NumBasesToCallHere = 1
    TheirFraction = freqs[0]
    while TheirFraction < MinFraction:
      try:
        TheirFraction += freqs[NumBasesToCallHere]
      except IndexError:
        print('Unexpected error: the sum of all freqencies on line', line, \
        'in', BaseFreqFile, 'is less than the user-specified threshold', \
        str(MinFraction)+'.\nQuitting', file=sys.stderr)
        exit(1)
      NumBasesToCallHere += 1

    if NumBasesToCallHere == 1:
      BaseHere = BasesHere[0]
    else:
      BasesToCallHere = BasesHere[:NumBasesToCallHere]
      # If a gap is one of the things most common at this position, call an 'N';
      # otherwise, find the ambiguity code for this set of bases.
      GapHere = False
      for GapChar in GapChars:
        if GapChar in BasesToCallHere:
          GapHere = True
          break
      if GapHere:
        BaseHere = 'N'
      else:  
        BasesToCallHere = ''.join(sorted(BasesToCallHere))
        try:
          AmbiguityCode = ReverseIUPACdict2[BasesToCallHere]
        except KeyError:
          print('Unexpected set of bases', BasesToCallHere, 'at line',line, \
          'in', BaseFreqFile, 'not found amonst those for which we have', \
          'ambiguity codes, namely:', ' '.join(ReverseIUPACdict2.keys()) +\
          '\nQuitting.', file=sys.stderr)
          exit(1)
        BaseHere = AmbiguityCode
    if MajorityCoverage < MinCovForUpper-0.5:
      BaseHere = BaseHere.lower()
    else:
      BaseHere = BaseHere.upper()
    consensus += BaseHere

if not FoundRefName:
  print('Did not find a line beginning', RefNamePrefix, 'in', BaseFreqFile +\
  '.\nQuitting.', file=sys.stderr)
  exit(1)

# Thanks Stackoverflow:
def insert_newlines(string, every=FastaSeqLineLength):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

# As a name for the sequence take the file name minus extension and path.
SeqName = os.path.basename(BaseFreqFile).rsplit('.',1)[0]
print('>'+SeqName)
print(insert_newlines(consensus))
if PrintRefToo:
  print('>'+RefName)
  print(insert_newlines(RefSeq))
