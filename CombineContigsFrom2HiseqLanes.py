#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script merges all contigs assembled from two Hiseq
lanes, each labelled with a Sanger Institute lane ID, into a single new fasta
file with uniquely labelled contigs (replacing _1_ and _2_ with _3_, and
appending _r for contigs coming from the _2_ lane).'''

import argparse
import os
import sys
import collections
from Bio import SeqIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile1', type=File)
parser.add_argument('FastaFile2', type=File)
args = parser.parse_args()

def FastaFileToSeqDict(FastaFile):
  '''Reads a fasta file into a dictionary.'''
  MyDict = collections.OrderedDict()
  for seq in SeqIO.parse(open(FastaFile),'fasta'):
    if seq.id in MyDict:
      print('Sequence', seq.id, 'occurs multiple times in', FastaFile+\
      '\nQuitting.', file=sys.stderr)
      exit(1)
    MyDict[seq.id] = seq
  return MyDict

def occurrences(string, sub):
  '''Counts the number of potentially overlapping occurences of a substring 
  inside a string.'''
  count = start = 0
  while True:
    start = string.find(sub, start) + 1
    if start > 0:
      count+=1
    else:
      return count

def RenameSeqs(SeqDict, WhichFasta):
  'Replaces _x_ by _3_ in seq names, where x is the second argument.'
  MyDict = collections.OrderedDict()
  substring = '_'+str(WhichFasta)+'_'
  for SeqName, seq in SeqDict.items():
    if occurrences(seq.id, substring) != 1:
      print('Sequence name', seq.id, 'contains more than one occurence of the'+\
      ' string "'+substring+'". Quitting, as it is ambiguous how to rename.',\
      file=sys.stderr)
      exit(1)
    seq.id = seq.id.replace(substring, '_3_')
    seq.description = ''
    MyDict[seq.id] = seq
  return MyDict


File1dict = FastaFileToSeqDict(args.FastaFile1)
File2dict = FastaFileToSeqDict(args.FastaFile2)
File1dict = RenameSeqs(File1dict, 1)
File2dict = RenameSeqs(File2dict, 2)

for seq in File2dict.values():
  seq.id += '_r'
  if seq.id in File1dict:
    print('Unexpected clash of names between files for', seq.id, '\nQuitting.',\
    file=sys.stderr)
    exit(1) 
  File1dict[seq.id] = seq

SeqIO.write(File1dict.values(), sys.stdout, "fasta")
