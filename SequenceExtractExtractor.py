#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview: run from the command line with one aligned fasta file as an 
## argument, and start & end positions specified internally below, this script
## prints only that part of the alignment between the start & end positions.
##

import os.path, sys
from Bio import SeqIO
from Bio import Seq


# Check this file is called from the command line with one argument
if len(sys.argv) != 4:
  print('Incorrect number of arguments given.\nUsage:\n', sys.argv[0], \
  'NameOfYourFastaFile.fasta StartPosition EndPosition', file=sys.stderr)
  exit(1)
FastaFile = sys.argv[1]
start=int(sys.argv[2])
end=int(sys.argv[3])

# Check that the argument exists and is a file
if not os.path.isfile(FastaFile):
  print(FastaFile, 'does not exist or is not a file.', file=sys.stderr)
  exit(1)

# Sanity checks on the start and end.
if start >= end:
  print('StartPosition should be less than EndPosition. Quitting.', \
  file=sys.stderr)
  exit(1)
if start < 1:
  print('StartPosition should be greater than zero. Quitting.', file=sys.stderr)
  exit(1)

seqs = []

for seq in SeqIO.parse(open(FastaFile),'fasta'):
  try:
    CutSeq = str(seq.seq)[start-1:end]
  except IndexError:
    print('The coordinates', start, ', ', end, 'go beyond the edges of', \
    seq.id+'.Quitting', file=sys.stderr)
    exit(1)
  seq.seq = Seq.Seq(CutSeq)
  seqs.append(seq)
NumSeqs = len(seqs)
if NumSeqs == 0:
  print('No sequences found in', FastaFile +'. Quitting.', file=sys.stderr)
  exit(1)

SeqIO.write(seqs, sys.stdout, "fasta")

  

