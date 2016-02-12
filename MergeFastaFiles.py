#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script takes one or more fasta files as arguments,
and reads the sequences contained therein. Sequences encountered a second time -
i.e. with the same name - are checked to be identical. All sequences are then
printed to stdout, in the same order in which they appeared in the files, with
duplicates skipped.'''

import argparse
import os
import sys
from Bio import SeqIO
from collections import OrderedDict

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File, nargs='+')
parser.add_argument('-C', '--case-sensitive', action='store_true', \
help='when comparing sequences that appear multiple times in the input files,'+\
'this switch causes differences in case to be considered as genuine'+\
'differences, thus causing the program to exit with an error.'+\
' (By default, case is ignored.)')
args = parser.parse_args()


AllUniqueSeqs = OrderedDict()
for FastaFile in args.FastaFile:
  for seq in SeqIO.parse(open(FastaFile),'fasta'):
    if seq.id in AllUniqueSeqs:
      ThisOccurence = str(seq.seq)
      FormerOccurrence = str(AllUniqueSeqs[seq.id].seq)
      if args.case_sensitive:
        DuplicateIsDifferent = (ThisOccurence != FormerOccurrence)
      else:
        DuplicateIsDifferent = \
        (ThisOccurence.upper() != FormerOccurrence.upper())
      if DuplicateIsDifferent:
        print('Sequence', seq.id, 'occurs multiple times in the input fasta'+\
        ' files and not all occurences are identical.\nQuitting.', \
        file=sys.stderr)
        exit(1)
    AllUniqueSeqs[seq.id] = seq

SeqIO.write(AllUniqueSeqs.values(), sys.stdout, "fasta")


