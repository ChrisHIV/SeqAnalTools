#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script converts sequence data between formats.
Output is printed to stdout.'''

import argparse
import os
import sys
from Bio import SeqIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('InputSeqFile', type=File)
parser.add_argument('InputFormat')
parser.add_argument('OutputFormat')
args = parser.parse_args()

SeqIO.convert(args.InputSeqFile, args.InputFormat, sys.stdout, args.OutputFormat)