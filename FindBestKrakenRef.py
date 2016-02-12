#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script reads a Kraken report file and returns the
name of the taxon that has the most reads rooted under it, only looking inside a
specified (larger) taxon. If no reads are rooted to any taxon inside the
specified larger taxon, N/A is printed.'''

import sys
import os
import argparse

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('KrakenReportFile', type=File)
parser.add_argument('TaxonName', help='The scientific name'\
+' of the taxon you want to consider (i.e. hits outside of this taxon will be'+\
'ignored). If it contains whitespace, you must surround it with quotation '+\
'marks.')
parser.add_argument('-D', '--dict', type=File, help='A csv file in which the '+\
'first column is the name of a taxon used by Kraken, and the second is '+\
'something you want to rename it to.')
args = parser.parse_args()

if args.dict != None:
  import csv

PossibleRankCodes = ['U', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S', '-']

DesiredTaxonNameIndent = -1
results = {}
with open(args.KrakenReportFile, 'r') as f:
  for LineNumberMin1,line in enumerate(f):

    # Split the line into fields the first 5 times white space is encountered
    # (not thereafter, as the taxon name could contain whitespace).
    fields = line.split(None, 5)
    if len(fields) < 6:
      print('Not enough fields on line', str(LineNumberMin1+1)+'. Quitting.', \
      file=sys.stderr)
      exit(1)

    # Check the fields look like they supposed to.
    try:
      percentage = float(fields[0])
      NumReadsRootedHere = int(fields[1])
      NumReadsDirectlyHere = int(fields[2])
      RankCode = fields[3]
      assert RankCode in PossibleRankCodes
      NCBItaxonomyID = int(fields[4])
    except (ValueError,AssertionError):
      print('Unexpected field format on line', str(LineNumberMin1+1) +\
      '. Quitting.', file=sys.stderr)
      exit(1)

    # Find the position of the taxon name in the line.
    TaxonName = fields[5].rstrip()
    TaxonNameIndent = line.find(TaxonName)

    # If this taxon name is the desired one, record the position of the taxon
    # name. If we've already recorded it, that means we've encountered the 
    # desired taxon name a second time.
    if TaxonName == args.TaxonName:
      if DesiredTaxonNameIndent > -1:
        print('Encountered taxon name', TaxonName, 'a second time on line', \
        str(LineNumberMin1+1)+'. Unexpected. Quitting.', file=sys.stderr)
        exit(1)
      DesiredTaxonNameIndent = TaxonNameIndent
      continue

    # If the desired taxon name has been found:
    if DesiredTaxonNameIndent > -1:

      # If the taxon name is not more deeply indented than the desired taxon
      # name, we're not inside the taxon any more.
      if TaxonNameIndent <= DesiredTaxonNameIndent:
        break

      # We're inside the desired taxon! Check the name hasn't already been
      # encountered, and record.
      if TaxonName in results:
        print('Encountered taxon name', TaxonName, 'a second time on line', \
        str(LineNumberMin1+1)+'. Unexpected. Quitting.', file=sys.stderr)
        exit(1)
      results[TaxonName] = NumReadsRootedHere

# Check we found something.
'''
if DesiredTaxonNameIndent == -1:
  print('Taxon name', args.TaxonName, 'was not encountered in', \
  args.KrakenReportFile +'. Quitting.', file=sys.stderr)
  exit(1)
if len(results) == 0:
  print('Nothing was found inside taxon', args.TaxonName +'. Quitting.', \
  file=sys.stderr)
  exit(1)
'''
if DesiredTaxonNameIndent == -1 or len(results) == 0:
  print('N/A')
  exit(0)

# The ref with the highest number of reads rooted there.
BestRef = max(results, key=results.get)
if args.dict == None:
  print(BestRef)

else:
  # Use the renaming csv file.
  BestRefRenamed = None
  with open(args.dict, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for LineNumberMin1, fields in enumerate(reader):
      if len(fields) != 2:
        print('Line', LineNumberMin1+1, 'in', args.dict, \
        'does not have two fields. Quitting.', file=sys.stderr)
        exit(1)
      KrakenName, rename = fields
      if KrakenName == BestRef:
        if BestRefRenamed != None:
          print('Encountered', KrakenName, 'a second time on line', \
          LineNumberMin1+1, 'in', args.dict +'. Quitting.', file=sys.stderr)
          exit(1)
        BestRefRenamed = rename
  print(BestRefRenamed)

