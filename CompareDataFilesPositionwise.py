#!/usr/bin/env python2

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script should be called with many (or at least 2)
data files as arguments. Each data file should consist of whitespace-delimited
columns, one of which is a monotonically increasing positive integer, and
another being a (float) datum of interest. Firstly, missing lines in a file -
i.e. if the integer starts at a value greater than the MinInteger specified,
and/or where it increases by more than 1 from one line to the next - have the
value of the float set to a user specified value. Then, the values associated
with the same value of the integer in different files are compared -
user-specified percentiles are calculated - and these are plotted as a function
of the integer.
An example: you have some files, each of which contains a set of integers; you
pipe each file separately into
$ sort -n MyFile.dat | uniq -c > MyFile_ValueCounts.dat
giving files where the first column is an (increasing) integer and the second is
a count of how many times that integer appears. This script gives the
distribution over those integers of the distribution of integer counts within
each file.
'''

################################################################################
# USER INPUT
percentiles=[25,50,75]
yaxislabel='fraction of inserts of size $i$, per million'
xaxislabel='insert size, $i$'
title='$x$ and $y$ axes both truncated for clarity'
xlabelfontsize=16
ylabelfontsize=18
ticksize=16
linewidth=2
IncludeLegend=True
LegendLocation='upper right'
UseLogScaleYaxis=False
OutfileName='/home/chris/Dropbox (Infectious Disease)/chris/BEEHIVE/InsertSizes_BEEHIVE_ReadLength300_foo.pdf'

# Set any non-integer value to have automatic limits
xmin = 220
xmax = 310
# Set any non-numerical value to have automatic limits
ymin = 0
ymax = 25
################################################################################

percentiles = sorted(percentiles, reverse=True)

import os
import string
import argparse
import collections
import numpy
import matplotlib.pyplot as plt

# Set up the arguments for this script
ExplanatoryMessage = string.replace(ExplanatoryMessage, '\n', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage, \
epilog = 'Options relating the plot output are contained inside the code.')
parser.add_argument('IntegerColumnNumber', type=int,\
help='Which column in the data files is the increasing integer? e.g. a value '\
+'of 1 would indicate the first column.')
parser.add_argument('DatumColumnNumber', type=int, \
help='Which column in the data files is the datum?')
parser.add_argument('MinInteger', type=int,\
help='The minimum value of the integer that we should encounter.')
parser.add_argument('MissingDatumValue', type=float, \
help='The value of the datum to be used for lines that are missing.')
parser.add_argument('DataFile', nargs='+')
parser.add_argument('-D', '--delimiter', default=None, help='''Used to specify a
different delimiter for the data files (the default is whitespace).''')
args = parser.parse_args()

# Relabelling for brevity
IntegerColumn = args.IntegerColumnNumber
DatumColumn = args.DatumColumnNumber
MinInt = args.MinInteger
MissingDatumValue = args.MissingDatumValue
DataFiles = args.DataFile

# Check that each data file was specified only once.
CounterObject = collections.Counter(DataFiles)
DuplicatedArgs = [arg for arg in CounterObject if CounterObject[arg] > 1]
if len(DuplicatedArgs) != 0:
  print 'The following arguments were duplicated:', ' '.join(DuplicatedArgs)
  print 'All arguments should be unique. Quitting.'
  exit(1)

# Check that the data files exist and are files
for DataFile in DataFiles:
  if not os.path.isfile(DataFile):
    print DataFile, 'does not exist or is not a file.'
    exit(1)

AllDataDict = {}

# Iterate through all data files
for DataFile in DataFiles:
  with open(DataFile, 'r') as f:
    DataForThisFile = []
    GreatestIntSoFar = MinInt-1
    for line in f:

      EndingForErrorMessage = 'on line\n' + line +'in '+ DataFile+'.\nQuitting'

      # split line into fields
      fields = line.split(args.delimiter)

      # Try to get the columns we want
      try:
        integer = fields[IntegerColumn-1]
        datum = fields[DatumColumn-1]
      except IndexError:
        print 'Not enough fields '+EndingForErrorMessage
        exit(1)

      # Try to convert to an int and float
      try:
        integer = int(integer)
      except ValueError:
        print 'Could not understand the integer ' +EndingForErrorMessage
        exit(1)
      try:
        datum = float(datum)
      except ValueError:
        print 'Could not understand the datum as a float '+EndingForErrorMessage
        exit(1)

      # Check the integer is increasing
      if integer < MinInt:
        print 'Encountered integer '+str(integer) + ' ' +EndingForErrorMessage
        exit(1)
      if integer <= GreatestIntSoFar:
        print 'Encountered integer', integer, 'after', str(LastRefPosition)+\
        ' ' +EndingForErrorMessage
        exit(1)

      # Fill in missing lines i.e. missing values of the integer
      if integer != GreatestIntSoFar+1:
        DataForThisFile += [MissingDatumValue] * (integer -GreatestIntSoFar -1)

      # Record this line
      DataForThisFile.append(datum)
      GreatestIntSoFar = integer

  AllDataDict[DataFile] = DataForThisFile

# Calculate the desired percentiles (over the distribution of values found in
# the files) for each integer value.
LengthOfLongestList = max([len(List) for List in AllDataDict.values()])
ListsOfPercentiles = [[] for percentile in percentiles]
for ListPosition in range(0,LengthOfLongestList):
  DataHere = []
  for List in AllDataDict.values():
    try:
      datum = List[ListPosition]
    except IndexError:
      datum = MissingDatumValue
    DataHere.append(datum)
  for i,percentile in enumerate(percentiles):
    ListsOfPercentiles[i].append(\
    numpy.percentile(DataHere,percentile))

#for List in AllDataDict.values():
#  print numpy.log(List[165]/List[164]), numpy.log(List[165]/List[166])


# Allow user-specified x- and y-ranges if appropriate
try:
  xmin, xmax = int(xmin), int(xmax)
except ValueError:
  pass
else:
  xmin = max(xmin, MinInt)
  xmax = min(xmax, MinInt+LengthOfLongestList)
  plt.xlim((xmin,xmax))
try:
  ymin, ymax = float(ymin), float(ymax)
except ValueError:
  pass
else:
  plt.ylim((ymin,ymax))

# Plot!
xValues=range(MinInt,MinInt+LengthOfLongestList)
for i,ListOfPercentiles in enumerate(ListsOfPercentiles):
  percentile=percentiles[i]
  plt.plot(xValues,[value*10**6 for value in ListOfPercentiles],linewidth=linewidth,\
  label=str(percentile)+'$^{th}$ percentile')

if IncludeLegend:
  plt.legend(loc=LegendLocation)
plt.xlabel(xaxislabel, fontsize=xlabelfontsize)
plt.ylabel(yaxislabel, fontsize=ylabelfontsize)
fig = plt.figure(1)
plot = fig.add_subplot(111)
if UseLogScaleYaxis:
  plot.set_yscale('log')
plot.tick_params(axis='both', which='major', labelsize=ticksize)
plot.tick_params(axis='both', which='minor', labelsize=ticksize)
plt.title(title)
plt.savefig(OutfileName)


