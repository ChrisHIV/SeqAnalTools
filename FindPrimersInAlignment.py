#!/usr/bin/env python2

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251 
##
## Overview: TODO
## Usage: call this script from the command line with the first argument being
## the alignment file, the second argument being the name of the chosen
## reference, and at least one primer: each primer coming after either -S
## (indicating that you want the position of the start of the primer) or -E
## (you want the position of the end of the primer).
## e.g. using the following options (without the line break of course)
## -E AAAATGATAGGRGGAATTGGAGG -S GGGAAGTGAYATAGCWGGAAC -E GAYTATGGAAAACAGATGGCAG
## -S TTAAAAGAAAAGGGGGGATTGGG -E CGCTGACGGTACARGCCA -S CCTATGGCAGGAAGAAGCG
## will return the start/end locations of these six primers, which are the
## primers of Gall et al. doi:10.1128/JCM.01516-12  
## Positions before a sequence begins are given as 1; positions after a sequence
## ends are given as the final position in that sequence.
##
################################################################################
## USER INPUT
# Characters that indicates a gap (missing base)
GapChars = '-.?'
################################################################################

# Import what's needed
from optparse import OptionParser
import os.path, collections
from AuxiliaryFunctions import ReadSequencesFromFile, IUPACdict, BaseMatch

# Define the arguments and options
parser = OptionParser()
parser.add_option("-S", "--start", action="append", type="str")
parser.add_option("-E", "--end", action="append", type="str")
parser.add_option("-A", "--alignment-coords", dest='AlignmentCoords', \
action="store_true", default=False, help='print out the primer positions just'+\
' as coordinates with respect to the alignment')
(options, args) = parser.parse_args()

# Check this file is called from the command line with the correct number of
# arguments, and that the specified file(s) exist.
if len(args) != 2:
  print 'Two arguments are required: firstly the alignment file, and secondly',\
  'the chosen reference therein.\nQuitting.'
  exit(1)
AlignmentFile = args[0]
ChosenRef     = args[1]
if not os.path.isfile(AlignmentFile):
  print AlignmentFile, 'does not exist or is not a file.'
  exit(1)

# If no start primers were specified or no end primers, make empty lists for
# convenience.
StartPrimers = options.start
EndPrimers = options.end
if StartPrimers == None:
  StartPrimers = []
if EndPrimers == None:
  EndPrimers = []

# Check that at least one primer was specified
if len(StartPrimers) == 0 and len(EndPrimers) == 0:
  print 'At least one primer is required, specified using the -S or -E',\
  'options.\nQuitting.'
  exit(1)

# Convert the primers to upper case, check they don't contain gaps (that would
# be weird) and check they are unique (though they may appear in both lists).
for PrimerList in [StartPrimers,EndPrimers]:
  for i in range(0,len(PrimerList)):
    PrimerList[i] = PrimerList[i].upper()
    for GapChar in GapChars:
      if GapChar in PrimerList[i]:
        print 'Primer', PrimerList[i], 'contains a gap. This is unexpected.'+\
        '\nQuitting.'
        exit(1)
  CounterObject = collections.Counter(PrimerList)
  DuplicatedPrimers = [i for i in CounterObject if CounterObject[i]>1]
  if len(DuplicatedPrimers) != 0:
    for DuplicatedPrimer in DuplicatedPrimers:
      print 'Primer', DuplicatedPrimer, 'was specified twice with the same',\
      'option.'
    print 'Quitting.'
    exit(1)

# Read in the sequences from the alignment file (into a dictionary)
SeqDict, AlignmentLength = ReadSequencesFromFile(AlignmentFile)

# Check the chosen reference is in the alignment
if not ChosenRef in SeqDict:
  print 'Could not find', ChosenRef, 'in', AlignmentFile+'.\nQuitting.'
  exit(1)
ChosenRefSeq = SeqDict[ChosenRef]

# Define the set of unique primers, i.e. StartPrimers+EndPrimers but not 
# double counting those that appear in both. Record their lengths in a dict.
AllUniquePrimers = StartPrimers + \
[primer for primer in EndPrimers if not primer in StartPrimers]
NumUniquePrimers = len(AllUniquePrimers)
PrimerLengths = {primer : len(primer) for primer in AllUniquePrimers} 

# Finds the position in the alignment, for each primer, after
# which we should stop checking for primer-reference matches because the number
# of reference bases (not including gaps) to the right is less than the length
# of the primer.
PrimerLastChancesForMatch = {}
NumRefBasesSoFar = 0
AllPrimersShorterThanRef = False
for PositionMin1 in range(AlignmentLength-1,-1,-1):
  if not ChosenRefSeq[PositionMin1] in GapChars:
    NumRefBasesSoFar += 1
    for primer,length in PrimerLengths.items():
      if NumRefBasesSoFar == length:
        PrimerLastChancesForMatch[primer] = PositionMin1+1
    if len(PrimerLastChancesForMatch) == NumUniquePrimers:
      AllPrimersShorterThanRef = True
      break
# Any primer for which we worked our way through the whole reference without
# finding enough bases to get to the length of the primer must be longer than
# the reference.
if not AllPrimersShorterThanRef:
  OverlyLongPrimers = [primer for primer in AllUniquePrimers if not primer in \
  PrimerLastChancesForMatch]
  print 'The following primers are longer than', ChosenRef+':', \
  ' '.join(OverlyLongPrimers) +'\nQuitting.'
  exit(1)

# Find the primers in the chosen reference sequence, converted to upper case.
StartPrimerPositions = {}
EndPrimerPositions = {}
ChosenRefSeq = ChosenRefSeq.upper()
for PositionMin1,base in enumerate(ChosenRefSeq):

  # Skip gaps
  if base in GapChars:
    continue

  # For each primer, work forward through both the primer and the chosen ref
  # while the two sequences agree, stopping as soon as they don't, ignoring
  # gaps in the reference.
  # When it equals the primer length, we have found that primer: we record the
  # primer position (start and/or end) in the alignment.
  for primer,length in PrimerLengths.items():
    if PositionMin1 > PrimerLastChancesForMatch[primer]-1:
      continue
    matches = 0
    StepsForward = 0
    while BaseMatch(ChosenRefSeq[PositionMin1+StepsForward],primer[matches]) \
    or ChosenRefSeq[PositionMin1+StepsForward] in GapChars:
      if not ChosenRefSeq[PositionMin1+StepsForward] in GapChars:
        matches += 1
        if matches == length:
          # Found the primer!
          if primer in StartPrimerPositions or primer in EndPrimerPositions:
            print 'Encountered primer', primer, 'a second time in',ChosenRef+ \
            '. Primers should be unique.\nQuitting.'
            exit(1)
          if primer in StartPrimers:
            StartPrimerPositions[primer] = PositionMin1 +1
          if primer in EndPrimers:
            EndPrimerPositions[primer] = PositionMin1 +StepsForward+1
          break
      StepsForward += 1
        
# Check that all primers were found
MissingPrimers = [primer for primer in AllUniquePrimers if \
(not primer in StartPrimerPositions) and (not primer in EndPrimerPositions)]
if len(MissingPrimers) != 0:
  print 'Unable to find', ' or '.join(MissingPrimers), 'in', ChosenRef+\
  '.\nQuitting.'
  exit(1)

# Merge the start and end positions into a single sorted list, each value being
# coupled to its primer name with 'start_of_' or 'end_of_' prepended.
SortedList = [['start_of_'+primer,value] for primer,value in \
StartPrimerPositions.items()] + [['end_of_'+primer,value] for primer,value in \
EndPrimerPositions.items()]
SortedList = sorted(SortedList, key=lambda item: item[1])

if options.AlignmentCoords:
  print ' '.join(str(value) for key,value in SortedList)
  exit(0)

# Now convert those primer positions, which were with respect to the alignment,
# into positions with respect to each reference.
# TODO: update comment
PositionsDict = {}
for SeqName, seq in SeqDict.items():
  PositionsWRTseq = []
  LastPositionWRTalignment = 0
  PositionWRTseq = 0
  for [name,PrimerPosition] in SortedList:
    for base in seq[LastPositionWRTalignment:PrimerPosition]:
      if not base in GapChars:
        PositionWRTseq += 1
    PositionsWRTseq.append(PositionWRTseq)
    LastPositionWRTalignment = PrimerPosition

  # Replace any zeroes by ones (i.e. map positions off to the left onto the
  # first position for this sequence).
  for i in range(0,len(PositionsWRTseq)):
    if PositionsWRTseq[i] == 0:
      PositionsWRTseq[i] = 1
    else:
      break

  PositionsDict[SeqName] = PositionsWRTseq

# Print output, with an explanatory header line beginning with a hash.
print '# name ', '  '.join([str(item[0]) for item in SortedList])
for SeqName in sorted(PositionsDict.keys()):
  print SeqName, ' '.join(map(str,PositionsDict[SeqName]))


