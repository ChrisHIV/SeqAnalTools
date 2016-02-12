#!/usr/bin/env python3

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview: this script calculates the linkage disequilibrium (quantified by D,
## Dprime and r - see http://en.wikipedia.org/wiki/Linkage_disequilibrium )
## between each possible pairing of positions specified by the user.
## Usage: it is called from the command line with the first argument being a bam
## file, and subsequent arguments being positions (with respect to the
## reference) that the user has identified as having enough diversity and
## coverage to be interesting for a linkage disequilibrium calculation. Such
## positions are trivially determined, for example, by running samtools mpileup
## on the bam file, feeding the output to AnalysePileup.py in this directory,
## and using e.g. awk '{ if ($5 < 0.95 && $3 > 19 && $1 != 'NA') print $1; }' 
## to print all positions with at least 5% diversity and at least 20X coverage.
## The fields in the output of this script are:
## the left position in the pair, the right position in the pair, the number of
## reads or read-pairs that span both positions, D, Dprime, r.
## Of the spanning reads, we keep (and count) only those that have either the
## most common or second most common base at both positions.

import pysam, os.path, sys
from collections import Counter

# Check this file is called from the command line with at least three arguments
if len(sys.argv) < 4:
  print('\nUsage',sys.argv[0], 'BamFile VariantPosition1 VariantPosition2',\
  '[VariantPosition3...]')
  print('Quitting.')
  exit(1)

BamFileName = sys.argv[1]
VariantPositions = sys.argv[2:]

# Check that the BamFileName exists and is a file
if not os.path.isfile(BamFileName):
  print(BamFileName, 'does not exist or is not a file.')
  exit(1)

# Make sure all the variant positions are integers.
for i in range(0,len(VariantPositions)):
  try:
    VariantPositions[i] = int(VariantPositions[i])
  except ValueError:
    print('The first argument should be a file and all subsequent arguments',\
    'should be integers.\nQuitting.')
    exit(1)

# Check all the variant positions specified are unique. If not, ignore duplicates.
if len(VariantPositions) > len(set(VariantPositions)):
  print('Warning: some of the variant positions were repeated amongst the',\
  'input. Skipping duplicates.')
  VariantPositions = set(VariantPositions)

# Put the variant positions in increasing order
VariantPositions = sorted(VariantPositions)
NumVarPos = len(VariantPositions)

BamFile = pysam.AlignmentFile(BamFileName, "rb")

# Find the reference in the bam file; there should only be one.
AllReferences = BamFile.references
if len(AllReferences) != 1:
  print('Expected exactly one reference in', BamFileName+'; found',\
  str(len(AllReferences))+'.\nQuitting.')
  exit(1)
TheRef = AllReferences[0]

# Get the length of the reference.
AllReferenceLengths = BamFile.lengths
if len(AllReferenceLengths) != 1:
  print('Pysam error: found one reference but', len(AllReferenceLengths), \
  'reference lengths.\nQuitting.')
  exit(1)
RefLength = AllReferenceLengths[0]

# Check all positions are between 1 and the ref genome length
if VariantPositions[0] < 1 or VariantPositions[-1] > RefLength: 
  print('Variant positions should be specified with respect to the reference,',\
  'and so should be between 1 and the length of the reference, inclusively.'+\
  '\nQuitting.')
  exit(1)


# TESTING: finding maximum fragment sizes
#AllReads = BamFile.fetch(TheRef,0,RefLength)
#LeftReadLeftEdges = {}
#RightReadRightEdges = {}
#for read in AllReads:
#  if read.is_paired:
#    if read.is_read1:
#      LeftEdge = min(read.get_reference_positions())
#      LeftReadLeftEdges[read.query_name] = LeftEdge
#    elif read.is_read2:
#      RightEdge = max(read.get_reference_positions())
#      RightReadRightEdges[read.query_name] = RightEdge
#for read,LeftEdge in LeftReadLeftEdges.items():
#  if read in RightReadRightEdges:
#    RightEdge = RightReadRightEdges[read]
#    if RightEdge - LeftEdge > 1000:
#      print(RightEdge - LeftEdge, read)
#exit(0)


# The nth element in this list will be a dictionary, for the nth variant
# position, in which the indices are the names of the reads and the values are 
# the bases there.
ListOfDictsOfBases = []

# Record all of the bases and their read names at each variant position.
for VarPos in VariantPositions:

  DictOfBases = {}
  reads = BamFile.fetch(TheRef,VarPos-1, VarPos)
  for read in reads:

    ReadName = read.query_name
    if read.is_read1:
      ReadName += '_1'
    elif read.is_read2:
      ReadName += '_2'
    else:
      print('Read', ReadName, 'is neither read 1 nor read 2. Unexpected.'+\
      '\nQuitting.')
      exit(1)

    if ReadName in DictOfBases:
      print('Multiple bases are called',ReadName+'. Unexpected.\nQuitting.')
      exit(1)
    # Where does the VarPos occur within this sequence?
    # (VarPos - read.reference_start - 1 + read.query_alignment_start) almost 
    # always gets it right, but can get it wrong if this sequence has an indel not  
    # present in the consensus. Instead we get the list of all reference positions 
    # mapped to, and find the position in that list of this column. This is two
    # times slower. Exceptionally the column may be missing, if the read spans
    # that position but has a deletion just at that point.
    try:
      WhereVarPos = \
      read.get_reference_positions(full_length=True).index(VarPos-1)
    except ValueError:
      continue
    BaseThere = read.seq[WhereVarPos]

    DictOfBases[ReadName] = BaseThere

  ListOfDictsOfBases.append(DictOfBases)

# TESTING:
#Dict1 = {}; Dict2 = {}
#Dict1['read1_1'] = 'a'; Dict2['read1_1'] = 'g'
#Dict1['read2_1'] = 'a'; Dict2['read2_1'] = 'g'
#Dict1['read3_1'] = 'a'; Dict2['read3_1'] = 'g'
#Dict1['read4_1'] = 'a'; Dict2['read4_1'] = 'g'
#Dict1['read5_1'] = 'a'; Dict2['read5_1'] = 't'
#Dict1['read6_1'] = 'a'; Dict2['read6_1'] = 't'
#Dict1['read7_1'] = 'a'; Dict2['read7_1'] = 't'
#Dict1['read8_1'] = 'a'; Dict2['read8_1'] = 't'
#Dict1['read9_1'] = 'a'; Dict2['read9_1'] = 'g'
#Dict1['read10_1'] = 'c'; Dict2['read10_1'] = 't'
#Dict1['read11_1'] = 'c'; Dict2['read11_1'] = 'a'
#Dict1['foo_1'] = 'Q'; Dict2['spam_1'] = 'W'
#ListOfDictsOfBases = [Dict1, Dict2]
    
# For each VarPos, iterate over all the other VarPos to the right of it. i.e.
# loop through all distinct pairings.
for i,DictOfBases1 in enumerate(ListOfDictsOfBases):
  for j in range(i+1,NumVarPos):

    DictOfBases2 = ListOfDictsOfBases[j]
    LeftVarPos  = VariantPositions[i]
    RightVarPos = VariantPositions[j]

    # TODO: if (RightVarPos - LeftVarPos) large enough: break

    # Loop through all of the reads at the left VarPos. For each, see if that
    # read and/or its mate covers the right VarPos; if so, record the left and
    # right bases.
    LeftBases  = []
    RightBases = []

    for read,LeftBase in DictOfBases1.items():
      if read[-1] == 1:
        MateRead = read[:-1]+'2'
      else:
        MateRead = read[:-1]+'1'
      ThisReadReachesNextVarPos = read in DictOfBases2
      MateReadReachesNextVarPos = MateRead in DictOfBases2

      # If both this read and its mate cover the right VarPos, and they disagree
      # about what base is there, skip! (This situation is ambiguous.)
      if ThisReadReachesNextVarPos and MateReadReachesNextVarPos:
        #print(read,LeftBase, DictOfBases2[read], MateRead, DictOfBases2[MateRead])
        if DictOfBases2[read] != DictOfBases2[MateRead]:
          continue
        else:
          RightBase = DictOfBases2[read]
      elif ThisReadReachesNextVarPos:
        RightBase = DictOfBases2[read]
        #print(read,LeftBase,RightBase)
      elif MateReadReachesNextVarPos:
        RightBase = DictOfBases2[MateRead]
        #print(read,LeftBase,MateRead,RightBase)
      else:
        continue

      
      if LeftBase in ['A','C','G','T'] and RightBase in ['A','C','G','T']:
        LeftBases.append(LeftBase)
        RightBases.append(RightBase)

    # Skip this pairing of VarPos if no reads (or read pairs) spanned the pair
    if len(LeftBases) == 0:
      continue

    # Find the most common two bases at the left VarPos and the right VarPos.
    MostCommon2bases_left  = Counter(LeftBases).most_common(2)
    MostCommon2bases_right = Counter(RightBases).most_common(2)

    # If either site has no minority bases at all (only one base is seen), there
    # can be no linkage disequilibrium. 'count' is the number of spanning reads
    # used to calculate linkage disequilibrium.
    # Use C for consensus, M for minority, for brevity.
    NoLD = False
    if len(MostCommon2bases_left) == 1 or len(MostCommon2bases_right) == 1:
      NoLD = True
      count = len(LeftBases)
    else:

      [(C_base_left, C_count_left), (M_base_left, M_count_left)] = \
      MostCommon2bases_left
      [(C_base_right, C_count_right), (M_base_right, M_count_right)] = \
      MostCommon2bases_right

      # Restrict our attention to the subset of (paired) reads having
      # either C or M at the first site and C or M at the second site (i.e.
      # ignore any cases containing the third-most-common base).
      # In the very unlikely event that there are no such reads, continue.
      CC_count = 0
      CM_count = 0
      MC_count = 0
      MM_count = 0
      for k,LeftBase in enumerate(LeftBases):
        RightBase = RightBases[k]
        if LeftBase == C_base_left:
          if RightBase == C_base_right:
            CC_count += 1
          elif RightBase == M_base_right:
            CM_count += 1
        elif LeftBase == M_base_left:
          if RightBase == C_base_right:
            MC_count += 1
          elif RightBase == M_base_right:
            MM_count += 1
      count = CC_count + CM_count + MC_count + MM_count
      #print(CC_count, CM_count, MC_count, MM_count)
      if count == 0:
        continue

      # We have checked that each of the sites has at least two bases present;
      # however it's possible (though very unlikely) that this is not true in
      # the restricted subset of reads we consider. If it is true, D = 0.
      Cx_count = CC_count + CM_count
      Mx_count = MC_count + MM_count
      xC_count = CC_count + MC_count
      xM_count = CM_count + MM_count
      if (Cx_count == 0 or Mx_count == 0 or xC_count == 0 or xM_count == 0):
        NoLD = True
      else:

        # Let x stand for either C or M. 'freq' means frequency.
        Cx_freq = float(Cx_count)/count
        Mx_freq = float(Mx_count)/count
        xC_freq = float(xC_count)/count
        xM_freq = float(xM_count)/count

        D = float(CC_count)/count - Cx_freq * xC_freq

        if D < 0:
          Dmax = min(Cx_freq*xC_freq, Mx_freq*xM_freq)
          Dprime = D/Dmax
        elif D > 0:
          Dmax = min(Cx_freq*xM_freq, Mx_freq*xC_freq)
          Dprime = D/Dmax
        else:
          NoLD = True

        r = D / (Cx_freq * Mx_freq * xC_freq * xM_freq)**0.5

    if NoLD:
      D = 0
      Dprime = 0
      r = 0

    print(LeftVarPos, RightVarPos, count, D, Dprime, r)


# TODO: calculate p value for at least that much LD in total. Multinomial?


exit(0)
#print(x.seq[x.query_alignment_start])
#print(x.seq)
#print(x.query_name)
#print(x.query_alignment_start,'query_alignment_start')
#print(x.reference_start,'reference_start')
#print()

#base = read.seq[ThisColInThisSeq]
#if base != 'G':
#  print(read.query_name, base)
#print()

iter = BamFile.pileup(TheRef, col-1, col)
for x in iter:
  #print(x.nsegments)
  RefPos = x.reference_pos
  if RefPos != 1001:
    continue
  print(RefPos)
  for pileup in x.pileups:
    read = pileup.alignment
    sequence = read.seq
    print(read.query_name)
  break



'''
    if len(MostCommon2bases_left) == 1:
      [(C_base_left,C_count_left)] = MostCommon2bases_left
      M_count_left = 0
    else:
      [(C_base_left, C_count_left), (M_base_left, M_count_left)] = \
      MostCommon2bases_left
    if len(MostCommon2bases_right) == 1:
      [(C_base_right, C_count_right)] = MostCommon2bases_right
      M_count_right = 0
    else:
      [(C_base_right, C_count_right), (M_base_right, M_count_right)] = \
      MostCommon2bases_right

    if not (M_count_left == 0 or M_count_right == 0):
'''
