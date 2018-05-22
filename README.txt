AddGapsToAlignment.py 

## Overview:
ExplanatoryMessage = '''This script adds gaps into an alignment. The coordinates
are moved through from smallest to largest; at each coordinate, that column of
the alignment is moved one position to the right and a column of gaps is
inserted. For example, specifying the coordinates 1, 3 and 4 would insert one 
column of gaps onto the left of the alignment, and two columns of gaps in
between the first and second columns.'''



AuxiliaryFunctions.py 

## Overview: here we define some functions: for reading in all sequences from
## a file into a dictionary, for reading in patient details from a file
## into a dictionary, and for checking if two bases match while allowing for one
## or both to be ambiguous.



CalculateEntropy.py 

## Overview: this script, called from the command line with a fasta file as an
## argument, calculates the entropy of the (nominally aligned) set of sequences
## contained therein. IUPAC ambiguity codes are interpreted as counting for half
## or one third of each of their associated letters. The characters '-', 'n' and
## 'N' are deemed to be equivalent and are treated as a fifth distinct letter.



CombineContigsFrom2HiseqLanes.py 

## Overview:
ExplanatoryMessage = '''This script merges all contigs assembled from two Hiseq
lanes, each labelled with a Sanger Institute lane ID, into a single new fasta
file with uniquely labelled contigs (replacing _1_ and _2_ with _3_, and
appending _r for contigs coming from the _2_ lane).'''



CompareDataFilesPositionwise.py 

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



FindBestKrakenRef.py 

## Overview:
ExplanatoryMessage = '''This script reads a Kraken report file and returns the
name of the taxon that has the most reads rooted under it, only looking inside a
specified (larger) taxon. If no reads are rooted to any taxon inside the
specified larger taxon, N/A is printed.'''



FindInsertSizeSpikes.py 

## Overview:
ExplanatoryMessage = '''This script finds those read pairs in a bam file for
 which the insert size is a user-specified value. Amongst these read pairs, the
 number of times each unique read occurs is counted. The most common reads (call
 with the -h flag for more detail) are printed to stdout in fasta format. The
intended purpose is to find cases where one read pair has been exactly
duplicated a large number of times, leading to a spike in an otherwise smooth
insert-size distribution.
'''



FindPairwiseGapFracsFromAlignment.py 

## Overview:
ExplanatoryMessage = '''This script, taking an sequence alignment in fasta
format as input, considers all possible pairs of sequences therein and
calculates the gap fraction of each sequence if those two sequences were aligned
on their own, i.e. not counting any position at which both sequences have a gap.
(Instead of actually removing such positions, the calculation is made using a
look-up table of hashed gap coordinates, easing the pain of the unavoidable
O(N^2) runtime for N sequences by improving speed 50-fold.)
'''



FindUniqueInsertions.py 

## Overview:
ExplanatoryMessage = '''Given an alignment of sequences, this script prints the
(one-based) coordinates of positions (if any exist) at which a named sequence
has unique insertions, i.e. positions where only that sequence has bases, all
other sequences have gaps. '''



LinkageDiseq.py 

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



MergeFastaFiles.py 

## Overview:
ExplanatoryMessage = '''This script takes one or more fasta files as arguments,
and reads the sequences contained therein. Sequences encountered a second time -
i.e. with the same name - are checked to be identical. All sequences are then
printed to stdout, in the same order in which they appeared in the files, with
duplicates skipped.'''



PropagateNoCoverageChar.py 

## Overview:
ExplanatoryMessage = '''This script replaces any "-" character (taken to mean a deletion at this position) that neighbours a "?" character (taken to mean unknown / missing data) by a "?" character, iteratively until no more replacements need to be made. The idea is that a deletion should not be called if it is not known what's next to it. Output is printed to stdout, suitable for redirection to a fasta file.'''



RankSequencesBySimilarityToOne.py 

## Overview:
ExplanatoryMessage = '''Calculates the similarity of every other sequence to
one named sequence in a fasta-format alignment. Comparing each other sequence in
turn we find the number of matching bases in the region of the alignment
covered by both sequences (ignoring missing coverage - "?" positions), the
length of that region, and the ratio of the former to the latter i.e. the
fractional identity in the overlap. Results are printed to stdout, in csv
format, sorted by this ratio. (For non-overlap this ratio is 0/0: we report 0.)
'''



RemoveColumnsFromAlignment.py 

## Overview:
ExplanatoryMessage = '''This script removes positions, specified with respect to
one sequence therein, from an alignment.'''



SeqFormatConverter.py 

## Overview:
ExplanatoryMessage = '''This script converts sequence data between formats; it
is a simple command-line wrapper for Bio.SeqIO's convert function. Output is
printed to stdout, suitable for redirection to a file.'''



SequenceExtractExtractor.py 

## Overview: run from the command line with one aligned fasta file as an 
## argument, and start & end positions specified internally below, this script
## prints only that part of the alignment between the start & end positions.



SummariseFastqQualities.py 

## Overview:
ExplanatoryMessage = '''Calculates the mean quality, as a function of position
within the read (since quality typically drops at the edges), of all reads found
in all fastq.gz files passed as arguments.'''



TrimGapsFromSeqEnds.py 

## Overview:
ExplanatoryMessage = '''This script trims gap characters "-" and "?" from the
ends of all sequences in a fasta file. (Note that if the fasta file is an
alignment, this operation will, in general, unalign the sequences.)'''



TrivialAlign.py 

## Overview: call this script from the command line with a single fasta file as
## an argument, that file containing two sequences. This script moves the
## shorter one along the longer one to find the best match.



UpdateIndex.bash 

## Overview: this script produces the index in this directory.



