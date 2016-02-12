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



CallConsensus.py 

## Overview:
ExplanatoryMessage = '''
This script interprets a base frequency file of the format produced
by AnalysePileup.py and calls the consensus sequence.
Usage: call it from the command line with the first argument being a base
frequency file, the second argument being the minimum number of reads in
agreement with each other before a base is called, and the third argument
being the minimum number of reads in agreement with each other before upper
case is used for the base present.'''



CheckFastaFileEquality.py 

## Overview:
ExplanatoryMessage = '''This script checks that all fasta files supplied as
arguments contain the same sequences, regardless of formatting. It exits with
exit code 111 if any difference is found.'''



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



FindContaminantReadPairs.py 

## Overview: this script reads in one blast file of hits for forward reads and
## one for backward reads (from paired-read data). It finds read pairs for which
## either a) both reads blast to something other than the named sequence, or b)
## one read blasts to the named sequence and its mate blasts to something else
## but that second blast has a better evalue. The names of the reads in these
## pairs are written, without preserving their original order, to two files: one
## for forward reads and one for reverse reads.
## The intended use is for processing the result of blasting paired-read data
## to a database consisting of something that looks like your sample plus other
## sequences that should be recognised as contamination. (For example, if de
## novo assembly has been done with these reads, some of the resulting contigs
## may be identified as contamination -- we want to find the reads that
## correspond to those contigs, in order to remove them.)



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



FindNamedReadsInSortedFastq.py 

## Overview: takes a fastq file sorted by read name, e.g. with
## $ cat MyReads.fq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > out.fq
## and a file of read names sorted in the same way (one read name per line), and
## finds those reads, printing them to the screen suitable for redirection into
## another fastq file.
## Usage:
## $ FindNamedReadsInSortedFastq.py AllReadsSorted.fastq ReadNamesIwant.txt



FindPrimersInAlignment.py 

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



FindSeqsInFasta.py 

## Overview:
ExplanatoryMessage = '''This script retrieves searched-for sequences from a
fasta file. Output is printed to stdout in fasta format, with options to invert the search, extract a window of alignment, and strip gaps (call with --help for details).'''



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



PrintSeqLengths.py 

## Overview:
ExplanatoryMessage = '''This script takes one fasta file as an argument and
prints the length of each sequence found therein, optionally ignoring gaps.'''



RemoveColumnsFromAlignment.py 

## Overview:
ExplanatoryMessage = '''This script removes positions, specified with respect to
one sequence therein, from an alignment.'''



SequenceExtractExtractor.py 

## Overview: run from the command line with one aligned fasta file as an 
## argument, and start & end positions specified internally below, this script
## prints only that part of the alignment between the start & end positions.



SummariseFastqQualities.py 

## Overview:
ExplanatoryMessage = '''Calculates the mean quality, as a function of position
within the read (since quality typically drops at the edges), of all reads found
in all fastq.gz files passed as arguments.'''



TranslateCoords.py 

## Overview: this script translates coordinates with respect to one sequence to
## coordinates with respect to all other sequences in an alignment.
## Usage: call it from the command line thus, for reference-based coordinates:
## ./TranslateCoords.py MyAlignmentFile MyChosenReference coord1 [coord2...]
## or thus, for alignment-based coordinates:
## ./TranslateCoords.py MyAlignmentFile -A coord1 [coord2...]
## The translated coordinates are reported in the order in which they were
## specified.



TrimGapsFromSeqEnds.py 

## Overview:
ExplanatoryMessage = '''This script trims gap characters "-" and "?" from the
ends of all sequences in a fasta file. (Note that if the fasta file is an
alignment, this operation will, in general, unalign the sequences.)'''



TrivialAlign.py 

## Overview: call this script from the command line with a single fasta file as
## an argument, that file containing two sequences. This script moves the
## shorter one along the longer one to find the best match.



UngapFasta.py 

## Overview:
ExplanatoryMessage = '''This script removes the gap character "-" from sequences
in a fasta file. Output is printed to stdout.'''



UpdateIndex.bash 

## Overview: this script produces the index in this directory.



