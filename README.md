# PairwiseSequenceAlignment
* Implementation of pairwise Sequence Alignment using dynamic programming approach

## User guide
* User guide: program argument sequence.fasta (or other txt files) containing 2 or more DNA / RNA / amino acid sequences with no degenerated code
* DNA / RNA / amino acid sequences are case insensitive
* Result is exported to a txt file and on screen

## Variants of alignment algorithms
* Global alignment (Needleman-Wunsch algorithm)
* Local alignment (Smith-Waterman algorithm)
* Semi-global alignment (Modified Needleman-Wunsch)
* Overlap alignment (Modified Smith-Waterman)
* Gap penalty at the end of sequence were removed in semi-global and overlap alignment

## Approaches used
* Traceback function used dynamic programming
  * Not required to store the traceback path --> only store the aligned part started from the end of sequences
* Function pointers were used to select different score function for DNA/RNA or protein
  * Users are allowed to alter Match, Mismatch and Gap score for DNA/RNA seuqneces
  * BLOSUM62 scoring matrix was used for amino acid sequences
    * Runtime of amino acid sequences would be a bit longer (Using 2D vector matrix and mapping to obtain score for each comparison)

## Example dataset and results
* Genomic sequence of human and house cat Cyclin B1 genes for DNA alignment
* Hepatitis B virus polymerase protein sequence between wild-type(YMDD motif) and partial sequence from a mutant causing lamivudine-resistance in clinical outcome
* GenBank accession number were included in header of each sequence in example dataset fasta file (Reference of each sequence could obtain from GenBank)

## Limitation of recursive approach
* Recursive approach to find all optimal global alignment was added in RecursiveFunction.cpp
* Common to obtain a large number of alignments with the same optimal score using default scoring setting (+1, -1, -2) --> Not suitable for long sequence in general
* Example demo dataset contained 15 bp --> 12 optimal alignment results
