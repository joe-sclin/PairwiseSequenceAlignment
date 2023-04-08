# PairwiseSequenceAlignment
* Implementation of pairwise Sequence Alignment using dynamic programming approach
* User guide: program argument sequence.fasta (or other txt files) containing 2 or more DNA / RNA sequences with no degenerated code (Only A / C / G / T / U)
* DNA / RNA sequences are case insensitive
* Result is exported to a txt file and on screen
* 4 variants were included
    * Global alignment (Needleman-Wunsch algorithm)
    * Local alignment (Smith-Waterman algorithm)
    * Semi-global alignment (Modified Needleman-Wunsch)
    * Overlap alignment (Modified Smith-Waterman)
    * Gap penalty at the end of sequence were removed in semi-global and overlap alignment
* Traceback function used dynamic programming
  * Not required to store the traceback path --> only store the aligned part started from the end of sequences