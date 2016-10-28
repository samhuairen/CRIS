![alt text][logo]

[logo]: http://i1.wp.com/www.artofthecell.com/wp-content/uploads/2014/06/Art-of-the-Cell-CRISPR-Cas9-in-Complex-with-Guide-RNA-and-target-DNA.jpg "CRISPR/Cas9"


# Find CRISPR/Cas9 binding sites in a genbank file. 
To get help:
```
python CRIS.py -h
usage: CRIS.py [-h] -s SEQ_INFILE [-l LENGTH_TARGET_SEQ]
               [-t THREE_PRIME_CLAMP] [-p PAM_SEQ] [-q FEATURE_QUALIFIER] [-v]
               [-n | -o]

Find CRISPR/Cas9 target sites in a genbank file. Accepts files with multiple
genbank records per file. Author: dr.mark.schultz@gmail.com. Acknowledgements:
Torsten Seemann, Ian Monk, Timothy P. Stinear.

optional arguments:
  -h, --help            show this help message and exit
  -s SEQ_INFILE, --seq_infile SEQ_INFILE
                        DNA seqs (contigs) to scan (Genbank format)
  -l LENGTH_TARGET_SEQ, --length_target_seq LENGTH_TARGET_SEQ
                        Set total length of target sequence NOT including the
                        PAM sequence. Default=20.
  -t THREE_PRIME_CLAMP, --three_prime_clamp THREE_PRIME_CLAMP
                        At the 3' end, how long do you want the clamp sequence
                        to be? Default=12.
  -p PAM_SEQ, --PAM_seq PAM_SEQ
                        Protospacer Adjacent Motif (PAM). Depends on Cas9
                        species. Default='NGG'.
  -q FEATURE_QUALIFIER, --feature_qualifier FEATURE_QUALIFIER
                        Genbank feature qualifier in which to find target
                        sites. Could be 'gene', 'CDS', 'mRNA' etc. Case-
                        sensitive, exact spelling required. Default='gene'.
  -v, --verbose         Verbose on. Default=False
  -n, --no_overwrite    Do not overwrite output file if it exists.
  -o, --overwrite       Overwrite output file if it exists, otherwise write
                        new. Default is to overwrite.
```

Example usage:
```
python CRIS.py -s test_multigbk
```

##Explanation of CRIS.py
Reads in multi-record genbank file. <br><br>
User sets up the PAM sequence, target length and 3'-clamp length and/or accepts the defaults. <br><br>
For each record, searches through the features of type requested on command line (e.g., gene, CDS or mRNA).<br><br>
Within the feature, CRIS.py searches for seqs matching the 3'-clamp sequence (defined by regex rule).<br><br>
After finding matches, CRIS.py assesses whether each of the full length CRISPR/Cas9 target sequences is unique to the feature and that it does not overlap with other features of the requested type.  If a sequence is unique and does not overlap other features, it is stored in the candidate list.  Within the candidate list, the GC content of each CRISPR/Cas9 target sequence is calculated.  Sequences with max GC content are retained.  Of those retained, the CRISPR/Cas9 target sequence closest to the 5' end of the gene is selected as the 'best' match.  The best match is reported.  A summary is printed now and at the end of the run.  In verbose mode, lots of statements are printed as the run progresses.  <br><br>
The output of the run is a copy of the input genbank file with all the best hits marked up in the file.  This annotated genbank can be viewed in Artemis etc.  
