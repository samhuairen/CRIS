#!/usr/bin/env python

'''
Search for CRISPR binding sites in an annotated genbank sequence.
'''

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import re
from collections import defaultdict
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Find CRISPR binding sites.')
PARSER.add_argument('-s', '--seq_infile', help='DNA seq to scan.',
                    required=True)
# PARSER.add_argument('-p', '--PAM_seq', help='5\'to 3\' sequence at the\
#                      3\' end of the genomic target', default='NGG', required=False)
# PARSER.add_argument('-t', '--length_target', help='Length of genomic target.',
#                     required=False)

ARGS = PARSER.parse_args()

#The CRISPR binding site rule
SEQS = '[ATGC]{22}GG'



def readseq(seqinfile, seqinformat):
    '''
    Read in the sequence from file.  Return it as the variable 'sequence'.
    '''
    #test if there is more than one record, will raise ValueError if True
    SeqIO.read(seqinfile, seqinformat)
    #Continue if not more than one record using parse
    sequence = SeqIO.parse(seqinfile, seqinformat)
    return sequence

def gene_locations(sequence):
    '''
    From the sequence, extract gene names and coordinates.
    '''
    locus_locations = defaultdict(list)
    for gb_record in sequence:
        for index, feature in enumerate(gb_record.features):
            if hasattr(feature, 'type'):
                if feature.type == 'gene':
                    #Capture the whole feature in gb_feature.
                    gb_feature = gb_record.features[index]
                    #Store a name for the feature if it is a gene
                    if 'gene' in feature.qualifiers:
                        locus_name = feature.qualifiers['gene'][0]
                    else:
                        locus_name = feature.qualifiers['locus_tag'][0]
                    print '\n', locus_name
                    #Store the location of the locus
                    locus_location = feature.location
                    print locus_location.start.position
                    locus_start = locus_location.start.position
                    locus_end = locus_location.end.position
                    locus_strand = locus_location.strand
                    print locus_location, locus_strand
                    locus_locations[locus_name].append([locus_start, locus_end, locus_strand])
    return locus_locations


def main():
    '''
    Get sequence, extract genes and coordinates.
    '''
    seqn = readseq(ARGS.seq_infile, 'genbank')

    locus_locations = gene_locations(seqn)
    print locus_locations
    for gb_record in seqn:
        for index, feature in enumerate(gb_record.features):
            if hasattr(feature, 'type'):
                if feature.type == 'gene':
                    #Capture the whole feature in gb_feature.
                    gb_feature = gb_record.features[index]
                    #Store a name for the feature if it is a gene
                    if 'gene' in feature.qualifiers:
                        locus_name = feature.qualifiers['gene'][0]
                    else:
                        locus_name = feature.qualifiers['locus_tag'][0]
                    print '\n', locus_name
                    #Store the location of the locus
                    locus_location = feature.location
                    locus_start = locus_location.start
                    locus_end = locus_location.end
                    locus_strand = locus_location.strand
                    print locus_location #, locus_strand
                    locus_location[locus_name].append([locus_start, locus_end])

                    #Store the DNA sequence of the feature
                    gene_seq = gb_feature.extract(gb_record.seq).upper()
                    gene_seq_rev = gb_feature.extract(gb_record.seq).upper().reverse_complement()
                    #Find the possible CRISPR sites in the gene, stored as a list
                    #SEQS is a regular expression to define the CRISPR rule
                    potential_CRISPR_seqs = re.findall(SEQS, str(gene_seq))
                    #Check if they hit elsewhere in the whole_gb_record_seq or its revcomp
                    for potential_CRISPR_seq in potential_CRISPR_seqs:
                        whole_gb_forward_CRISPR_hits = re.findall(potential_CRISPR_seq, str(whole_gb_record_seq))
                        whole_gb_rev_CRISPR_hits = re.findall(potential_CRISPR_seq, str(whole_gb_record_seq_rev))
                        fwd_hits = [i for i in whole_gb_forward_CRISPR_hits]
                        rev_hits = [i for i in whole_gb_rev_CRISPR_hits]
                        print fwd_hits, rev_hits
                        if len(fwd_hits) + len(rev_hits) == 1:
                            if locus_strand > 0:
                                CRISPR_pos_iterobj = re.finditer(potential_CRISPR_seq, str(whole_gb_record_seq))
                                pos = [(i.start(), i.end()) for i in CRISPR_pos_iterobj]
                                recrd = SeqFeature(FeatureLocation(pos[0][0], pos[0][1]), type='misc_binding')
                                gb_record.features.append(recrd)
                                #got the feature, now break for next gene
                                break
                            if locus_strand < 0:
                                CRISPR_pos_iterobj = re.finditer(potential_CRISPR_seq, str(whole_gb_record_seq_rev))
                                pos = [(i.start(), i.end()) for i in CRISPR_pos_iterobj]
                                strt = whole_gb_record_seq_len - pos[0][0]
                                stp = whole_gb_record_seq_len - pos[0][1]
                                recrd = SeqFeature(FeatureLocation(strt, stp), strand=-1, type='misc_binding')
                                gb_record.features.append(recrd)
                                #got the feature, now break for next gene
                                break
    with open(ARGS.seq_infile+'edit.gbk', 'w') as outhandle:
        SeqIO.write(gb_record, outhandle, 'genbank')


if __name__ == '__main__':
    main()
    
    
