#!/usr/bin/env python

'''
Search for CRISPR binding sites in an annotated genbank sequence.
email dr.mark.schultz@gmail.com
20161025_YYYYMMDD
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
PARSER.add_argument('-s', '--seq_infile', help='DNA seqs (contigs) to scan \
                    (Genbank format)', required=True)

PARSER.add_argument('-l', '--length_CRISPR_seq', help='Set total length \
                    of CRISPR binding site sequence NOT including the PAM \
                    sequence.', default=20, type=int, required=False)
PARSER.add_argument('-c', '--three_prime_clamp', help='At the 3\' end, how \
                    long do you want the clamp sequence to be? \
                    Default=\'12\'.', default=12, type=int, required=False)
PARSER.add_argument('-p', '--PAM_seq', help='Protospacer Adjacent Motif \
                    (PAM).  Depends on Cas9 species. Default=\'NGG\'.',
                    default='NGG', required=False)
PARSER.add_argument('-k', '--gb_feature_key', help='Genbank feature key. \
                    Typically \'gene\', but could be \'CDS\' etc.  Case \
                    sensitive, exact spelling required.', default='gene',
                    required=False)
#Add option for type selection (gene, cds etc)
#Add options for PAM site, and total length of CRISPR seq
#From this extract the total length
#plot the pairwise SNP distance between the CRISPR seqs
ARGS = PARSER.parse_args()

#The CRISPR binding site rule, define it as a regex
PAM = list(ARGS.PAM_seq.upper())
IUPAC_DICT = {'A':'A',
              'C':'C',
              'T':'T',
              'G':'G',
              'R':'[AG]{1}',
              'Y':'[CT]{1}',
              'S':'[GC]{1}',
              'W':'[AT]{1}',
              'K':'[GT]{1}',
              'M':'[AC]{1}',
              'B':'[CGT]{1}',
              'D':'[AGT]{1}',
              'H':'[ACT]{1}',
              'V':'[ACG]{1}',
              'N':'[ACGT]{1}'}
for i in range(0,len(PAM)):
    PAM[i] = IUPAC_DICT[PAM[i]]
SEQS = '[ATGC]{'+str(ARGS.length_CRISPR_seq)+'}'+''.join(PAM)

def gene_locations(gb_record):
    '''
    From the gb_record, extract gene names and coordinates.
    '''
    locus_locations = defaultdict(list)
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
#                    print '\n', locus_name
                #Store the location of the locus
                locus_location = feature.location
#                    print locus_location.start.position
                locus_start = locus_location.start.position
                locus_end = locus_location.end.position
                locus_strand = locus_location.strand
#                    print locus_location, locus_strand
                locus_locations[locus_name].append([locus_start, locus_end, locus_strand])
    return locus_locations


def main():
    '''
    Get sequence, extract genes and coordinates.
    '''
    #Store these objects as a list.  The file handle moves on accession. 
    #Can't reset file handle with stdin
    seqn = SeqIO.parse(open(ARGS.seq_infile, 'r'), 'genbank')
    in_file_recs = [gb_record for gb_record in seqn]
    seqs_in_file_str = [str(i.seq) for i in in_file_recs]
    #seqn = SeqIO.parse(open(ARGS.seq_infile, 'r'), 'genbank')
    seqs_in_file_revcomp_str = [str(i.seq) for i in in_file_recs]
    #Create a super-contig of forward and reverse
    all_gb_records_seq = 'N'.join(seqs_in_file_str).upper()
    all_gb_records_seq_rev = 'N'.join(seqs_in_file_revcomp_str).upper()
    for gb_record in in_file_recs:
        print gb_record
        locus_locs = gene_locations(gb_record)
        print locus_locs
        gb_record_seq = gb_record.seq
        gb_record_seq_len = len(gb_record_seq)
        gb_record_seq_rev = gb_record.seq.reverse_complement()
        for index, feature in enumerate(gb_record.features):
            if hasattr(feature, 'type'):
                if feature.type == ARGS.gb_feature_key:
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
                    #Store the DNA sequence of the feature
                    gene_seq = gb_feature.extract(gb_record.seq).upper()
                    gene_seq_rev = gb_feature.extract(gb_record.seq).upper().reverse_complement()
                    #Find the possible CRISPR sites in the gene, stored as a list
                    #SEQS is a regular expression to define the CRISPR rule
                    potential_CRISPR_seqs = re.findall(SEQS, str(gene_seq))
                    #Check if they hit elsewhere in the all_gb_records_seq or its revcomp
                    if len(potential_CRISPR_seqs) == 0:
                        print 'No CRISPR hits.'
                    finalists = []
                    for potential_CRISPR_seq in potential_CRISPR_seqs:
                        #check firstly for matches at the 12 bases at 3' end
                        whole_gb_forward_CRISPR_hits = re.findall(potential_CRISPR_seq[-ARGS.three_prime_clamp:], str(all_gb_records_seq))
                        whole_gb_rev_CRISPR_hits = re.findall(potential_CRISPR_seq[-ARGS.three_prime_clamp:], str(all_gb_records_seq_rev))
                        fwd_hits = [i for i in whole_gb_forward_CRISPR_hits]
                        rev_hits = [i for i in whole_gb_rev_CRISPR_hits]
                        print 'Sequence', potential_CRISPR_seq, 'had', str(len(fwd_hits)), 'forward and', str(len(rev_hits)), 'reverse hits'

                        if 0 < (len(fwd_hits) + len(rev_hits)) <= 1:
                            if locus_strand > 0:
                                CRISPR_pos_iterobj = re.finditer(potential_CRISPR_seq, str(gb_record_seq))
                                pos = [(i.start(), i.end()) for i in CRISPR_pos_iterobj]
                                strt = pos[0][0]
                                stp = pos[0][1]
                                #overlaps will store a value if site overlaps two genes
                                overlaps = []
                                for key, value in locus_locs.items():
                                    if key == locus_name:
                                        if len(value) > 1:
                                            print 'Note: '+locus_name+' present in '+str(len(value))+' copies.'
                                    if key != locus_name:
                                        if value[0][0] <= strt <= value[0][1]:
                                            print 'This binding seq is in two genes. '
                                            overlaps.append(1)
                                if len(overlaps) > 0:
                                    break
                                else:
                                    recrd = SeqFeature(FeatureLocation(strt, stp), strand=1, type='misc_binding')
                                    gb_record.features.append(recrd)
                                    #got the feature, now break for next gene
                                    break
                            if locus_strand < 0:
                                CRISPR_pos_iterobj = re.finditer(potential_CRISPR_seq, str(gb_record_seq_rev))
                                pos = [(i.start(), i.end()) for i in CRISPR_pos_iterobj]
                                strt = gb_record_seq_len - pos[0][0]
                                stp = gb_record_seq_len - pos[0][1]
                                overlaps = []
                                for key, value in locus_locs.items():
                                    if key == locus_name:
                                        if len(value) > 1:
                                            print 'Note: '+locus_name+' present in '+str(len(value))+' copies.'
                                    if key != locus_name:
                                        if value[0][0] <= strt <= value[0][1]:
                                            print 'This binding seq is in two genes. '
                                            overlaps.append(1)
                                if len(overlaps) > 0:
                                    break
                                else:
                                    recrd = SeqFeature(FeatureLocation(strt, stp), strand=-1, type='misc_binding')
                                    gb_record.features.append(recrd)
                                    #got the feature, now break for next gene
                                    break
        with open(gb_record.id+'_CRISPRsites.gbk', 'w') as outhandle:
            SeqIO.write(gb_record, outhandle, 'genbank')


if __name__ == '__main__':
    main()
    print '\nCRISPR binding seq searched for as regex: '+SEQS+'\n'
