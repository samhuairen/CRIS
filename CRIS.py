#!/usr/bin/env python

'''
Search for CRISPR binding sites in an annotated genbank sequence.
email dr.mark.schultz@gmail.com
20161026_YYYYMMDD
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
import os

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Find CRISPR binding sites \
                                 in a genbank file.  Accepts files with \
                                 multiple genbank records per file.  \
                                 Author: dr.mark.schultz@gmail.com.  \
                                 Acknowledgements: Torsten Seemann, Ian Monk, \
                                 Timothy P. Stinear.')
PARSER.add_argument('-s', '--seq_infile', help='DNA seqs (contigs) to scan \
                    (Genbank format)', required=True)
PARSER.add_argument('-l', '--length_CRISPR_seq', help='Set total length \
                    of CRISPR binding site sequence NOT including the PAM \
                    sequence. Default=20.', default=20, type=int,
                    required=False)
PARSER.add_argument('-t', '--three_prime_clamp', help='At the 3\' end, how \
                    long do you want the clamp sequence to be? \
                    Default=12.', default=12, type=int, required=False)
PARSER.add_argument('-p', '--PAM_seq', help='Protospacer Adjacent Motif \
                    (PAM).  Depends on Cas9 species. Default=\'NGG\'.',
                    default='NGG', required=False)
PARSER.add_argument('-q', '--feature_qualifier', help='Genbank feature \
                    qualifier in which to find binding sites. Could be \
                    \'gene\', \'CDS\', \'mRNA\' etc.  Case-sensitive, exact \
                    spelling required.  Default=\'gene\'.',
                    default='gene', required=False)
PARSER.add_argument('-x', '--suppress_screen_output', help='Verbose off.  \
                    Default=False', default=False, action='store_true',
                    required=False)
PARSER.add_argument('-c', '--circularise_off', help='CRIS.py assumes that \
                    the genbank file contains complete circular contigs from \
                    chromosomes and/or plasmids, with the start of the \
                    contig not artificially duplicated at the end of the \
                    contig (as sometimes happens with assembly of PacBio \
                    seqs).  If the genbank file contains only a draft \
                    assembly, it probably does not make sense to treat each \
                    contig as circular.  For draft assemblies, use this \
                    switch.  Default=\'False\'.', default=False,
                    action='store_true', required=False)

GROUP = PARSER.add_mutually_exclusive_group(required=False)
GROUP.add_argument('-n', '--no_overwrite', help='Do not overwrite output file \
                   if it exists.', dest='feature', action='store_false')
GROUP.add_argument('-o', '--overwrite', help='Overwrite output file if it \
                   exists, otherwise write new. Default is to overwrite.',
                   dest='feature', action='store_true')
GROUP.set_defaults(feature=True)
#plot the pairwise SNP distance between the CRISPR seqs
ARGS = PARSER.parse_args()

VERSION = '0.1.1'

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
SITE_LENGTH = ARGS.length_CRISPR_seq+len(PAM)

def gene_locations(gb_record):
    '''
    From the gb_record, extract locus names, coordinates and strands.
    '''
    locus_locations = defaultdict(list)
    features = []
    for index, feature in enumerate(gb_record.features):
        if hasattr(feature, 'type'):
            if feature.type == 'gene':
                #Capture the whole feature in gb_feature.
                features.append(gb_record.features[index])
    for i in features:
        locus_name = feature.qualifiers.get(ARGS.feature_qualifier)[0]
        if locus_name == None:
            locus_name = feature.qualifiers.get('locus_tag')[0]
        locus_location = feature.location
        loc_strt = locus_location.start.position
        loc_end = locus_location.end.position
        loc_strnd = locus_location.strand
        locus_locations[locus_name].append([loc_strt, loc_end, loc_strnd])
    return locus_locations


def main():
    '''
    Get sequence, extract genes and coordinates.
    '''
    seqn = SeqIO.parse(open(ARGS.seq_infile, 'r'), 'genbank')
    infile_recs = [gb_record for gb_record in seqn]
    seqs_infile_str = [str(i.seq) for i in infile_recs]
    seqs_infile_revcomp_str = [str(i.seq) for i in infile_recs]
    #Create a super-contig of forward and reverse
    #Need to circularise the contigs or not depending on ARGS
    all_gb_records_seq = 'N'.join(seqs_infile_str).upper()
    all_gb_records_seq_rev = 'N'.join(seqs_infile_revcomp_str).upper()
    updated_gb_records = []
    for gb_record in infile_recs:
        if not ARGS.suppress_screen_output:
            print gb_record
        locus_locs = gene_locations(gb_record)
        gb_record_seq = gb_record.seq
        gb_record_seq_len = len(gb_record_seq)
        gb_record_seq_rev = gb_record.seq.reverse_complement()
        n_qualifiers_found = []
        for index, feature in enumerate(gb_record.features):
            if hasattr(feature, 'type'):
                if feature.type == ARGS.feature_qualifier:
                    n_qualifiers_found.append(1)
                    #Capture the whole feature in gb_feature.
                    gb_feature = gb_record.features[index]
                    #Store a name for the feature
                    if ARGS.feature_qualifier in feature.qualifiers:
                        locus_name = feature.qualifiers[ARGS.feature_qualifier][0]
                    else:
                        locus_name = feature.qualifiers['locus_tag'][0]
                    if not ARGS.suppress_screen_output:
                        print '\n', locus_name
                    #Store the location of the locus
                    locus_location = feature.location
                    locus_start = locus_location.start
                    locus_end = locus_location.end
                    locus_strand = locus_location.strand
                    if not ARGS.suppress_screen_output:
                        print locus_location #, locus_strand
                    #Store the DNA sequence of the feature
                    gene_seq = gb_feature.extract(gb_record.seq).upper()
                    gene_seq_rev = gb_feature.extract(gb_record.seq).upper().reverse_complement()
                    #Find the possible CRISPR sites in the gene, stored as a list
                    #SEQS is a regular expression to define the CRISPR rule
                    potential_CRISPR_seqs = re.findall(SEQS, str(gene_seq))
                    #Check if they hit elsewhere in the all_gb_records_seq or its revcomp
                    if not ARGS.suppress_screen_output:
                        if len(potential_CRISPR_seqs) == 0:
                            print 'No CRISPR hits.'
                    finalists = []
                    for potential_CRISPR_seq in potential_CRISPR_seqs:
                        #check firstly for matches at the 12 bases at 3' end
                        whole_gb_forward_CRISPR_hits = re.findall(potential_CRISPR_seq[-ARGS.three_prime_clamp:], str(all_gb_records_seq))
                        whole_gb_rev_CRISPR_hits = re.findall(potential_CRISPR_seq[-ARGS.three_prime_clamp:], str(all_gb_records_seq_rev))
                        fwd_hits = [i for i in whole_gb_forward_CRISPR_hits]
                        rev_hits = [i for i in whole_gb_rev_CRISPR_hits]
                        if not ARGS.suppress_screen_output:
                            print 'Sequence', potential_CRISPR_seq, 'had', str(len(fwd_hits)), 'forward and', str(len(rev_hits)), 'reverse hits'
                        if 0 < (len(fwd_hits) + len(rev_hits)) <= 1:
                            if locus_strand > 0:
                                CRISPR_pos_iterobj = re.finditer(potential_CRISPR_seq, str(gb_record_seq))
                                pos = [(i.start(), i.end()) for i in CRISPR_pos_iterobj]
                                strt = pos[0][0]
                                stp = pos[0][1]
                                #overlaps will store a value if site overlaps two genes
                                overlaps = []
                                #use list comprehension instead of nested dict loop
                                #k=[True for i in locus_locs.values() if i[0] <= strt <=i[1]]
                                for key, value in locus_locs.items():
                                    if key == locus_name:
                                        if len(value) > 1:
                                            if not ARGS.suppress_screen_output:
                                                print 'Note: '+locus_name+' present in '+str(len(value))+' copies.'
                                    if key != locus_name:
                                        if value[0][0] <= strt <= value[0][1]:
                                            if not ARGS.suppress_screen_output:
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
                                            if not ARGS.suppress_screen_output:
                                                print 'Note: '+locus_name+' present in '+str(len(value))+' copies.'
                                    if key != locus_name:
                                        if value[0][0] <= strt <= value[0][1]:
                                            if not ARGS.suppress_screen_output:
                                                print 'This binding seq is in two genes. '
                                            overlaps.append(1)
                                if len(overlaps) > 0:
                                    break
                                else:
                                    recrd = SeqFeature(FeatureLocation(strt, stp), strand=-1, type='misc_binding')
                                    gb_record.features.append(recrd)
                                    #got the feature, now break for next gene
                                    break
        updated_gb_records.append(gb_record)
        #Tell the user how many features were processed
        if len(n_qualifiers_found) == 0:
            print '\nNo', ARGS.feature_qualifier, 'feature found in genbank record.'
        else:
            print '\n', str(len(n_qualifiers_found)), ARGS.feature_qualifier, 'features processed.\n'
    return updated_gb_records


if __name__ == '__main__':
    outfile = os.path.splitext(ARGS.seq_infile)[0]+'_CRISPRsites.gbk'
    #Is overwrite True or False? Print the ARGS.feature answer here
    print '\nOverwrite', outfile, '==', ARGS.feature
    if ARGS.suppress_screen_output:
        print '\nVerbose off ==', ARGS.suppress_screen_output
        print 'Search is now running ...\n'
    if ARGS.feature == False:
        if os.path.exists(outfile):
            print outfile, 'already exists. Use \'overwrite\' or move', outfile
            sys.exit()
    gb_data_to_write = main()
    #As overwrite defaults to true, just overwrite
    with open(outfile, 'w') as output_handle:
        for i in gb_data_to_write:
            SeqIO.write(i, output_handle, 'genbank')
    print '\nWritten data to', outfile
    print '\nCRISPR binding seq searched for as regex: '+SEQS+'.'
    print 'CRISPR target length was', str(SITE_LENGTH), 'nt'
    print 'Thank you for using CRIS.py version', VERSION
    print 'email: dr.mark.schultz@gmail.com; github: \'schultzm\'.'
