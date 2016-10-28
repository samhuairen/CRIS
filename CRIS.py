#!/usr/bin/env python

'''
Search for CRISPR/Cas9 target sites in an annotated genbank sequence.
email dr.mark.schultz@gmail.com
20161028_YYYYMMDD

To do: fix layout of script, move code from main() into 
functions, parallelise searching, pylint, remove redundancy in code
'''

import time
start_time = time.time()
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import re
from collections import defaultdict
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import os

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Find CRISPR/Cas9 target sites \
                                 in a genbank file.  Accepts files with \
                                 multiple genbank records per file.  \
                                 Author: dr.mark.schultz@gmail.com.  \
                                 Acknowledgements: Torsten Seemann, Ian Monk, \
                                 Timothy P. Stinear.')
PARSER.add_argument('-s', '--seq_infile', help='DNA seqs (contigs) to scan \
                    (Genbank format)', required=True)
PARSER.add_argument('-l', '--length_target_seq', help='Set total length \
                    of target sequence NOT including the PAM \
                    sequence. Default=20.', default=20, type=int,
                    required=False)
PARSER.add_argument('-t', '--three_prime_clamp', help='At the 3\' end, how \
                    long do you want the clamp sequence to be? \
                    Default=12.', default=12, type=int, required=False)
PARSER.add_argument('-p', '--PAM_seq', help='Protospacer Adjacent Motif \
                    (PAM).  Depends on Cas9 species. Default=\'NGG\'.',
                    default='NGG', required=False)
PARSER.add_argument('-q', '--feature_qualifier', help='Genbank feature \
                    qualifier in which to find target sites. Could be \
                    \'gene\', \'CDS\', \'mRNA\' etc.  Case-sensitive, exact \
                    spelling required.  Default=\'gene\'.',
                    default='gene', required=False)
PARSER.add_argument('-v', '--verbose', help='Verbose on.  \
                    Default=False', default=False, action='store_true',
                    required=False)
#need to implement this feature
# PARSER.add_argument('-c', '--circularise_off', help='CRIS.py assumes that \
#                     the genbank file contains complete circular contigs from \
#                     chromosomes and/or plasmids, with the start of the \
#                     contig not artificially duplicated at the end of the \
#                     contig (as sometimes happens with assembly of PacBio \
#                     seqs).  If the genbank file contains only a draft \
#                     assembly, it probably does not make sense to treat each \
#                     contig as circular.  For draft assemblies, use this \
#                     switch.  Default=\'False\'.', default=False,
#                     action='store_true', required=False)

GROUP = PARSER.add_mutually_exclusive_group(required=False)
GROUP.add_argument('-n', '--no_overwrite', help='Do not overwrite output file \
                   if it exists.', dest='feature', action='store_false')
GROUP.add_argument('-o', '--overwrite', help='Overwrite output file if it \
                   exists, otherwise write new. Default is to overwrite.',
                   dest='feature', action='store_true')
GROUP.set_defaults(feature=True)
#to do: plot the pairwise SNP distance between the target seqs
ARGS = PARSER.parse_args()

VERSION = '0.1.1'

#The target site rule, define it as a regex
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
SEQS = '[ATGC]{'+str(ARGS.length_target_seq)+'}'+''.join(PAM)
SITE_LENGTH = ARGS.length_target_seq+len(PAM)

def loci_locations(gb_record):
    '''
    From the gb_record, extract locus names, coordinates and strands.
    '''
    locus_locations = defaultdict(list)
    features = []
    for index, feature in enumerate(gb_record.features):
        if hasattr(feature, 'type'):
            if feature.type == ARGS.feature_qualifier:
                #Capture the whole feature in gb_feature.
                features.append(gb_record.features[index])
    for feature in features:
        locus_name = feature.qualifiers.get(ARGS.feature_qualifier,
                                            feature.qualifiers.get('locus_tag'))[0]
        loc_strt = feature.location.start.position
        loc_end = feature.location.end.position
        loc_strnd = feature.location.strand
        locus_locations[locus_name].append([loc_strt, loc_end, loc_strnd])
    return locus_locations

def main():
    '''
    Get sequence, extract genes and coordinates.
    '''
    seqn = SeqIO.parse(open(ARGS.seq_infile, 'r'), 'genbank')
    infile_recs = [gb_record for gb_record in seqn]
    seqs_infile_str = [str(i.seq) for i in infile_recs]
    seqs_infile_revcomp_str = [str(i.seq.reverse_complement()) for i in infile_recs]
    #Create a super-contig of forward and reverse
    #Need to circularise the contigs or not depending on ARGS
    all_gb_records_seq = 'N'.join(seqs_infile_str).upper()
    all_gb_records_seq_rev = 'N'.join(seqs_infile_revcomp_str).upper()
    updated_gb_records = []
    print 'Searching for CRISPR/Cas9 target seqs...\n'
    for gb_record in infile_recs:
        did_not_hit = []
        total_hits = []
        if ARGS.verbose:
            print gb_record
        locus_locs = loci_locations(gb_record)
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
                    if ARGS.verbose:
                        print '\n____\n----\n', locus_name
                    #Store the location of the locus
                    locus_location = feature.location
                    locus_start = locus_location.start
                    locus_end = locus_location.end
                    locus_strand = locus_location.strand
                    if ARGS.verbose:
                        print locus_location #, locus_strand
                    #Store the DNA sequence of the feature
                    gene_seq = gb_feature.extract(gb_record.seq).upper()
                    #Find the possible target sites in the feature, store as a list
                    #SEQS is the regular expression to define the target rule
                    potential_target_seqs = re.findall(SEQS, str(gene_seq))
                    if ARGS.verbose:
                        if len(potential_target_seqs) == 0:
                            print 'No target hits.'
                    if ARGS.verbose:
                        print potential_target_seqs
                    finalist_SeqFeatureObj = []
                    for potential_target_seq in potential_target_seqs:
                        target = potential_target_seq[-ARGS.three_prime_clamp:]
                        whole_gb_forward_target_hits = re.findall(target, str(all_gb_records_seq))
                        whole_gb_rev_target_hits = re.findall(target, str(all_gb_records_seq_rev))
                        fwd_hits = [i for i in whole_gb_forward_target_hits]
                        rev_hits = [i for i in whole_gb_rev_target_hits]
                        if ARGS.verbose:
                            print 'Sequence', target, 'had', str(len(fwd_hits)), 'forward and', str(len(rev_hits)), 'reverse hits'
                        if (len(fwd_hits) + len(rev_hits)) == 1:
                            if locus_strand > 0:
                                #Use finditer to find its position
                                target_pos_iterobj = re.finditer(potential_target_seq, str(gb_record_seq))
                                pos = [(i.start(), i.end()) for i in target_pos_iterobj]
                                strt = pos[0][0]
                                stp = pos[0][1]
                                #overlaps will store a value if start is bounded by two features
                                #use list comprehension instead of nested dict loop
                                overlaps=[True for i in locus_locs.values() if (i[0][0] <= strt <=i[0][1]) or (i[0][0] <= stp <=i[0][1])]
                                if len(overlaps) > 1:
                                    if ARGS.verbose:
                                        print 'Within target-feature, hit overlaps off-target features.'
                                else:
                                    finalist_SeqFeatureObj.append((potential_target_seq, SeqFeature(FeatureLocation(strt, stp), strand=locus_strand, type='misc_binding')))
                            if locus_strand < 0:
                                target_pos_iterobj = re.finditer(potential_target_seq, str(gb_record_seq_rev))
                                pos = [(i.start(), i.end()) for i in target_pos_iterobj]
                                strt = gb_record_seq_len - pos[0][0]
                                stp = gb_record_seq_len - pos[0][1]
                                overlaps=[True for i in locus_locs.values() if (i[0][0] <= (strt-1) <=i[0][1]) or (i[0][0] <= (stp+1) <=i[0][1])]
                                if len(overlaps) > 1:
                                    if ARGS.verbose:
                                        print 'Within target-feature, hit overlaps off-target features.'
                                else:
                                    finalist_SeqFeatureObj.append((potential_target_seq, SeqFeature(FeatureLocation(strt-1, stp+1), strand=locus_strand, type='misc_binding')))
                    if ARGS.verbose:
                        print 'Candidate search complete.'
                        print 'There are', str(len(finalist_SeqFeatureObj)), 'potential sites in', locus_name
                    if len(finalist_SeqFeatureObj) == 0:
                        did_not_hit.append(locus_name)
                    else:
                        #each SeqFeatureObj is a tuple containing the target seq, and the SeqFeature
                        maxGC = max([GC(i[0]) for i in finalist_SeqFeatureObj])
                        highestGC_SeqFeatureObjs = [i for i in finalist_SeqFeatureObj if GC(i[0]) == maxGC]
                        upstream_most_highestGC_SeqFeatureObjs = min([i[1].location.start for i in highestGC_SeqFeatureObjs])
                        #grab just the SeqFeature
                        best_SeqFeatureObj = [i for i in highestGC_SeqFeatureObjs if i[1].location.start == upstream_most_highestGC_SeqFeatureObjs][0]
                        total_hits.append(locus_name)
                        #append the SeqFeature to the gb_record
                        gb_record.features.append(best_SeqFeatureObj[1])
                        print 'Best hit for '+ARGS.feature_qualifier+'\t'+locus_name+'\twas\t'+best_SeqFeatureObj[0]+'\tat pos\t'+str(best_SeqFeatureObj[1].location)
        updated_gb_records.append(gb_record)
        #Tell the user how many features were processed
        if len(n_qualifiers_found) == 0:
            print '\nNo', ARGS.feature_qualifier, 'feature found in genbank record.'
        else:
            print '\nFound hits for', str(len(total_hits)), 'of', str(len(n_qualifiers_found)), ARGS.feature_qualifier, 'features:'
            if ARGS.verbose:
                print ', '.join(total_hits)
            print '\nDid not find hits for:', ', '.join(did_not_hit),'\n'
    return updated_gb_records


if __name__ == '__main__':
    outfile = os.path.splitext(ARGS.seq_infile)[0]+'_CRISPRCas9targetsites.gbk'
    #Print whether overwrite is True or False
    print '\nOverwrite', outfile, '==', ARGS.feature, '\n'
    if ARGS.verbose:
        print '\nSuppress screen printout ==', ARGS.verbose
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
    print '\ntarget seq searched for as regex: '+SEQS+'.'
    print 'target length was', str(SITE_LENGTH), 'bp'
    print '3\' clamp length was', ARGS.three_prime_clamp, 'bp'
    print 'Thank you for using CRIS.py version', VERSION
    print 'email: dr.mark.schultz@gmail.com; github: \'schultzm\'.\n'
    print '--- %s seconds ---' % (time.time() - start_time)


