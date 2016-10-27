#!/usr/bin/env python

'''
Search for 16S rRNA annotations in a genbank sequence. Print '.fasta'
email dr.mark.schultz@gmail.com
20161026_YYYYMMDD
'''

import argparse
from Bio import SeqIO

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Pull out 16S rRNA annotations')
PARSER.add_argument('-s', '--seq_infile', help='Genbank file', required=True)
ARGS = PARSER.parse_args()

def main():
    '''
    Get sequence, extract 16S rRNA annotations.
    '''
    seqn = SeqIO.parse(open(ARGS.seq_infile, 'r'), 'genbank')
    for gb_record in seqn:
        for index, feature in enumerate(gb_record.features):
            if feature.type == 'rRNA':
                gb_feature = gb_record.features[index]
                gene_seq = gb_feature.extract(gb_record.seq).upper()
                if '16S' in feature.qualifiers['product'][0]:
                    print '>'+gb_record.id+'_'+\
                          feature.qualifiers['product'][0].replace(' ','_')+\
                          ' | position '+str(gb_feature.location)
                    print str(gene_seq)

if __name__ == '__main__':
    main()
