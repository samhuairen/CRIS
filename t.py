from Bio import SeqIO
seqn = SeqIO.parse('/Users/mschultz/Desktop/MDU_DAMG/MDU/CPE_ongoing/IMP4/PROKKA_07042016.gbk', 'genbank')
for record in seqn:
    seq_rev = record.seq.reverse_complement()

