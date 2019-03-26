#!/usr/bin/python
# Samuel Moijueh
# July 2016

import sys
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_dna
import time

###########################################################
'''
== SYNOPOSIS ==
This Program reads takes a list of HGMD transcriptIDs, and for each gene
it gets the sequence for every exon less than 400bp in length from Alamut.
Skip transcriptIDs that we have already created primers for.

OUTPUT: Alamut Chromosomal file with sequence:
a) NM_transcript_primers.bed
(chromosomal file with exon sequences that are >=400bps and from coding transcripts.
	This is the file of interest)

b) NR_transcript_primers.bed
(chromosomal file with exon sequences that are >=400bps AND
	either 1. contain unknown bases, or 2. are from noncoding transcripts)


USAGE: python 1_Extract_ExonSeq.py

'''
#########################################################


def main():
    start = time.time()
    genomeSeq = get_entire_sequence()
    assign_sequence(genomeSeq)
    end = time.time()
    print str(end - start) + " seconds"
    # sys.exit();


# read in the hg19 reference genome file (hg19.fa) into a dictionary where the
# KEY is the chromosome and the VALUE is the genomic sequence
# this function takes ~33 seconds
def get_entire_sequence():
    handle = open("input/hg19.fa", "rU")
    genomeSeq = dict()
    for record in SeqIO.parse(handle, "fasta"):
        genomeSeq[record.id] = record.seq
    return genomeSeq


# extracts exon sequences from the human geonme that we would like to create primers for
def assign_sequence(genomeSeq):

    # list of transcriptIDs that we have already created primers for
    HGMD_transcripts_to_skip = [line.strip().split()[3] for line in open("input/AlamutHGMD_lessthan400_noN_noNR_notinTULIP.bed", 'r')]

    noncoding_primers_fh = open('NR_transcript_primers.txt', 'w')
    coding_primers_fh = open('NM_transcript_primers.txt', 'w')

    # extract the sequence of each gene/transcriptID listed in target HGMD transcripts
    with open('input/alamut_feb_1_2015_exons.txt', 'r') as alamut_fh:
        next(alamut_fh)  # skip header line
        for line in alamut_fh:
            column = line.split()

            if (column[7] not in HGMD_transcripts_to_skip):
                #column[7] is transcriptID
                chromosome = str(column[3])
                chrStart = int(column[13]) - 1
                chrEnd = int(column[14])
                geneName = column[1]
                exonNum = column[12]
                strand = str(column[6])
                transcriptID = column[7]
                # skip over mitochrondrial chromosomes
                if (chromosome == 'M'):
                    continue

                # skip exons greater than 400 bps; for now we're only interested in
                # creating primers for exons less than or equal to 400
                if (chrEnd - chrStart) >= 400:
                    # print "transcriptID = " + transcriptID + " ; exon num = " + str(exonNum) + "\n";
                    # print "exon length = " + str(chrEnd - chrStart) + "\n";
                    continue

                # write noncoding transcript exons to file
                if 'NR' in transcriptID:
                    sequence = str(genomeSeq[chromosome][chrStart - 320: chrEnd - 320]);
                    newChrStart = str(chrStart - 320)
                    newChrEnd = str(chrEnd + 320)
                    noncoding_primers_fh.write(chromosome + '\t' + newChrStart + '\t' + newChrEnd + '\t' + transcriptID +
                                               '\t' + strand + '\t' + sequence + '\n')
                    continue

                # ASSERT: transcripts beyond this point are coding AND less than or equal to 400 base pairs
                sequence = genomeSeq[chromosome][chrStart - 320: chrEnd + 320];
                sequence_string = str(sequence)
                # if the sequence has an unknown nucleotide base, we consider the exon for
                # that transcript as 'noncoding', and it's sequence is written to file
                if sequence.find('N') != -1:
                    newChrStart = str(chrStart - 320)
                    newChrEnd = str(chrEnd + 320)
                    noncoding_primers_fh.write(chromosome + '\t' + newChrStart + '\t' + newChrEnd + '\t' + transcriptID +
                                               '\t' + exonNum + '\t' + geneName + '\t' + strand + '\t' + sequence_string + '\n')
                else:  # if the sequence does not have an N in it then proceed
                    chrStart = str(chrStart)
                    chrEnd = str(chrEnd)
                    if strand == '1':  # for genes on the positive strand
                        coding_primers_fh.write(chromosome + '\t' + chrStart + '\t' + chrEnd + '\t' + transcriptID +
                                                '\t' + exonNum + '\t' + geneName + '\t' + '+' + '\t' + sequence_string + '\n')
                    else:  # for genes on the negative strand
                        my_seq = Seq(sequence_string, generic_dna)
                        reverse_comp_seq = str(my_seq.reverse_complement())
                        chrStart = str(chrStart)
                        chrEnd = str(chrEnd)
                        coding_primers_fh.write(chromosome + '\t' + chrStart + '\t' + chrEnd + '\t' + transcriptID +
                                                '\t' + exonNum + '\t' + geneName + '\t' + '-' + '\t' + reverse_comp_seq + '\n')
    # close the output file
    noncoding_primers_fh.close()
    coding_primers_fh.close()

# end of assign_sequence();
#########################################################



# Program entry point
if __name__ == "__main__":
    main()
