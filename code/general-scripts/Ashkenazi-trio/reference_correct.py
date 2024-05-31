#! /usr/bin/env python


import parasail
import argparse
import edlib
import pysam
import sys
import os

def get_sequences_from_fasta(fn):
    out = list()
    fh = pysam.FastxFile(fn)
    for s in fh:
        out.append(s)
    return out

def reverse_complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(s))
    return reverse_complement

parser = argparse.ArgumentParser()
parser.add_argument('--left-flank', type=str, required=True)
parser.add_argument('--right-flank', type=str, required=True)
parser.add_argument('--alleles', type=str, required=True)
parser.add_argument('--reads', type=str, required=True)
parser.add_argument('--flank-length', type=int, default=200)
parser.add_argument('--debug', type=int, default=0)
parser.add_argument('--max-edit-distance', type=int, default=20)
args = parser.parse_args()

left_flank = get_sequences_from_fasta(args.left_flank)[0]
right_flank = get_sequences_from_fasta(args.right_flank)[0]
alleles = get_sequences_from_fasta(args.alleles)

# make the alleles with flanking sequence
for a in alleles:
    a.sequence = left_flank.sequence[-args.flank_length:] + a.sequence + right_flank.sequence[0:args.flank_length]

fh = pysam.FastxFile(args.reads)
for read in fh:
    # align this read against all alleles
    best_distance = args.max_edit_distance
    best_allele_idx = -1
    best_strand = 0
    best_correction = read.sequence
    for strand_idx, strand_sequence in enumerate([read.sequence, reverse_complement(read.sequence)]):
        for allele_idx, allele in enumerate(alleles):
            result = edlib.align(strand_sequence, allele.sequence, mode = "HW", task = "path")
            if result['editDistance'] < best_distance:

                best_distance = result['editDistance']
                best_strand = strand_idx
                best_idx = allele_idx
                start = result['locations'][0][0]
                end = result['locations'][0][1]
                best_correction = allele.sequence[start:end]

                if args.debug:
                    out_alignment = edlib.getNiceAlignment(result, strand_sequence, allele.sequence)
                    print(allele.name, strand_idx, allele_idx)
                    print(result)
                    print(out_alignment['query_aligned'])
                    print(out_alignment['matched_aligned'])
                    print(out_alignment['target_aligned'])
                    print(best_correction)

    out_sequence = read.sequence
    if best_correction != "":
        if best_strand == 0:
            out_sequence = best_correction
        else:
            out_sequence = reverse_complement(best_correction)
    print(">%s\n%s" % (read.name, out_sequence))
