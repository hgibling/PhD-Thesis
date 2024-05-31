#! /usr/bin/env python

import argparse
import pysam
import sys
import os

from collections import defaultdict

class ReferenceKmer:
    def __init__(self):
        self.positions = list()
        self.read_count = 0

class KmerSegment:
    def __init__(self, name, fasta, coordinates, k):
        self.name = name
        self.sequence = reference.fetch(coordinates[0], coordinates[1], coordinates[2])
        self.coordinates = coordinates
        self.k = k
        self.aligned_reads = 0

        # for "position" analysis
        self.depth_count = [0] * len(self.sequence)
        self.exact_kmer_count = [0] * len(self.sequence)
        self.inexact_kmer_count = [0] * len(self.sequence)

        # for "kmer" analysis
        self.kmer_map = defaultdict(ReferenceKmer)
        for idx in range(len(self.sequence) - k + 1):
            kmer = self.sequence[idx:(idx+args.k)]
            self.kmer_map[kmer].positions.append(idx + self.coordinates[1])

    def consume(self, alignment, debug_fh=None):
        #print(alignment)
        for (query_pos, reference_pos) in alignment.get_aligned_pairs():
            if reference_pos is None or query_pos is None:
                continue

            flank_pos = reference_pos - self.coordinates[1]
            if flank_pos >= 0 and flank_pos < len(self.depth_count):
                
                self.depth_count[flank_pos] += 1
                read_kmer = alignment.query_sequence[query_pos:(query_pos + self.k)]
                ref_kmer = self.sequence[flank_pos:flank_pos+self.k]

                if debug_fh is not None and len(read_kmer) == self.k:
                    correct = ref_kmer == read_kmer
                    debug_fh.write("%s\t%d\t%s\t%s\t%d\n" % (alignment.query_name, reference_pos, ref_kmer, read_kmer, correct))

                if len(read_kmer) == self.k:
                    self.exact_kmer_count[flank_pos] += read_kmer == ref_kmer
                    self.inexact_kmer_count[flank_pos] += read_kmer != ref_kmer

    def consume_kmers(self, alignment):
        for idx in range(len(alignment.query_sequence) - self.k + 1):
            kmer = alignment.query_sequence[idx:(idx+args.k)]
            self.kmer_map[kmer].read_count += 1

def run_position_analysis():

    debug_fh = None
    if args.debug_filename != "":
        debug_fh = open("read_kmer_position.tsv", 'w')
    save = pysam.set_verbosity(0)
    bam_file = pysam.AlignmentFile(args.bam, "rb")
    pysam.set_verbosity(save)
    for alignment in bam_file.fetch(until_eof=True):
        for r in regions:
            r.consume(alignment, debug_fh)

    print("chromosome\tposition\ttype\tbase_depth\texact_kmer_count\tinexact_kmer_count\tsum_kmer_count")
    for r in regions:
        for (idx, depth) in enumerate(r.depth_count):
            ekc = r.exact_kmer_count[idx]
            ikc = r.inexact_kmer_count[idx]
            print("%s\t%d\t%s\t%d\t%d\t%d\t%d" % (r.coordinates[0], idx + r.coordinates[1], r.name, depth, ekc, ikc, ekc + ikc))

def run_kmer_analysis():
    save = pysam.set_verbosity(0)
    bam_file = pysam.AlignmentFile(args.bam, "rb")
    pysam.set_verbosity(save)
    for alignment in bam_file.fetch(until_eof=True):
        for r in regions:
            r.consume_kmers(alignment)

    print("kmer\ttype\tpositions\tread_count")
    for r in regions:
        for (kmer, o) in r.kmer_map.items():
            if len(o.positions) > 0:
                print("%s\t%s\t%s\t%d" % (kmer, r.name, ",".join([str(x) for x in o.positions]), o.read_count))

parser = argparse.ArgumentParser()
parser.add_argument('--reference', type=str, required=True)
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--k', type=int, default=71)
parser.add_argument('--analysis', type=str, default="positions")
parser.add_argument('--debug-filename', type=str, default="")

args = parser.parse_args()

reference = pysam.FastaFile(args.reference)

left_region = KmerSegment("left", reference, ("5", 23516784, 23526780), args.k)
znf_region = KmerSegment("znf", reference, ("5", 23526788, 23527872), args.k)
right_region = KmerSegment("right", reference, ("5", 23527880, 23537869), args.k)
regions = [ left_region, znf_region, right_region ]

if args.analysis == "positions":
    run_position_analysis()
elif args.analysis == "kmer":
    run_kmer_analysis()
