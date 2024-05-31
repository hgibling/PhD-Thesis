#! /usr/bin/env python

import argparse
import math
import pysam
import sys
import os
import pandas

from collections import defaultdict

class ReferenceSegment:
    def __init__(self, name, fasta, coordinates):
        self.name = name
        self.sequence = reference.fetch(coordinates[0], coordinates[1], coordinates[2]).upper()
        self.coordinates = coordinates
        self.sum_align_length = 0
        self.sum_differences = 0

    def consume(self, alignment, debug_fh=None):
        # ignore reads with mapq of 0 (no read sequence in bam)
        if alignment.query_length == 0:
            return
        aligned_pairs = alignment.get_aligned_pairs()
        # we need to trim the pairs to this region
        # we can't do this the trivial way using the coordinates
        # get_aligned_pairs returns since reference_pos may be None
        # for inserted bases
        first_idx = None
        last_idx = None
        for idx, (query_pos, reference_pos) in enumerate(aligned_pairs):
            if reference_pos is not None and reference_pos >= self.coordinates[1] and reference_pos <= self.coordinates[2]:
                if first_idx == None or idx < first_idx:
                    first_idx = idx
                if last_idx == None or idx > last_idx:
                    last_idx = idx
        
        # ignore reads not in region being tested
        if first_idx == None:
            return

        # TODO debug output
        #     if debug_fh is not None:
        #         debug_fh.write("%s\t%s\t%s\t%d\n" % (alignment.query_name, reference_pos, query_pos, int(is_error)))

        idx = 0
        use_aligned_pairs = aligned_pairs[first_idx:last_idx] # all aligned bases in region of interest
        all_var_positions = list(germline_variants.Position)

        while idx < len(use_aligned_pairs):
            query_pos = use_aligned_pairs[idx][0]
            reference_pos = use_aligned_pairs[idx][1]
            is_error = False

            # check for deletion 
            # (should only reach here if previous refpos was not in germline list)
            if query_pos is None:
                is_error = True
                self.sum_align_length += 1
                idx += 1

            # check for insertion 
            # (should only reach here if previous refpos was not in germline list)
            elif reference_pos is None:

                # insertion at start of alignment. find leftmost germline variant and check if insertion within range
                if idx == 0:
                    # get first instance of ref_pos
                    list_of_ref_pos_none = list(map(lambda i: i[1], use_aligned_pairs))
                    # convert None to inf to use min func
                    list_of_ref_pos = [math.inf if j == None else j for j in list_of_ref_pos_none]
                    rightmost_ref_pos = min(list_of_ref_pos)
                    rightmost_ref_idx = list_of_ref_pos.index(rightmost_ref_pos)

                    # find closest upstream germline variant position
                    var_pos_differences = [(rightmost_ref_pos + 1) - i for i in all_var_positions]
                    leftmost_var_idx = var_pos_differences.index(min([j for j in var_pos_differences if j > 0]))
                    leftmost_var_pos = all_var_positions[leftmost_var_idx]
                    leftmost_var_length_diff = germline_variants.loc[germline_variants["Position"] == leftmost_var_pos, "MaxLengthDiff"].values[0]

                    if leftmost_var_length_diff > 0 and rightmost_red_pos - leftmost_var_pos == 0:
                        # rightmost pos immediately after closest insertion
                        # now check that lengths match
                        
                        if leftmost_var_length_diff >= rightmost_ref_idx:
                            # skip over full insertion at beginning of read 
                            idx = rightmost_ref_idx
                            continue

                        else:
                            # germline insertion plus more inserted bases, so add single error 
                            # could further identify which read bases match insertion, but fringe case so not going to figure that out now
                            is_error = True
                            self.sum_align_length += 1
                            # skip over full insertion at beginning of read 
                            idx = rightmost_ref_idx
                        
                    
                    else:
                        # read starts with non-germline insertion, so treat as usual
                        is_error = True
                        self.sum_align_length += 1
                        idx += 1
                        # rest of insertion will be handled in next iters

                # insertion not at beginning of read
                else:
                    is_error = True
                    self.sum_align_length += 1
                    idx += 1

            else:
                # check if at site of germline variant and skip over
                if reference_pos + 1 in all_var_positions:
                    var_length_diff = germline_variants.loc[germline_variants["Position"] == reference_pos + 1, "MaxLengthDiff"].values[0]
                    
                    # check for SNP or MNP
                    if var_length_diff == 0:
                        idx += germline_variants.loc[germline_variants["Position"] == reference_pos + 1, "RefLength"].values[0]
                        continue
                
                    # check for indels if current query pos isn't the last one
                    elif idx < len(use_aligned_pairs) - 1:

                        # check for deletion (next query pos None) or insertion (next ref pos None)
                        if var_length_diff != 0 and (use_aligned_pairs[idx + 1][0] == None or use_aligned_pairs[idx + 1][1] == None):
                            idx += abs(var_length_diff) + 1
                            continue
                        
                        # match or mismatch at same site as germline indel
                        else:
                            idx += 1
                            continue
                    
                    else:
                        # site of germline indel at the last base of the read, so can ignore
                        break
                    
                else:
                    # check for seq error mismatch
                    ref_base = self.sequence[reference_pos - self.coordinates[1]]
                    query_base = alignment.query_sequence[query_pos]

                    if ref_base != query_base:
                        is_error = True

                    self.sum_align_length += 1 
                    idx += 1

            self.sum_differences += int(is_error)


def run_error_analysis():

    debug_fh = None
    if args.debug_filename != "":
        debug_fh = open("read_error_position.tsv", 'w')
    save = pysam.set_verbosity(0)
    bam_file = pysam.AlignmentFile(args.bam, "rb")
    pysam.set_verbosity(save)

    for alignment in bam_file.fetch(until_eof=True):
        for r in regions:
            r.consume(alignment, debug_fh)

    # print results
    for r in regions:
        if r.sum_align_length != 0:
            print(f"{r.name}\t{r.sum_align_length}\t{r.sum_differences}\t{r.sum_differences/r.sum_align_length:.3}")
        else:
            print(f"{r.name}\t{r.sum_align_length}\t{r.sum_differences}\t{r.sum_differences}")


### arguments
parser = argparse.ArgumentParser()
parser.add_argument('--reference', type=str, required=True)
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--variants', type=str, required=True)
parser.add_argument('--chrom-name', type=str, default="chr5")
parser.add_argument('--split-region', type=bool, default=False)
parser.add_argument('--debug-filename', type=str, default="")

args = parser.parse_args()


### set variables
reference = pysam.FastaFile(args.reference)
chromosome_name = args.chrom_name

# PRDM9 coordinates
if args.split_region:
    left_region = ReferenceSegment("left", reference, (chromosome_name, 23516673, 23526672))
    znf_region = ReferenceSegment("znf", reference, (chromosome_name, 23526673, 23527764))
    right_region = ReferenceSegment("right", reference, (chromosome_name, 23527765, 23537764))
    regions = [ left_region, znf_region, right_region ]
else:
    region = ReferenceSegment("znf+10k", reference, (chromosome_name, 23516673, 23537764))
    regions = [ region ]

### prep germline variants table
germline_variants = pandas.read_csv(args.variants, delimiter="\t", header=None, names=["Position", "Ref", "Alt"])
# filter to PRDM9 coordinates and get variant lengths
germline_variants = (germline_variants.query('Position >= 23516673 and Position <=23537764')
  .assign(RefLength = germline_variants.Ref.str.len())
  .assign(MaxAltLength = germline_variants.Alt.apply(lambda x: max(len(alt) for alt in x.split(',')))))
germline_variants = germline_variants.assign(MaxLengthDiff = germline_variants.MaxAltLength - germline_variants.RefLength)
# merge different variants at same position and take the longest
germline_variants = (germline_variants.groupby('Position')
  .apply(lambda y: y.loc[y['MaxAltLength'].idxmax()]))
germline_variants.reset_index(drop=True, inplace=True)

run_error_analysis()























