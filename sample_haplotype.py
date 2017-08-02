from __future__ import print_function
from __future__ import division

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import sys
import pandas as pd

import pysam
import random
import argparse

def get_allele_seq(record, indiv, fixed=None):
    """Take a pysam VariantFile record and an individual
    and return the sequence of a random allele from that indiv.
    Also can take a DataFrame fixed with CHROM, POS, and ALT columns
    to force ALT alleles for particular records."""

    if fixed is not None:
        fix = fixed.loc[(fixed.CHROM == record.chrom) & (fixed.POS == record.pos) & (fixed.ALT == record.alts[0])]
        if fix.shape[0] == 1:
            return record.alts[0]

    gt = record.samples[indiv]["GT"]
    allele_num = gt[random.randrange(len(gt))]
    if allele_num is None or allele_num == 0:
        print(record.pos, record.alts, gt)
        return record.ref
    else:
        return record.alts[allele_num - 1]

def update_seq(sequence, pos, ref, allele, region):
    if ref == allele:
        return sequence
    return sequence[:pos - region[1] - 1] + allele + sequence[pos + len(ref) - region[1] - 1:]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("indiv")
    parser.add_argument("reference_fasta")
    parser.add_argument("outfasta")
    parser.add_argument("region", help="Region chr:start-end")
    parser.add_argument("--fixed_variants", help="Table of CHROM, POS, REF, and ALT for variants to force the ALT allele")

    args = parser.parse_args()

    region = args.region.replace(":", " ").replace("-", " ").split()
    region[1] = int(region[1])
    region[2] = int(region[2])

    fixed = pd.DataFrame(columns=["CHROM", "POS", "ALT"])
    if args.fixed_variants is not None:
        fixed_vcf = pysam.VariantFile(args.fixed_variants)
        for i, record in enumerate(fixed_vcf.fetch(*region)):
            fixed.loc[i] = [record.chrom, record.pos, record.alts[0]]
    sequence = pysam.FastaFile(args.reference_fasta).fetch(*region)
    vcf_reader = pysam.VariantFile(args.vcf)
    for variant in vcf_reader.fetch(*region):
        allele = get_allele_seq(variant, args.indiv, fixed)
        sequence = update_seq(sequence, variant.pos, variant.ref, allele, region)
    seq_record = SeqRecord(Seq(sequence), id=args.indiv, description="")
    SeqIO.write([seq_record], args.outfasta, "fasta")
