from __future__ import print_function
from __future__ import division

import matplotlib.pyplot as plt
import matplotlib_venn

import pandas as pd

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("A", help="First table to intersect")
    parser.add_argument("B", help="Second table to intersect")
    parser.add_argument("plot_file")
    parser.add_argument("--input_mode", choices=["tab", "vcf"], default="tab", help="Input file type (Default: %(default)s)")
    parser.add_argument("--A_label", default="A", help="Label for first table")
    parser.add_argument("--B_label", default="B", help="Label for second table")
    parser.add_argument("--columns", nargs="+", default=["CHROM", "POS", "REF", "ALT"], help="List of column names to intersect (Default: %(default)s)")

    args = parser.parse_args()

    if args.input_mode == "tab":
        A = pd.read_table(args.A)
        B = pd.read_table(args.B)
    else:
        cols = ["CHROM", "POS", "ID", "REF", "ALT"]
        A = pd.read_table(args.A, comment="#", names=cols, usecols=cols, header=None)
        B = pd.read_table(args.B, comment="#", names=cols, usecols=cols, header=None)

    intersect = A.merge(B, how="inner", on=args.columns)

    subsets = (
                A.shape[0] - intersect.shape[0],
                B.shape[0] - intersect.shape[0],
                intersect.shape[0],
              )

    ven = matplotlib_venn.venn2(subsets=subsets, set_labels=(args.A_label, args.B_label))
    c = matplotlib_venn.venn2_circles(subsets=subsets)
    plt.savefig(args.plot_file)
