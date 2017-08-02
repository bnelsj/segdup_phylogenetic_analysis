from __future__ import print_function
from __future__ import division

import vcf

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_vcf")
    parser.add_argument("output_vcf")

    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.input_vcf, "r"))
    vcf_writer = vcf.Writer(open(args.output_vcf, "w"), vcf_reader)
    for record in vcf_reader:
        no_call = [sn["GT"] for sn in record.samples if sn["GT"] in ["././.", "./././."]]
        if len(no_call) > 0:
            continue
        if all([sn["GT"] in ["0/0/1", "././.", "0/0/1/1", "./././."] for sn in record.samples]) or all([sn["GT"] in ["1/1/1", "././.", "1/1/1/1", "./././."] for sn in record.samples]):
            vcf_writer.write_record(record)
