#!/usr/bin/env python3

import csv
import argparse

desired_info_fields = ["CIEND", "CIPOS", "CHR2", "END", "MAPQ", "RE", "IMPRECISE", "PRECISE", "SVLEN", "SVMETHOD", "SVTYPE", "SUPP_VEC", "SUPP", "STRANDS"]

desired_format_fields = ["GT", "PSV", "LN", "DR", "ST", "QV", "TY", "ID", "RAL", "AAL", "CO"]

allowed_chromosomes = {str(i) for i in range(1, 23)} | {'X', 'Y'}

def parse_vcf_to_tsv(vcf_file, tsv_file):
    with open(vcf_file, 'r') as vcf, open(tsv_file, 'w', newline='') as tsv:
        tsv_writer = csv.writer(tsv, delimiter='\t')
        
        header_written = False
        sample_names = []

        for line in vcf:
            line = line.strip()
            if line.startswith("##"):
                continue  # <- Skip meta-information lines
            elif line.startswith("#CHROM"):
                # Parse the header line to get sample names
                header_parts = line.split('\t')
                sample_names = header_parts[9:]
                
                tsv_header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + desired_info_fields + sample_names
                tsv_writer.writerow(tsv_header)
                header_written = True
            else:
                if not header_written:
                    raise ValueError("VCF header not found before data lines.")
                
                parts = line.split('\t')
                chrom, pos, vid, ref, alt, qual, filt = parts[:7]
                if chrom not in allowed_chromosomes:
                    continue
                info = parts[7]
                format_keys = parts[8].split(':') if len(parts) > 8 else []
                sample_values = parts[9:] if len(parts) > 9 else []

                info_dict = {key: "." for key in desired_info_fields}
                for entry in info.split(';'):
                    if '=' in entry:
                        key, value = entry.split('=', 1)
                        if key in info_dict:
                            info_dict[key] = value
                    else:
                        if entry in info_dict:
                            info_dict[entry] = "."

                format_dicts = []
                for sample_value in sample_values:
                    format_dict = {key: "." for key in desired_format_fields}
                    values = sample_value.split(':')
                    for key, value in zip(format_keys, values):
                        if key in format_dict:
                            format_dict[key] = value
                    format_dicts.append(format_dict)

                row = [chrom, pos, vid, ref, alt, qual, filt] + [info_dict[key] for key in desired_info_fields]
                for format_dict in format_dicts:
                    row.extend(format_dict[key] for key in desired_format_fields)
                tsv_writer.writerow(row)

def load_chromosome_sizes(file_name):
    chromosome_sizes = {}
    with open(file_name, 'r') as file:
        for line in file:
            chrom, size = line.strip().split('\t')
            chromosome_sizes[chrom] = int(size)
    return chromosome_sizes

def parse_vcf_to_bed(vcf_file, bed_file, chr_size_file=None):
    chromosome_sizes = {}
    if chr_size_file:
        chromosome_sizes = load_chromosome_sizes(chr_size_file)

    with open(vcf_file, 'r') as vcf, open(bed_file, 'w', newline='') as bed:
        bed_writer = csv.writer(bed, delimiter='\t')

        for line in vcf:
            line = line.strip()
            if line.startswith("#"):
                continue
            else:
                parts = line.split('\t')
                chrom, pos, vid, ref, alt, qual, filt = parts[:7]
                if chrom not in allowed_chromosomes:
                    continue
                pos = int(pos)

                if pos >= 5001:
                    start = pos - 5000
                else:
                    start = 1

                end = pos + 5000
                
                if chrom in chromosome_sizes:
                    chr_size = chromosome_sizes[chrom]
                    if end > chr_size:
                        end = chr_size
                
                # Write the bed entry
                row = [chrom, start, end]
                bed_writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Parse VCF to TSV or BED format.")
    parser.add_argument("vcf_file", help="Path to the input VCF file")
    parser.add_argument("output_file", help="Path to the output file (TSV or BED)")
    parser.add_argument("output_type", choices=["tsv", "bed"], help="Type of the output file (tsv or bed)")
    parser.add_argument("--chr_size_file", help="Optional file with chromosome sizes", default=None)

    args = parser.parse_args()

    if args.output_type == "tsv":
        parse_vcf_to_tsv(args.vcf_file, args.output_file)
    elif args.output_type == "bed":
        parse_vcf_to_bed(args.vcf_file, args.output_file, args.chr_size_file)

if __name__ == "__main__":
    main()
