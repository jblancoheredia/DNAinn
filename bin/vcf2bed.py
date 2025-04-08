import sys
import argparse

def vcf_to_bed(vcf_file, bed_file):
    with open(vcf_file, "r") as vcf, open(bed_file, "w") as bed:
        for line in vcf:
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 5:
                print(f"Skipping malformed line: {line.strip()}", file=sys.stderr)
                continue
            
            chrom, pos, _, ref, alt = fields[:5]

            start = int(pos) - 1
            end = start + max(len(ref), len(alt))
            
            start_expanded = max(0, start - 1000)
            end_expanded = end + 1000

            bed.write(f"{chrom}\t{start_expanded}\t{end_expanded}\n")

            if end != start + 1:
                start_expanded_end = max(0, end - 1000)
                end_expanded_end = end + 1000
                bed.write(f"{chrom}\t{start_expanded_end}\t{end_expanded_end}\n")

    print(f"BED file saved: {bed_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a VCF file to a BED file with expanded regions for GATK4_BEDTOINTERVALLIST.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    args = parser.parse_args()

    vcf_to_bed(args.input, args.output)
