#!/usr/bin/env python3

import pandas as pd
import sys
import os

def main():
    # Define input/output files
    input_counts = 'HLA_COUNTS.tsv'
    input_greps = 'HLA_GREPNG.tsv'
    output_file = 'imgt.tsv'

    # Check input files exist
    for file in [input_counts, input_greps]:
        if not os.path.isfile(file):
            print(f"Error: Input file '{file}' not found.")
            sys.exit(1)

    # Load data
    df_c = pd.read_csv(input_counts, sep='\t', header=None)
    df_g = pd.read_csv(input_greps, sep='\t', header=None)

    # Build map from greps
    grep_map = {row[0]: row[1] for _, row in df_g.iterrows()}

    # Collect values
    list_imgt = []
    list_hlas = []
    list_cont = []

    for _, row in df_c.iterrows():
        list_imgt.append(row[0])
        list_cont.append(row[1])
        list_hlas.append(grep_map.get(row[0]))

    # Build output DataFrame
    df_out = pd.DataFrame({
        'IMGTID': list_imgt,
        'HLA_ID': list_hlas,
        'COUNTS': list_cont    
    })

    # Sort and reset index
    df_out = df_out.sort_values(by='HLA_ID', key=lambda col: col.astype(str))
    df_out.reset_index(drop=True, inplace=True)

    # Write to output
    df_out.to_csv(output_file, sep='\t', index=False)
    print(f"Output written to: {output_file}")

if __name__ == "__main__":
    main()
