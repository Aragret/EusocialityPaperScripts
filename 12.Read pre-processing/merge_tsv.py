# -*- coding: utf-8 -*-

"""
Created on Fri May 31 06:26:33 2024.

@author: caumont
"""

import os
import argparse
import pandas as pd


def merging(folder, outfile):
    """Merge multiple tsv file into a single one on gene name."""
    file_list = [
        os.path.join(folder, cfile)
        for cfile in os.listdir(folder)
        if cfile.endswith(".tsv")
    ]

    df_merge = pd.DataFrame()
    df_merge["gene_count"] = pd.read_table(file_list[0],
                                           sep="\t",
                                           header=None)[0]
    for file in file_list:
        header = file.split("/")[1].strip(".count.tsv")
        print(header)
        df_txt = pd.read_table(file,
                               sep="\t",
                               header=None,
                               names=["gene_count", header])

        df_merge = df_merge.merge(df_txt,
                                  on="gene_count",
                                  how="outer")

        df_final = df_merge.drop([df_merge.index[-1],df_merge.index[-2],df_merge.index[-3],df_merge.index[-4],df_merge.index[-5]])

    with open(outfile, "w") as of:
        df_final.to_csv(of, sep="\t", index=False)


def main():
    """Start main function of the program."""
    parser = argparse.ArgumentParser(
        description="merge the sample counts into a species count"
    )

    parser.add_argument(
        "-f",
        "--finput",
        help="The input folder path containing the gene_count.tsv files.",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-o",
        "--output",
        help="The ouput species.readcount.tsv file, path included.",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    with open(args.output, "w"):
        print()

    merging(
        args.finput,
        args.output,
    )


if __name__ == "__main__":
    main()
