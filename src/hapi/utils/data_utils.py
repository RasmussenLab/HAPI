import csv
from pathlib import Path
import argparse
from typing import List, Tuple
import pandas as pd
import pysam


OpenFilesArgsOutput = Tuple[pysam.libcalignmentfile.AlignmentFile,
                            pysam.libcalignmentfile.AlignmentFile,
                            pysam.libcfaidx.FastaFile,
                            pysam.libcfaidx.FastaFile]

############## FILES OPENING FUNCTION DECLARATION ##############

def open_files_args(args: argparse.Namespace,
                    sample: str) -> OpenFilesArgsOutput:
    """
    Open files Bam and Fasta files based on CLI arguments and SNP names

    :param args: CLI script arguments
    :param sample: SNP name
    :return: bamvsref, BAM file of the GRCh37 
    :return: bamsvdel, BAM file of the GRCh37 with deletion 
    :return: fasta_ref, Fasta File of GRCh37
    :return: fasta_coll, Fasta File of GRCh37 with deletion
    """    
    bamvsref_path = Path.joinpath(args.folder_ref,
                                  sample + args.files_extension)

    bamvsdel_path = Path.joinpath(args.folder_coll,
                                  sample + args.files_extension)

    # Loading the bam file aligned vs the reference GRCh37 genome
    bamvsref = pysam.AlignmentFile(bamvsref_path, "rc",
                                   reference_filename=str(args.fasta_ref_file))

    # Loading the bam file aligned vs the coll reference 32del bamvsdel =
    # pysam.AlignmentFile(bamvsdel_file, "rc", reference_filename =
    # "/home/projects/cpr_10006/projects/ccr5/refs/CCR5_del32_120b.fasta")

    bamvsdel = pysam.AlignmentFile(bamvsdel_path, "rc", reference_filename=str(
        args.fasta_coll_file))

    # Loading the reference GRCh37 fasta file
    fasta_ref = pysam.FastaFile(args.fasta_ref_file)

    # Loading the coll GRCh37 fasta file
    fasta_coll = pysam.FastaFile(args.fasta_coll_file)

    return bamvsref, bamvsdel, fasta_ref, fasta_coll


def snp_haplo_list(snp_file: str) -> List[List[str]]:
    """Open file containing the list of the SNPs to analyze and save it in a
    list. The file was found at this path in computerome 2:
    /home/projects/cpr_10006/people/s162317/ancient/samtools_mpileup_LD
    /NyVikinSimon/nucleotideCount/ceu_haplotype86snps_nucleotides.txt """

    with open(snp_file, "r") as file:
        snp_list = [line.strip().split(sep="\t") for line in file]
    return snp_list


def dict_to_list(reads_dict: dict) -> list:
    """
    Convert the reads dictionary to list containing only the average minimum
    overlapping lengths without read names
    :param reads_dict: 
    :return: reads_list: list containing the minimum overlapping lengths of the
    reads from the dictionary
    """
    return [reads_dict[key] for key in sorted(reads_dict.keys())]


def write_probdf(prob_df: pd.DataFrame, outdir: Path, sample: str) -> None:
    """
    Function to write to file the probability dataframe of the 4 TOP SNPs
    along with their coverages etc

    :param prob_df: probability Dataframe
    :param outdir: Pathway to the folder where to save results
    :param sample: the SNP of which prob_df is saved
    """    

    prob_df.to_csv(Path.joinpath(outdir, sample + "top4SNPs_prob_df.tsv"),
                   sep="\t")


# I write a file containing the settings that I used to run the script
def write_settings(args: argparse.Namespace) -> None:
    """
    Write CLI arguments used in the script to the .tsv file 

    :param args: CLI arguments
    """    
    settings_dict = {arg: str(getattr(args, arg)) for arg in vars(args)}

    with open(args.output_folder / "settings.tsv", "w") as settings_file:
        writer = csv.writer(settings_file, delimiter="\t")
        for row in settings_dict.items():
            writer.writerow(row)


def write_results(results_filepath: Path, records: dict, header: bool) -> bool:
    """
    Write the results dict into the .tsv file 

    :param results_filepath: Pathway to the file where save results
    :param records: dictionary with records to save
    :param header: bool if use header when saving results
    :return: bool if use header when saving results next time
    """    
    if records:

        records_df = pd.DataFrame.from_records(records)

        if header:
            records_df.to_csv(results_filepath, sep="\t", header=header,
                              mode='w', index=False)
            return False

        else:
            records_df.to_csv(results_filepath, sep="\t", header=header,
                              mode='a', index=False)

    return header


def averaging_df_column(df: pd.DataFrame, cols_to_group: List[str],
                        avg_df_column: str) -> pd.DataFrame:
    """
    Average a column based on grouped columns

    :param df: pd.DataFrame of which column we want to average
    :param cols_to_group: a list of column names we want to group by
    :param avg_df_column: a column we want want to average
    :return: updated pd.DataFrame
    """
    df = df.assign(
        average_min_over=
        lambda x: x.groupby(cols_to_group)[avg_df_column].transform("mean"))

    df = df.drop(columns=[avg_df_column])
    df = df.drop_duplicates()

    return df
