import csv
from pathlib import Path

import pandas as pd
import pysam

############## FILES OPENING FUNCTION DECLARATION ##############

# TODO: change to regex to create a string containing the sample and then everything that is afterwards, but needs to end with cram or bam
def open_files_args(args, sample):

    bamvsref_path = Path.joinpath(args.folder_ref, sample + args.files_extension)

    bamvsdel_path = Path.joinpath(args.folder_fake, sample + args.files_extension)

    # Loading the bam file aligned vs the reference GRCh37 genome
    bamvsref = pysam.AlignmentFile(bamvsref_path, "rc", reference_filename = str(args.fasta_ref_file))

    # Loading the bam file aligned vs the fake reference 32del
    # bamvsdel = pysam.AlignmentFile(bamvsdel_file, "rc", reference_filename = "/home/projects/cpr_10006/projects/ccr5/refs/CCR5_del32_120b.fasta")

    bamvsdel = pysam.AlignmentFile(bamvsdel_path, "rc", reference_filename = str(args.fasta_fake_file))

    # Loading the reference GRCh37 fasta file
    fasta_ref = pysam.FastaFile(args.fasta_ref_file)

    # Loading the fake GRCh37 fasta file
    fasta_fake = pysam.FastaFile(args.fasta_fake_file)

    return bamvsref, bamvsdel, fasta_ref, fasta_fake


def snp_haplo_list(snp_file):
    """Open file containing the list of the SNPs to analyze and save it in a list. The file was found at this path in computerome 2:
    /home/projects/cpr_10006/people/s162317/ancient/samtools_mpileup_LD/NyVikinSimon/nucleotideCount/ceu_haplotype86snps_nucleotides.txt"""

    with open(snp_file, "r") as file:
        snp_list = [line.strip().split(sep="\t") for line in file]
    return snp_list

def dict_to_list(reads_dict):
    """
    Convert the reads dictionary to list containing only the average minimum overlapping lengths without read names
    :param reads_dict:
    :return: reads_list: list containing the minimum overlapping lengths of the reads from the dictionary
    """
    return [reads_dict[key] for key in sorted(reads_dict.keys())]

def write_probdf(prob_df, outdir, sample):
    """Function to write to file the probability dataframe of the 4 TOP SNPs along with their coverages etc"""    
    prob_df.to_csv(Path.joinpath(outdir, sample + "top4SNPs_prob_df.tsv"), sep="\t")

# I write a file containing the settings that I used to run the script
def write_settings():
    settings_dict = collections.OrderedDict()
    settings_dict["samples"] = str(args.samples)
    settings_dict["files-extension"] = args.files_extension
    settings_dict["folder-ref"] = args.folder_ref
    settings_dict["folder-fake"] = args.folder_fake
    settings_dict["fasta-ref"] = str(args.fasta_ref)
    settings_dict["fasta-fake"] = str(args.fasta_fake)
    settings_dict["snps"] = str(snp_file)
    settings_dict["length-threshold"] = str(length_threshold)
    settings_dict["overlapping-length-threshold"] = str(ol_threshold)
    settings_dict["perfect-match"] = str(perfect_match)
    settings_dict["output-folder"] = str(args.output_folder)


    with open(args.output_folder + "/settings_" + results_filename , "w") as settings_file:
        writer = csv.DictWriter(settings_file, fieldnames = ["Option","Chosen"], delimiter="\t")
        writer.writerows(settings_dict)


def write_results(results_filepath, records, header):
    
    if records:
        
        records_df = pd.DataFrame.from_records(records)
        
        if header:
            records_df.to_csv(results_filepath, sep="\t", header=header, mode='w', index=False)
            return False
            
        else:
            records_df.to_csv(results_filepath, sep="\t", header=header, mode='a', index=False)
    
    return header
    
        
        
    