""" The script is divided in two parts:
A) CEU rs333 haplotype probability calculation

For each SNP of the snp_file, I extract the pileup of their position and calculate the probability of each possible
genotype (Ref/Ref, Ref/Alt, Alt/Alt) given the data, i.e. p(G|D), where:
- G stands for Genotype of the alleles
- D stands for Data.

After calculating it for each SNP, I do a weighted average/multiplication across all the SNPs.
In this way, for each sample I'll have the probability of the haplotype.
I'll use this as a Prior probability, i.e. p(D|G) for the calculation of the next step

B) Deletion and Reference sequence probability calculation

Each ancient sample DNA has been aligned against the reference genome GRCh37 and against a fake reference containing the
32 bp deletion in the ccr5 gene. For both of these bams I want to calculate the probability of having the deletion and
of having the reference sequence, following the genotype strategy of point A. In particular, I'll calculate the probability
of each possible "genotype"  (Ref/Ref, Ref/Del, Del/Del) given the data, i.e. p(G|D), where:
- G stands for Genotype of the Reference or Deleted sequence
- D stands for Data

p(G|D) = p(G) p(D|G) / p(D)

Note: Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.
Here I converted the 0-based coordinates to 1-based coordinates, simply adding +1 to the 0-based coordinate.
E.g. of these variables:
reference_start: converted with +1
query_position: converted with +1
reference_end: not converted because it's in 0-based exclusive
"""

# Libraries loading
import sys
import numpy as np
from numpy import prod
import pandas as pd
import pysam
import math
import copy
from time import time
import os
import collections
from collections import defaultdict
from collections import OrderedDict
from statistics import mean
import argparse
import re
import csv




# Start time initiation
start = time()

# ArgParse to parse the arguments given in the command-line
parser = argparse.ArgumentParser()

parser.add_argument("--samples", required=True, help = "File containing list of samples to process")
# parser.add_argument("--bam-ref", required=True,
#                     help="write the input bam file aligned against the GRCH37 reference genome")
# parser.add_argument("--bam-del", required=True, help="write the input bam file aligned against the 32del fake genome")

parser.add_argument("--files-extension", required=True, help = "String describing the extension of the file, e.g. .rmdup.realign.md.cram or .rmdup.realign.bam")
parser.add_argument("--folder-ref", required=True,
                    help="Folder path of the bam files aligned against the GRCH37 reference genome")
parser.add_argument("--folder-fake", required=True,
                    help="Folder path of the bam files aligned against the 32del fake genome")


parser.add_argument("--fasta-ref", required=True, help="fasta file containing the reference genome")
parser.add_argument("--fasta-fake", required=True, help="fasta file containing the fake reference genome")

parser.add_argument("--snps", required=True, help="text file containing list of the 4 top SNPs")
parser.add_argument("--haplotype", required=False, help="text file containing list of the 86 SNPs")
parser.add_argument("--output-folder", required=True,
                    help="write the output folder in which to append the results of the probability calculations")
parser.add_argument("--baq-snps", choices=["yes", "no"], default="no", required=False)
parser.add_argument("--baq-deletion", choices=["yes", "no"], default="no", required=False)

parser.add_argument("--length-threshold", required=True)
parser.add_argument("--overlapping-length-threshold", default=4)
parser.add_argument("--perfect-match", required=True)
parser.add_argument("--adjustment-threshold", required=True)
# baq_snp = "no"
# baq_deletion = "no"
args = parser.parse_args()

############## FILES OPENING FUNCTION DECLARATION ##############

def samples_to_list():
    with open(args.samples, "r") as samples_file:
        samples_list = samples_file.read().splitlines()

    return samples_list


# TODO: change to regex to create a string containing the sample and then everything that is afterwards, but needs to end with cram or bam
def open_files_args(sample):
    bamvsref_file = str(args.folder_ref) + "/" + sample + str(args.files_extension)

    bamvsdel_file = str(args.folder_fake) + "/" + sample + str(args.files_extension)

    # Loading the bam file aligned vs the reference GRCh37 genome
    bamvsref = pysam.AlignmentFile(bamvsref_file, "rc", reference_filename = str(args.fasta_ref))

    # Loading the bam file aligned vs the fake reference 32del
    # bamvsdel = pysam.AlignmentFile(bamvsdel_file, "rc", reference_filename = "/home/projects/cpr_10006/projects/ccr5/refs/CCR5_del32_120b.fasta")

    bamvsdel = pysam.AlignmentFile(bamvsdel_file, "rc", reference_filename = str(args.fasta_fake))

    # Loading the reference GRCh37 fasta file
    fasta_ref = pysam.FastaFile(args.fasta_ref)

    # Loading the fake GRCh37 fasta file
    fasta_fake = pysam.FastaFile(args.fasta_fake)

    # Loading the file containing the list of SNPs to use to calculate the haplotype probability
    snp_file = args.snps

    # Loading the file containing the list of the 86 SNPs that form the haplotype

    haplotype_file = args.haplotype

    # if args.haplotype:
    #     print("--haplotype given")

    #     haplotype_file = args.haplotype

    # else:

    #     print("--haplotype not given")
    #     haplotype_file = None

    baq_snp = args.baq_snps

    baq_deletion = args.baq_deletion

    length_threshold = int(args.length_threshold)

    ol_threshold = int(args.overlapping_length_threshold)

    perfect_match = str(args.perfect_match)

    adjustment_threshold = int(args.adjustment_threshold)

    return bamvsref, bamvsdel, fasta_ref, fasta_fake, snp_file, haplotype_file, baq_snp, baq_deletion, length_threshold, ol_threshold, perfect_match, adjustment_threshold


def open_files():
    # Loading the bam file aligned vs the reference GRCh37 genome
    bamvsref_file = pysam.AlignmentFile(
        "/Users/lmz306/OneDrive - Københavns Universitet/delta_ccr5/analyses/del_neg_id/VK548.sort.rmdup.realign.md.bam", "rb")

    # Loading the bam file aligned vs the fake reference 32del
    bamvsdel_file = pysam.AlignmentFile(
        "/Users/lmz306/OneDrive - Københavns Universitet/delta_ccr5/analyses/del_neg_id/library_VK548.bam", "rb")

    # Loading the reference GRCh37 fasta file
    fasta_ref = pysam.FastaFile("/Users/lmz306/OneDrive - Københavns Universitet/delta_ccr5/refs/hs.build37.1.fa")

    # Loading the fake GRCh37 fasta file
    fasta_fake = pysam.Fastafile(
        "/Users/lmz306/OneDrive - Københavns Universitet/delta_ccr5/refs/CCR5_del32_120b.fasta")

    # Loading the file containing the list of SNPs to use to calculate the haplotype probability
    snp_file = "/Users/lmz306/OneDrive - Københavns Universitet/delta_ccr5/analyses/del_neg_id/ceu_haplo_82snps.txt"

    baq_snp = "no"

    baq_deletion = "no"
    return bamvsref_file, bamvsdel_file, fasta_ref, fasta_fake, snp_file, baq_snp, baq_deletion


############## Part A FUNCTIONS DECLARATION ##############
def phred2prob(x):
    """Convert Phred score to probability"""
    return 10 ** (-x / 10)


def prob2phred(x):
    """Convert probability to Phred score"""
    return -10 * math.log10(x)


def snp_haplo_list(snp_file):
    """Open file containing the list of the SNPs to analyze and save it in a list. The file was found at this path in computerome 2:
    /home/projects/cpr_10006/people/s162317/ancient/samtools_mpileup_LD/NyVikinSimon/nucleotideCount/ceu_haplotype86snps_nucleotides.txt"""

    with open(snp_file, "r") as file:
        snp_list = []
        for line in file:
            snp_list.append(line.strip().split(sep="\t"))
    return snp_list


# def extr_rbases_bam(bamvsref_file, chrom, coordinate, ref, alt, baq, adjustment_threshold, min_base_quality=30,
#                     min_mapping_quality=30):

# Trying to do the same as Kirstine
def extr_rbases_bam(bamvsref_file, chrom, coordinate, ref, alt, baq, adjustment_threshold, min_base_quality=0,
                    min_mapping_quality=0):
    """
    Extract the read bases, i.e. the bases of the reads that map to a specific position in the bam file.
    bq and mq are base quality and mapping quality. Change them if I want to filter.
    :param bamvsref_file: bam file aligned vs the GRCh37 reference; pysam.AlignmentFile
    :param chrom: chromosome number of the position to extract; string
    :param coordinate: coordinate of the position to extract; integer
    :param ref: reference allele at the position; string
    :param alt: alternate allele at the position; string
    :param min_base_quality: base quality; integer
    :param adjust_capq_threshold:
    :param min_mapping_quality: mapping quality; integer
    :param baq: base alignment quality filtering; string
    :return:
    """

    # List containing the reads overlapping to the position when the base found corresponds to either the Ref or the Alt
    reads_list = []

    # If the base found is different from either the Ref or the Alt, I'll put the read in this other list
    other_list = []

    ref_list, alt_list = [], []
    if baq == "no":
        for pileupcolumn in bamvsref_file.pileup(chrom, coordinate - 1, coordinate, truncate=True,
                                                 min_base_quality=min_base_quality,
                                                 adjust_capq_threshold=adjustment_threshold,
                                                 min_mapping_quality=min_mapping_quality):
            reads_list, other_list, ref_list, alt_list = extract_lists(pileupcolumn, ref, alt)

            


    elif baq == "yes":
        for pileupcolumn in bamvsref_file.pileup(chrom, coordinate - 1, coordinate, truncate=True, stepper="samtools",
                                                 fastafile=fasta_ref, compute_baq=True,
                                                 min_base_quality=min_base_quality,
                                                 adjust_capq_threshold=adjustment_threshold,
                                                 min_mapping_quality=min_mapping_quality):

            reads_list, other_list, ref_list, alt_list = extract_lists(pileupcolumn, ref, alt)

    return reads_list, other_list, ref_list, alt_list


def extract_lists(pileupcolumn, ref, alt):
    """
    Function to extract the reads mapping to a certain position from the PileupColumn object. If the read base at the
    position corresponds to either the Reference or Alternate allele of the SNP as written in the file, the read will be
    put in the object "reads_list". If not, the read will be put in the "other_list"
    :param pileupcolumn:
    :return: reads_list
    :return: other_list
    """
    # I initialize the empty lists
    reads_list, other_list, ref_list, alt_list = [], [], [], []

    # For each read mapping to the PileupColumn position
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:

            # Uncomment if want to print the coverage at the position
            # print('\tbase in read %s = %s' %
            #       (pileupread.alignment.query_name,
            #        pileupread.alignment.query_sequence[pileupread.query_position]))

            # I save the read name, the read base that map to the position, and the base quality
            read_name = pileupread.alignment.query_name
            base = pileupread.alignment.query_sequence[pileupread.query_position]

            quality = pileupread.alignment.query_alignment_qualities[pileupread.query_position]

            read_length = pileupread.alignment.query_length

            if read_length <= length_threshold:
            # if True:
                # If the read base corresponds to the Reference or to the Alternate allele, put it in the reads_list
                if base == ref or base == alt:
                    reads_list.append([base, quality, read_name, read_length])

                if base == ref:
                    ref_list.append([base, quality, read_name, read_length])

                if base == alt:
                    alt_list.append([base, quality, read_name, read_length])

                # If the read base is not like the reference nor like the alternate allele, put it in the other_list
                if base != ref and base != alt:
                    other_list.append([base, quality, read_name, read_length])

    return reads_list, other_list, ref_list, alt_list


# SNPs genotype calculation
# pD_G = p(D|G) likelihood
# pG_D = p(G|D)
# pD_ = p(D)
# p(G|D) = p(G) * p(D|G) / p(D)

# p(D|G)
def pD_G_(list, ref, alt):
    """
    Calculate step 1, 2, and 3 of genotype likelihood for each SNP
    :param list: list of read bases mapping each of the CEU rs33 haplotype SNPs; list
    :param ref: reference allele; string
    :param alt: alternate allele; string
    :return: pD_RR, pD_RA, pD_AA: genotypes likelihoods, i.e. probability of the data given each genotype; float
    """

    # Here with RR, RA, and AA I actually mean RR,RD,DD, so the deletion genotypes
    # p(D|RR), p(D|RA), p(D|AA)
    pD_RR, pD_RA, pD_AA = 0, 0, 0

    # If the list is not empty
    if not list == []:
        # Iterate over each read base of the list
        for base in list:

            # Following Simon's slides sequence:
            # Step 1: calculate prob of each base given the observed allele, given that the true base is the observed
            # p_b_A = p(b|A)
            p_b_A_main = (1 - phred2prob(base[1]))
            p_b_A_rest = (phred2prob(base[1])) / 3

            # Step 2: calculate prob of each observed base given EACH POSSIBLE genotype, so Ref/Ref, Ref/Alt, Alt/Alt.
            # Step 3: multiply (in this case sum, since they are logarithms), over all bases to calculate the Likelihood
            #         of each genotype p(D|G)
            if base[0] == ref:
                pD_RR += math.log10(p_b_A_main)  # Would be (p_b_A_main + p_b_A_main)/2, so it simplifies to p_b_A_main
                pD_RA += math.log10((p_b_A_main + p_b_A_rest) / 2)
                pD_AA += math.log10(p_b_A_rest)

            elif base[0] == alt:
                pD_RR += math.log10(p_b_A_rest)
                pD_RA += math.log10((p_b_A_main + p_b_A_rest) / 2)
                pD_AA += math.log10(p_b_A_main)

    else:
        # If the reads_list is empty, i.e. there are no reads mapping that position, set the probability of the data
        # given the genotype (the genotype likelihood) as 0.33, i.e. as random
        pD_RR = math.log10(0.33)
        pD_RA = math.log10(0.33)
        pD_AA = math.log10(0.33)

    return 10 ** pD_RR, 10 ** pD_RA, 10 ** pD_AA


# p(D)
def pD_(pD_RR, pD_RA, pD_AA):
    """
    p(D) = SUM p(Gi) p(D|Gi), for each genotype RR, RA, AA. Here we assume a uniform prior distribution, giving it
    a value of 0.33
    :param pD_RR: likelihood of getting the observed Data given the Genotype Ref Ref, calculated in the function pD_G_ ; float
    :param pD_RA: likelihood of getting the observed Data given the Genotype Ref Alt, calculated in the function pD_G_ ; float
    :param pD_AA: likelihood of getting the observed Data given the Genotype Alt Alt, calculated in the function pD_G_ ; float
    :return: pD: p(D)= SUM p(Gi) p(D|Gi)
    """

    pD = (0.33 * pD_RR) + (0.33 * pD_RA) + (0.33 * pD_AA)
    return pD


# p(G|D)
def pG_D_(pD_RR, pD_RA, pD_AA, pD, pG=0.33):
    """
    Posterior probability of the specific base
    :param pD_RR: likelihood of getting the observed Data given the Genotype Ref Ref, calculated in the function pD_G_ ; float
    :param pD_RA: likelihood of getting the observed Data given the Genotype Ref Alt, calculated in the function pD_G_ ; float
    :param pD_AA: likelihood of getting the observed Data given the Genotype Alt Alt, calculated in the function pD_G_ ; flo
    :param pD: p(D)
    :param pG: prior probability. Here we assume a uniform prior, so 0.33
    :return: pRR_D, pRA_D, pAA_D
    """
    # These are the probabilities of each genotype given the data FOR THIS BASE.
    pRR_D = (pD_RR * pG) / pD
    pRA_D = (pD_RA * pG) / pD
    pAA_D = (pD_AA * pG) / pD


    return pRR_D, pRA_D, pAA_D


def prob_to_weighted(probability_dataframe):
    """
    Add, to the dataframe prob_df, columns containing the probabilities of each SNP multiplied by each relative R squared value
    :param probability_dataframe:
    :return: prob_df
    """

    prob_df["P(RR|D)w"] = prob_df["P(RR|D)"] * prob_df["rsquared"]
    prob_df["P(RA|D)w"] = prob_df["P(RA|D)"] * prob_df["rsquared"]
    prob_df["P(AA|D)w"] = prob_df["P(AA|D)"] * prob_df["rsquared"]

    # Convert to logarithms
    prob_df["logP(RR|D)w"] = np.log10(prob_df["P(RR|D)w"])
    prob_df["logP(RA|D)w"] = np.log10(prob_df["P(RA|D)w"])
    prob_df["logP(AA|D)w"] = np.log10(prob_df["P(AA|D)w"])

    return prob_df

def coverage_dict(dict_snps_cov, reads_list, id, other_list, ref_list, alt_list):
    """
    Updates the dictionary counter "dict_snps_cov" with the number of reads overlapping the 4 TOP SNPs
    :param dict_snps_cov:
    :return:
    """

    if id == "rs113341849":
        dict_snps_cov["rs113341849"] = len(alt_list)

    elif id == "rs113010081":
        dict_snps_cov["rs113010081"] = len(alt_list)

    elif id == "rs11574435":
        dict_snps_cov["rs11574435"] = len(alt_list)

    elif id == "rs79815064":
        dict_snps_cov["rs79815064"] = len(alt_list)

    return dict_snps_cov


def calc_snps_posteriors(snp_list):
    """
    Calculate the Posterior Probabilities for each SNP, storing the information in a dataframe and saving the statistics
    relative to the total coverage on the SNPs and to the coverage of each of the 4 TOP SNPs

    :return: prob_df, dataframe containing the probabilities
    :return: coverage_ref and coverage_alt, counter containing the total number of bases mapping to all the SNPs, a total coverage
    :return: dict_snps_cov, dictionary containing the coverage for each of the 4 TOP SNPs
    """

    # I initialize an empty dataframe where I'll store the Posterior probabilities of each genotype at each position
    columns_df = ["id", "chrom", "ref", "alt", "rsquared", "other", "n_total_reads", "n_ref_reads","n_alt_reads", "P(RR|D)", "P(RA|D)", "P(AA|D)"]
    prob_df = pd.DataFrame(columns=columns_df)

    # I initialize an empty dictionary and an empty counter to save the coverage statistics of the SNPs
    coverage_ref, coverage_alt, coverage_other = 0, 0, 0
    dict_snps_cov = {}

    # Iterate through each SNP 
    for snp in range(len(snp_list)):

        # Extract the different fields of the file
        id = snp_list[snp][0]

        chrom = str(3)
        coordinate = int(snp_list[snp][1])
        ref = snp_list[snp][2]
        alt = snp_list[snp][3]
        rsquared = snp_list[snp][4]

        # 1 - Extract all the reads bases mapping to each snp
        reads_list, other_list, ref_list, alt_list = extr_rbases_bam(bamvsref, chrom, coordinate, ref, alt, baq_snp, adjustment_threshold)

        # 2 - Update SNPs coverage statistics
        dict_snps_cov = coverage_dict(dict_snps_cov, reads_list, id, other_list, ref_list, alt_list)

        coverage_ref += len(ref_list)
        coverage_alt += len(alt_list)
        coverage_other += len(other_list)

        # 3 - Calculate p(D|G) for each genotype, so for Ref/Ref (pD_RR), Ref/Del (pD_RA), Del/Del (pD_AA)
        pD_RR, pD_RA, pD_AA = pD_G_(reads_list, ref, alt)

        # 4 - Calculate p(D), i.e. denominator of the final equation. 
        pD = pD_(pD_RR, pD_RA, pD_AA)

        # 5 - Calculate p(G|D), the posterior probability of each deletion genotype. 
        pRR_D, pRA_D, pAA_D = pG_D_(pD_RR, pD_RA, pD_AA, pD)

        # If reads_list is empty it means that there were no reads overlapping the position, so I'll put 0.33
        if not reads_list:
            row = pd.DataFrame([[id, chrom, ref, alt, float(rsquared), 0, 0, 0, 0, 0.33, 0.33, 0.33]],
                               columns=columns_df,
                               index=[coordinate])
        else:
            row = pd.DataFrame([[id, chrom, ref, alt, float(rsquared), 0, len(reads_list), len(ref_list), len(alt_list), pRR_D, pRA_D, pAA_D]],
                               columns=columns_df, index=[coordinate])

        # Whatever is the case, I'll append the row with this SNP to the probability dataframe
        prob_df = prob_df.append(row)

        # If there are reads with a base other than the REF or the ALT, add them as a list in the column "other"
        if not other_list == []:
            prob_df = prob_df.astype({"other": object})

            prob_df.at[coordinate, "other"] = other_list


    return prob_df, coverage_ref, coverage_alt, coverage_other, dict_snps_cov





def calc_prob_joint(prob_df):
    """
    Calculate Joint Posterior Probabilities of all the SNPs and normalize them
    :param prob_df:
    :return: pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm
    """

    log_list = np.array([prob_df["logP(RR|D)w"].sum(), prob_df["logP(RA|D)w"].sum(), prob_df["logP(AA|D)w"].sum()])

    # I calculate the maximum
    maximum = np.float64(max(log_list))

    # I substract the maximum 
    max_substracted = log_list - maximum

    # I exponentiate the logarithms
    exponentiated = 10 ** (max_substracted)

    # I sum them
    sum_exponents = np.sum(exponentiated)

    # I divide each value by the sum
    normalized = np.true_divide(exponentiated, sum_exponents)

    # These will be used as the prior probability in the formula to calculate the Posterior of the deletion and the reference
    pRR_D_joint_norm = normalized[0]
    pRA_D_joint_norm = normalized[1]
    pAA_D_joint_norm = normalized[2]

    return pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm


############## Part B FUNCTIONS DECLARATION ##############

def minimum_overlap(bam_file, chrom, position_list, adjustment_threshold, df_mapping_all, baq, min_base_quality=30, min_mapping_quality=30):
    """
    Function to calculate the minimum length that a read overlaps the:
    - Starting and Ending position of the 32deletion in the bam aligned against the reference genome GRCh37
      --> To find reads NOT HAVING the deletion, i.e. having the reference 32bp sequence
    - Position breakpoint of the 32deletion in the bam aligned against the fake 32del reference
      --> To find reads HAVING the deletion, i.e. not having the reference 32bp sequence

    :param bam_file: bam file; pysam.AlignmentFile
    :param chrom: chromosome number; string
    :param position_list: list containing the 4 different Starting and Ending coordinates or the Position coordinate; list
    :param min_base_quality: base quality; integer
    :param adjust_capq_threshold:
    :param min_mapping_quality: mapping quality; integer
    :param baq: base alignment quality filtering; string
    :return: reads_dict: containing the dictionary of the reads with their minimum overlaps; dict

    """
    lengths_dict = {}
    nm_tags_dict = {}
    # If the list contains 4 elements --> Bam aligned vs Reference, to detect reads NOT having the deletion
    if len(position_list) == 4:
        # I initialize an empty dictionary in which the values of the keys will be of type list
        reads_dict = defaultdict(list)
        
        # For each coordinates couples
        for start_end in position_list:

            # I save the starting and ending coordinates of the 32bp sequence
            S = start_end[0]
            E = start_end[1]

            # For each of the S and E
            for pos in start_end:
                if baq == "no":
                    for pileupcolumn in bam_file.pileup(chrom, pos - 1, pos, truncate=True,
                                                        min_base_quality=min_base_quality,
                                                        adjust_capq_threshold=adjustment_threshold,
                                                        min_mapping_quality=min_mapping_quality):
                        reads_dict, lengths_dict, df_mapping_all, nm_tags_dict = min_over_reference(pileupcolumn, S, E, pos, reads_dict, lengths_dict, df_mapping_all, nm_tags_dict)

                elif baq == "yes":
                    for pileupcolumn in bam_file.pileup(chrom, pos - 1, pos, truncate=True, stepper="samtools",
                                                        fastafile=fasta_ref, compute_baq=True,
                                                        min_base_quality=min_base_quality,
                                                        adjust_capq_threshold=adjustment_threshold,
                                                        min_mapping_quality=min_mapping_quality):
                        reads_dict, lengths_dict, df_mapping_all, nm_tags_dict = min_over_reference(pileupcolumn, S, E, pos, reads_dict, lengths_dict, df_mapping_all, nm_tags_dict)

        return reads_dict, lengths_dict, df_mapping_all, nm_tags_dict

    # If the list contains 1 element, i.e. P --> Bam aligned vs fake Reference, to detect reads HAVING the deletion
    elif len(position_list) == 1:

        reads_dict = {}
        if baq == "no":
            for pileupcolumn in bam_file.pileup(chrom, position_list[0] - 1, position_list[0], truncate=True,
                                                min_base_quality=min_base_quality,
                                                adjust_capq_threshold=adjustment_threshold,
                                                min_mapping_quality=min_mapping_quality):
                reads_dict, lengths_dict, df_mapping_all, nm_tags_dict = min_over_32del(pileupcolumn, position_list, lengths_dict, df_mapping_all, nm_tags_dict)


        elif baq == "yes":
            for pileupcolumn in bam_file.pileup(chrom, position_list[0] - 1, position_list[0], truncate=True,
                                                stepper="samtools", fastafile=fasta_fake, compute_baq=True,
                                                min_base_quality=min_base_quality,
                                                adjust_capq_threshold=adjustment_threshold,
                                                min_mapping_quality=min_mapping_quality):
                reads_dict, lengths_dict, df_mapping_all, nm_tags_dict = min_over_32del(pileupcolumn, position_list, lengths_dict, df_mapping_all, nm_tags_dict)

        return reads_dict, lengths_dict, df_mapping_all, nm_tags_dict


def min_over_reference(pileupcolumn, S, E, pos, reads_dict, lengths_dict, df_mapping_all, nm_tags_dict):
    """
    Calculate the minimum overlapping length for each read in the bam file mapped against the GRCh37
    :param pileupcolumn:
    :param S:
    :param E:
    :param pos:
    :param reads_dict:
    :return:
    """

    for pileupread in pileupcolumn.pileups:
        
        # If the read is not a deletion
        if not pileupread.is_del and not pileupread.is_refskip:

            read_name = pileupread.alignment.query_name

            # Position of the starting base of the read on the genome
            reference_start = pileupread.alignment.reference_start + 1

            # Position of the ending base of the read on the genome
            reference_end = pileupread.alignment.reference_end

            # Position in the read overlapping the pileup position
            query_position = pileupread.query_position + 1

            # print("query position:", query_position)
            read_length = pileupread.alignment.query_length

            read_sequence = pileupread.alignment.query_sequence

            nm_tag = pileupread.alignment.get_tag("NM")




            # I filter for only the reads under the threshold that I set (80 bp)
            if read_length <= length_threshold:
   
                # If the starting and ending position of the read are to the left and right of the deletion,
                # I'll have a read overlapping across both the S and E.
                if reference_start <= S and reference_end >= E:
                    # Minimum overlapping length is 32
                    min_over = 32
                else:
                    # I calculate the left and right overlap of the read
                    left = query_position
                    right = reference_end - pos + 1

                    # The minimum overlapping length is the minimum between these two
                    min_over = min(left, right)


                # I save all in a row
                row_to_add = [sample, read_name, reference_start, reference_end, read_sequence, read_length, min_over, nm_tag, "ref"]

                # print("row_to_add", row_to_add)
                df_length = len(df_mapping_all)
                
                # print("df_length", df_length)
                # That I add to a dataframe so I can analyse it afterwards
                # df_mapping_all.loc[df_length] = row_to_add

                if min_over >= ol_threshold:
                    # I append the minimum overlapping length to the reads dictionary in a list

                    # print("read_name, min_over, nm_tag")
                    # print(read_name, min_over, nm_tag)
                    # to_add = [int(min_over), int(nm_tag)]
                    to_add = [int(min_over)]
                    reads_dict[read_name].append(int(min_over))
                    lengths_dict[read_name] = int(read_length)


                    nm_tags_dict[read_name] = int(nm_tag)
                    df_mapping_all.loc[df_length] = row_to_add

    return reads_dict, lengths_dict, df_mapping_all, nm_tags_dict


def min_over_32del(pileupcolumn, position_list, lengths_dict, df_mapping_all, nm_tags_dict):
    """
    Calculate the minimum overlapping length for each read in the bam file mapped against the fake reference
    :param position:
    :param position_list:
    :return:
    """
    reads_dict = {}

    for pileupread in pileupcolumn.pileups:

        # If the read is not a deletion
        if not pileupread.is_del and not pileupread.is_refskip:
            read_name = pileupread.alignment.query_name

            # Position of the starting base of the read on the genome
            reference_start = pileupread.alignment.reference_start + 1

            # Position of the ending base of the read on the genome
            reference_end = pileupread.alignment.reference_end

            # Position in the read overlapping the pileup position
            query_position = pileupread.query_position + 1

            # print("query position:", query_position)
            read_length = pileupread.alignment.query_length

            read_sequence = pileupread.alignment.query_sequence

            nm_tag = pileupread.alignment.get_tag("NM")


            # I calculate the left and right overlaps of the read
            left = query_position
            right = reference_end - position_list[0] + 1

            # I assign the minimum overlapping length and I add it with the relative read name to the dictionary
            min_over = min(left, right)

            read_length = pileupread.alignment.query_length

            # I filter for only the reads under the theshold that I set
            if read_length <= length_threshold:

                # I save all in a row
                row_to_add = [sample, read_name, reference_start, reference_end, read_sequence, read_length, min_over, nm_tag, "del"]

                df_length = len(df_mapping_all)
                
                # That I add to a dataframe so I can analyse it afterwards
                # df_mapping_all.loc[df_length] = row_to_add

                if min_over >= ol_threshold:
                    
                    # to_add = [int(min_over), int(nm_tag)]
                    to_add = [int(min_over)]

                    reads_dict[read_name] = int(min_over)
                    lengths_dict[read_name] = int(read_length)

                    nm_tags_dict[read_name] = int(nm_tag)


                    df_mapping_all.loc[df_length] = row_to_add


    return reads_dict, lengths_dict, df_mapping_all, nm_tags_dict


def average_minimum_overlap(reads_dict):
    """
    Function to take the average of the minimum overlapping lengths over the 4 different coordinates couples of the
    deletion.

    :param reads_dict: dictionary containing read names and the relative minimum overlapping lengths for each read
    for each coordinate couple; dict
    :return: dictionary containing read name - average minimum overlap on the genome as key,values; dict
    """
    # I create an empty dictionary
    reads_dict_minimum = {}



    if reads_dict != {}:

        for key, value in reads_dict.items():

            # print("key, value")
            # print(key, value)
            # print(type(key), type(value))
            # If the read overlaps across both breakpoints for at least one coordinate couple, set min_over to 32
            if 32 in value:
                min_over = 32
                # nm_tag = value[1]

            # Else, min_over will be calculate as the average of the minimum overlaps for each coordinate couple
            else:
                min_over = mean(value)
                # nm_tag = value[1]

            reads_dict_minimum[key] = min_over

    return reads_dict_minimum


def pD_RR_b_():
    """
    Joint likelihoods from both the bam files, vs GRCH37 and vs Fake
    :return: pD_RR_b, pD_RD_b, pD_DD_b
    """
    pD_RR_b = pD_RR_g * pD_RR_d
    pD_RD_b = pD_RD_g * pD_RD_d
    pD_DD_b = pD_DD_g * pD_DD_d

    return pD_RR_b, pD_RD_b, pD_DD_b

def dict_to_list(reads_dict):
    """
    Convert the reads dictionary to list containing only the average minimum overlapping lengths without read names
    :param reads_dict:
    :return: reads_list: list containing the minimum overlapping lengths of the reads from the dictionary
    """
    reads_list = []

    # I convert the reads_dict from a normal dict to an Ordered Dict
    reads_dict = OrderedDict(reads_dict)
    if reads_dict != {}:
        for key, value in reads_dict.items():
            reads_list.append(value)

    return reads_list


# 32bp sequence posterior probabilities genotype calculation
# pD_G_2 = p(D|G) likelihood
# pD_2 = p(D)
# pG_D_2 = p(G|D)
# p(G|D) = p(G) * p(D|G) / p(D)

# p(D|G)
def p_D_G_2(reads_dict, which_bam):
    """
    Calculate p(D|G), i.e. the probability of the data (the reads) given the Genotype, i.e. that the sample has RR, RD,
    or DD genotype
    :param reads_dict:
    :param which_bam: to specify which bam I want to analyze; string
    :return:
    """
    # p_ref_r = p(ref|r), i.e. the probability of having the reference sequence given the observed read
    # p_del_r = p(del|r), i.e. the probability of having the deleted sequence given the observed read


    pD_RR_list, pD_RD_list, pD_DD_list = [], [], []

    # If the dict is not empty, i.e. if there is at least one read overlapping the region
    if reads_dict != {}:

        # Bam vs GRCh37 --> finding reference sequence
        if which_bam == "GRCh37":

            for key, value in reads_dict.items():
                # Step 1: calculate prob of reference or deleted sequence given the observed read
                p_ref_r = 1 - (1 / value) ** 2
                p_del_r = (1 - p_ref_r) / 2


                # Step 2: calculate prob of Data, of reads, given each possible genotype, so Ref/Ref, Ref/Del, Del/Del
                # Step 3: multiply (in this case sum, since they are logarithms), over all reads to calculate the
                # Likelihood of each genotype p(D|G)
                # this would be like doing p_ref_r/2 + p_ref_r/2, and it will just give p_ref_r so I simplified
                pD_RR_list.append(p_ref_r)
                pD_RD_list.append((p_ref_r + p_del_r) / 2)
                pD_DD_list.append(p_del_r)


            pD_RR = np.prod(pD_RR_list)
            pD_RD = np.prod(pD_RD_list)
            pD_DD = np.prod(pD_DD_list)

        # Bam vs fake 32del --> finding deletion
        if which_bam == "del":

            for key, value in reads_dict.items():
                # Step 1:calculate prob of reference or deleted sequence given the observed read
                p_del_r = 1 - (1 / value) ** 2
                p_ref_r = (1 - p_del_r) / 2


                # Step 2: calculate prob of Data, of reads, given each possible genotype, so Ref/Ref, Ref/Del, Del/Del
                # Step 3: multiply (in this case sum, since they are logarithms), over all reads to calculate the
                # Likelihood of each genotype p(D|G)
                pD_RR_list.append(p_ref_r)
                pD_RD_list.append((p_ref_r + p_del_r) / 2)
                pD_DD_list.append(p_del_r)

            pD_RR = np.prod(pD_RR_list)
            pD_RD = np.prod(pD_RD_list)
            pD_DD = np.prod(pD_DD_list)

    # If the dict is empty, there are no reads mapping to the region
    else:

        # Set the probabilities to be random, as 0.33
        pD_RR = 0.33
        pD_RD = 0.33
        pD_DD = 0.33

    return pD_RR, pD_RD, pD_DD
    # return math.log10(pD_RR), math.log10(pD_RD), math.log10(pD_DD)


# p(D)
def pD_2_(pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm, pD_RR_b, pD_RD_b, pD_DD_b):
    """
    p(D) = SUM p(Gi) p(D|Gi), for each genotype RR, RD, DD. As a prior I'll use the posterior probability for the
    rs333 haplotype that I calculated earlier"""

    # Calculate with normalized
    pD_2_norm = (pRR_D_joint_norm * pD_RR_b) + (pRA_D_joint_norm * pD_RD_b) + (pAA_D_joint_norm * pD_DD_b)

    # Calculate it when it's random
    pD_2_r = (0.33 * pD_RR_b) + (0.33 * pD_RD_b) + (0.33 * pD_DD_b)

    return pD_2_norm, pD_2_r


# p(G|D)
def pG_D_2(pG, pD_G, pD):

    pG_D_2 = (pG * pD_G) / pD
    # print("pG_D_2 = (pG * pD_G) / pD")
    # print("{} = ({} * {}) / {}".format(str(pG_D_2), str(pG), str(pD_G), str(pD)))

    return pG_D_2

# def pG_D_2_random(pG, pD_G, pD_r):
#     pG_D_2_r = (0.33 * pD_G) / pD_r

#     # print("pG_D_2_r = (pG * pD_G) / pD")
#     # print("{} = ({} * {}) / {}".format(str(pG_D_2_r), str(0.33), str(pD_G), str(pD_r)))
#     return pG_D_2_r


def tag_filtering(bamfile, reads_dict, lengths_dict):
    aligned_list = []

    # I create a list containing all the reads aligning in this region.
    # Each element of the list is a pysam.AlignmentSegment object
    for read in bamfile.fetch("3", 46414943, 46414980):
        aligned_list.append(read)
        
        
    # For each read name in reads_dict:
    for key in list(reads_dict):

        # For each read in the above list
        for read in aligned_list:

            # If the read from reads_dict is in the list AND does not have a perfect match, I'll remove it from the dictionary
            if key == read.query_name and read.get_tag("NM") != 0:
                del reads_dict[key]
                del lengths_dict[key]

    return reads_dict

def tag_filtering_ccr5kirstine(bamfile, reads_dict, lengths_dict):
    aligned_list = []

    # I create a list containing all the reads aligning in this region.
    # Each element of the list is a pysam.AlignmentSegment object
    for read in bamfile.fetch("CCR5_del32_120b.fasta", 30, 90):
        aligned_list.append(read)
        
        
    # For each read name in reads_dict:
    for key in list(reads_dict):

        # For each read in the above list
        for read in aligned_list:

            # If the read from reads_dict is in the list AND does not have a perfect match, I'll remove it from the dictionary
            if key == read.query_name and read.get_tag("NM") != 0:
                del reads_dict[key]
                del lengths_dict[key]

    return reads_dict

def write_probdf(prob_df, sample):
    """Function to write to file the probability dataframe of the 4 TOP SNPs along with their coverages etc"""    
    prob_df.to_csv(outdir + sample + "top4SNPs_prob_df.tsv", sep="\t")

# I write a file containing the settings that I used to run the script
def write_settings():
    settings_dict = collections.OrderedDict()
    settings_dict["samples"] = str(args.samples)
    settings_dict["files-extension"] = str(args.files_extension)
    settings_dict["folder-ref"] = str(args.folder_ref)
    settings_dict["folder-fake"] = str(args.folder_fake)
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


def write_results(coverage_ref, coverage_alt, dict_snps_cov, sample, N_reads_mapping_both):

    # results_filename = re.search("[\w\d]+$", args.output_folder).group() + ".tsv"
    # results_filepath = output_folder + results_filename
    with open(results_filepath, "a") as output:
        output.write(str(sample) + "\t" +  # Sample name
                     # Posterior probabilities of the haplotype
                     str(pRR_D_2_norm) + "\t" +  # p(RR|Data) normalized
                     str(pRD_D_2_norm) + "\t" +  # p(RD|Data) normalized
                     str(pDD_D_2_norm) + "\t" +   # p(DD|Data) normalized
                     # N of reads mapping to either S or E in the bamVSGrch37 AND mapping to Pos breakpoint in bamVSFake
                     str(len(reads_dict_ref)) + "\t" +
                     str(len(reads_dict_del)) + "\t" +
                     # Minimum overlapping lengths of each read mapping to the above mentioned bams
                     str(reads_list_ref) + "\t" +
                     str(reads_list_del) + "\t" +
                     # Lengths of the reads
                     str(lengths_list_ref) + "\t" +
                     str(lengths_list_del) + "\t" +
                     # Number of reads overlapping all the 82 SNPs, i.e. total coverage, not averaged
                     str(coverage_ref) + "\t" +
                     str(coverage_alt) + "\t" +
                     # Number of reads overlapping each of the 4 TOP SNPs, i.e. coverage at that base
                     str(dict_snps_cov["rs113341849"]) + "\t" + 
                     str(dict_snps_cov["rs113010081"]) + "\t" +
                     str(dict_snps_cov["rs11574435"]) + "\t" +  # rs11574435
                     str(dict_snps_cov["rs79815064"]) + "\t" +  # rs11574435
                    # Prior probabilities of the haplotype
                     str(pRR_D_joint_norm) + "\t" +  # p(RR)
                     str(pRA_D_joint_norm) + "\t" +  # p(RA)
                     str(pAA_D_joint_norm) + "\t" +  # p(AA)
                     # Reference or Deletion genotype Likelihoods from both bam files
                     str(pD_RR_b) + "\t" +  # p(Data|RR)
                     str(pD_RD_b) + "\t" +  # p(Data|RD)
                     str(pD_DD_b) + "\t" +  # p(Data|DD)
                     # p(D)
                     str(pD_2_norm) + "\t" + # p(D)
                     str(pRR_D_2_r) + "\t" +  # p(RR|Data) with pG = 0.33
                     str(pRD_D_2_r) + "\t" +  # p(RD|Data) with pG = 0.33
                     str(pDD_D_2_r) + "\t" +  # p(DD|Data) with pG = 0.33
                     str(N_reads_mapping_both) + "\n")  


def some_function(haplotype_list, haplo_df, ref_haplo_count_df):

    # I initialize dictionary where I'll store the reference and alternate bases called for each SNP
    dict_snps = collections.OrderedDict()

    dict_ref_haplo_count = collections.OrderedDict()
    referencecount, purereferencecount, haplocount, purehaplocount, notavail = 0, 0, 0, 0, 0
    # Iterate through each SNP 
    for snp in range(len(haplotype_list)):

        # Extract the different fields of the file
        id = haplotype_list[snp][0]
        chrom = str(3)
        coordinate = int(haplotype_list[snp][1])
        ref = haplotype_list[snp][2]
        alt = haplotype_list[snp][3]
        rsquared = haplotype_list[snp][4]

        # 1 - Extract all the reads bases mapping to each snp
        reads_list, other_list, ref_list, alt_list = extr_rbases_bam(bamvsref, chrom, coordinate, ref, alt, baq_snp, adjustment_threshold)

        # 2 - Calculate the number of reference and alternate bases called for each SNP
        ref_bases_n = len(ref_list)
        alt_bases_n = len(alt_list)

        if ref_bases_n > 0:
            referencecount +=1

        if ref_bases_n > 0 and alt_bases_n == 0:
                purereferencecount +=1

        if alt_bases_n > 0:
                haplocount +=1
        
        if alt_bases_n > 0 and ref_bases_n == 0:
                purehaplocount +=1



        # 3 - Save these in a dictionary 
        dict_snps["Sample"] = sample

        if ref_bases_n == 0 and alt_bases_n == 0:
            notavail +=1


            dict_snps[id + "_ref"] = [math.nan]
            dict_snps[id + "_alt"] = [math.nan]

        else:

            dict_snps[id + "_ref"] = [ref_bases_n]
            dict_snps[id + "_alt"] = [alt_bases_n]

    # print("referencecount:", referencecount)
    # print("purereferencecount:", purereferencecount)
    # print("haplocount:", haplocount)
    # print("purehaplocount:", purehaplocount)
    # print("notavail:", notavail)

    dict_ref_haplo_count["Sample"] = [sample]
    dict_ref_haplo_count["referencecount"] = [referencecount]
    dict_ref_haplo_count["purereferencecount"] = [purereferencecount]
    dict_ref_haplo_count["haplocount"] = [haplocount]
    dict_ref_haplo_count["purehaplocount"] = [purehaplocount]
    dict_ref_haplo_count["notavail"] = [notavail]

    row_sample = pd.DataFrame.from_dict(dict_snps)
    haplo_df = haplo_df.append(row_sample)

    row_ref_haplo_row_sample = pd.DataFrame.from_dict(dict_ref_haplo_count)
    ref_haplo_count_df = ref_haplo_count_df.append(row_ref_haplo_row_sample)
    return haplo_df, ref_haplo_count_df

############## Execution #################

# I write the header of the output file
results_filename = re.search("[\w\d]+$", args.output_folder).group() + ".tsv"

results_filepath = args.output_folder + "/" + results_filename

outdir = args.output_folder + "/prob_dfs/"

if not os.path.exists(outdir):
    os.makedirs(outdir)

with open(results_filepath, 'w') as output_header: output_header.write("Sample\tpRR_Data_n\tpRD_Data_n\tpDD_Data_n\tN_reads_ref\tN_reads_del\tMin_over_ref\tMin_over_del\tLengths_ref\tLengths_del\tCoverage_ref\tCoverage_alt\tSNP_1_rs113341849\tSNP_2_rs113010081\tSNP_3_rs11574435\tSNP_4_rs79815064\tp(RR)\tp(RA)\tp(AA)\tpData_RR\tpData_RD\tpData_DD\tpD_norm\tpRR_Data_r\tpRD_Data_r\tpDD_Data_r\tN_reads_mapping_both\n")


# Open the samples 
samples_list = samples_to_list()


# Initialize empty dataframe
df_mapping_all = pd.DataFrame(dtype=float, columns = ["sample", "read_name", "reference_start", "reference_end", "read_sequence", "read_length", "min_over", "n_mismatches", "alignment"])

df_ref_haplo_counts = pd.DataFrame(dtype=float, columns = ["sample", "Referencecount", "Purereferencecount", "haplocount", "Purehaplocount", "notAvail"])

# In these dataframes I will store, for each sample, which reads were assigned to ref and which to del, according to the current criteria (overlapping lenght)
df_reads_ref = pd.DataFrame(dtype=float, columns = ["sample", "read_name", "class"])

df_reads_del = pd.DataFrame(dtype=float, columns = ["sample", "read_name", "class"])

# If --haplotype option is activated --> write a table to file containing the reporting of all the 86 SNPs
if args.haplotype:

    haplotype_list = snp_haplo_list(args.haplotype)

    # e.g. [['rs58697594', '46275570', 'G', 'A', '0.8602'], ['rs73833032', '46276490', 'T', 'C', '0.8602']]
    # Need to convert this list of lists in another list of just the column names with the format
    # rs58697594_ref, rs58697594_alt, rs73833032_ref, rs73833032_alt
    
    rs_list_ref = [rs_id[0] + "_ref" for rs_id in haplotype_list]
    rs_list_alt = [rs_id[0] + "_alt" for rs_id in haplotype_list]

    # I contruct the rs_list as a list containing all the 82 SNPs _alt and _ref

    rs_list = []
    rs_list.extend(rs_list_ref)
    rs_list.extend(rs_list_alt)

    rs_list.sort()

    # print(rs_list)
    # I initialize an empty dataframe where I'll store the individual 86 SNPs reporting
    haplotype_df = pd.DataFrame(dtype=float)
    ref_haplo_count_df =pd.DataFrame(dtype=float)

# For each sample to analyse
for sample in samples_list:
    print(sample)

    # I parse the arguments given when executing the script
    bamvsref, bamvsdel, fasta_ref, fasta_fake, snp_file, haplotype_file, baq_snp, baq_deletion, length_threshold, ol_threshold, perfect_match, adjustment_threshold = open_files_args(sample)


    ############## Part 0: If --haplotype option is activated --> write a table to file containing the reporting of all the 86 SNPs
    if args.haplotype:

        haplotype_list = snp_haplo_list(haplotype_file)
        # I report all the SNPs called of the haplotype
        haplotype_df, ref_haplo_count_df = some_function(haplotype_list, haplotype_df, ref_haplo_count_df)


    ############## Part A: Prior Probability calculated as joint probability of top 4 SNPs' posterior probabilities - Execution #################

    # 1 - Extract SNPs from file
    # List of lists containing the 4 SNPs of the CEU rs333 haplotype with coordinates and R squared value to the rs333
    # e.g. [['rs58697594', '46275570', 'G', 'A', '0.8602'], ['rs73833032', '46276490', 'T', 'C', '0.8602']]
    snp_list = snp_haplo_list(snp_file)

    
    # 2 - Calculate Posterior probability of each SNP given each possible Genotype
    prob_df, coverage_ref, coverage_alt, coverage_other, dict_snps_cov = calc_snps_posteriors(snp_list)



    # 2.5 - I calculate the average coverage haplotype
    coverage_ref = coverage_ref / len(snp_list)
    coverage_alt = coverage_alt / len(snp_list)
    coverage_other = coverage_other / len(snp_list)

    # 3 - Weighted probs calculation: multiply each SNP's probability by the relative Rsquared and store in new columns.
    # 3 - Plus, add a column with their logarithm10 conversion
    prob_df = prob_to_weighted(prob_df)

    # I write the prob_df as a tab separated file
    write_probdf(prob_df, sample)

    # 4 - Joint probs calculation + normalization.
    pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm = calc_prob_joint(prob_df)


    ### IF I WANT TO PRINT THE DATAFRAME
    # Uncomment to show the entire dataframe
    # pd.set_option("display.max_columns", None)
    # pd.set_option("display.max_rows", None)
    # print(prob_df)

    ############## Part B: 32bp sequence posterior probabilities - Execution #################

    # 1 - I declare the lists containing the positions I want to check for overlapping reads
    position_list_reference = [[46414944, 46414975], [46414945, 46414976], [46414946, 46414977], [46414947, 46414978]]

    position_list_deletion = [46414943]

    # 2 - Calculation of the minimum overlapping lengths of the reads
    # In the dataframe df_mapping_all I put all the reads mapping, so both those that map vs reference and those that map vs collapsed
    reads_dict_ref, lengths_dict_ref, df_mapping_all, nm_tags_dict_ref = minimum_overlap(bamvsref, "3", position_list_reference, adjustment_threshold, df_mapping_all, baq = baq_deletion)
    reads_dict_del, lengths_dict_del, df_mapping_all, nm_tags_dict_del = minimum_overlap(bamvsdel, "3", position_list_deletion, adjustment_threshold, df_mapping_all, baq = baq_deletion)

    # 3 - Average of the overlapping lengths of all the 4 coordinates couples in the bam vs GRCh37
    reads_dict_ref = average_minimum_overlap(reads_dict_ref)




    # In case there are reads that overlap both the reference and the collapsed genome, I'll keep only the one that
    # has the lowest number of mismatches and the highest overlapping length

    N_reads_mapping_both = 0

    list_reads_mapping_both  = []

    for key in list(reads_dict_del.keys()):

        if key in reads_dict_ref:
            N_reads_mapping_both +=1
            # print("##########################")
            # print("nm_tags_dict_del[key], nm_tags_dict_ref[key]")
            # print(nm_tags_dict_del[key], nm_tags_dict_ref[key])
            # If the read vs del has a lower number of mismatches than the read vs ref
            if nm_tags_dict_del[key] < nm_tags_dict_ref[key]:

                # Assign the read to del, i.e. remove the read vs ref from its dictionary
                del reads_dict_ref[key]
                del lengths_dict_ref[key]


             # If the read vs del has a higher number of mismatches than the read vs ref
            elif nm_tags_dict_del[key] > nm_tags_dict_ref[key]:

                # Assign the read to ref, i.e. remove the read vs del from its dictionary
                del reads_dict_del[key]
                del lengths_dict_del[key]

            # If the number of mismatches is the same among the two:
            else:
 
                # If the overlapping length of the read vs del is lower than the one of the read vs ref:
                if reads_dict_del[key] < reads_dict_ref[key]:

                    # Assign the read to ref, i.e. remove the read vs del from its dictionary
                    del reads_dict_del[key]
                    del lengths_dict_del[key]

                elif reads_dict_ref[key] < reads_dict_del[key]:

                    del reads_dict_ref[key]
                    del lengths_dict_ref[key]

                # If also this is the same, remove from both of them
                else:

                    del reads_dict_ref[key]
                    del lengths_dict_ref[key]

                    del reads_dict_del[key]
                    del lengths_dict_del[key]




            
    # 4 - Filter by the XM:i:0 tag --> keep only the perfect matching reads

    if perfect_match == "yes":
        reads_dict_ref = tag_filtering(bamvsref, reads_dict_ref, lengths_dict_ref)
        reads_dict_del = tag_filtering(bamvsdel, reads_dict_del, lengths_dict_del)
        
    
    # 5 - Convert the dicts to lists, so it's easier to write in the output file
    reads_list_ref = dict_to_list(reads_dict_ref)
    reads_list_del = dict_to_list(reads_dict_del)

    
    lengths_list_ref = dict_to_list(lengths_dict_ref)
    lengths_list_del = dict_to_list(lengths_dict_del)

    # 6 - I calculate p(D|G) for both the bam vs GRCh37 and vs 32del
    pD_RR_g, pD_RD_g, pD_DD_g = p_D_G_2(reads_dict_ref, "GRCh37")
    pD_RR_d, pD_RD_d, pD_DD_d = p_D_G_2(reads_dict_del, "del")


    # 7 - I calculate the JOINT p(D|G) from the 2 bams
    pD_RR_b, pD_RD_b, pD_DD_b = pD_RR_b_()


    # 8 - p(D) calculation
    pD_2_norm, pD_2_r = pD_2_(pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm, pD_RR_b, pD_RD_b, pD_DD_b)

    # 9 - Posterior Probabilities p(G|D) for each "sequence genotype" using the normalized likelihoods

    pRR_D_2_norm = pG_D_2(pRR_D_joint_norm, pD_RR_b, pD_2_norm)
    pRD_D_2_norm = pG_D_2(pRA_D_joint_norm, pD_RD_b, pD_2_norm)
    pDD_D_2_norm = pG_D_2(pAA_D_joint_norm, pD_DD_b, pD_2_norm)

    # 10 - Posterior Probabilities p(G|D) for each "sequence genotype" considering the RANDOM haplotype
    pRR_D_2_r = pG_D_2(0.33, pD_RR_b, pD_2_r)
    pRD_D_2_r = pG_D_2(0.33, pD_RD_b, pD_2_r)
    pDD_D_2_r = pG_D_2(0.33, pD_DD_b, pD_2_r)

    # 11 - I append the results to the output file
    write_results(coverage_ref, coverage_alt, dict_snps_cov, sample, N_reads_mapping_both)

    # TODO: fix the write_settings
    # write_settings()

    # I need to write a new dataframe in which I put the reads_dict_ref and reads_dict_del, so these are where the model put the reads

    # I save the reads mapped to ref and those to del in 2 dataframes that I save

    for key in reads_dict_ref.keys():
        row_ref = [sample, key, "ref"]
        df_reads_ref.loc[len(df_reads_ref)] = row_ref

    for key in reads_dict_del.keys():
        row_del = [sample, key, "del"]
        df_reads_del.loc[len(df_reads_del)] = row_del
        


# print("haplotype_df")
# print(haplotype_df)
# print(haplotype_df.ndim)
# print(haplotype_file)
# print(haplotype_list)
# print(len(haplotype_list))
# I average the overlapping lengths 
df_mapping_all = (df_mapping_all
.assign(average_min_over = 
    lambda x: 
        x.groupby(["sample","read_name","alignment"])
        ["min_over"].transform("mean")))

df_mapping_all = df_mapping_all.drop(columns = ["min_over"])

df_mapping_all = df_mapping_all.drop_duplicates()



# I need to average the overlapping lengths of the ref
df_mapping_all.to_csv(args.output_folder + "/all_reads_mapping.tsv", sep="\t", quoting=csv.QUOTE_NONE, index = False)


df_reads_ref.to_csv(args.output_folder + "/reads_assigned_ref.tsv", sep="\t", quoting=csv.QUOTE_NONE, index = False)
df_reads_del.to_csv(args.output_folder + "/reads_assigned_del.tsv", sep="\t", quoting=csv.QUOTE_NONE, index = False)


if args.haplotype:
    haplotype_df.to_csv(args.output_folder + "/SNPS_reporting.tsv", sep = "\t", quoting=csv.QUOTE_NONE, index = False)
    ref_haplo_count_df.to_csv(args.output_folder + "/ref_haplo_counts.tsv", sep = "\t", quoting=csv.QUOTE_NONE, index = False)
end = time()
length = end - start
print("Time:", length)
