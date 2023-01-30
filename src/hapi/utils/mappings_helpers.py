from collections import defaultdict
from collections import OrderedDict
from statistics import mean
from typing import Union, Any, Tuple, List

import pandas as pd
import pysam

from hapi.utils.probabilities_helpers import (
    pD_G_,
    pD_,
    pG_D_
)


############## Part A FUNCTIONS DECLARATION ##############

def get_pilecolumns(bam_file: pysam.libcalignmentfile.AlignmentFile, baq: bool,
                    chrom: str, min_base_quality: int,
                    adjustment_threshold: int, min_mapping_quality: int,
                    fastafile: pysam.libcfaidx.FastaFile,
                    start: int,
                    end: int) -> pysam.libcalignmentfile.IteratorColumnRegion:
    """
    Makes an object with the extracted reads mapping to the genomic coordinate
    :param bam_file: bam file 
    :param baq: whether to perform BAQ (base alignment quality) calculation
    :param chrom: chromosome number of the position to extract
    :param min_base_quality: base quality
    :param adjustment_threshold: adjust mapping quality.
    :param min_mapping_quality: mapping quality
    :param fastafile: Fasta File of the genome with deletion
    :param start: start position were pileup is performed
    :param end: start position were pileup is performed
    :raises ValueError: Raises error if baq is not True or False
    :return: Iterator object pysam pileup object: an iterable which represents all the reads in the SAM file that map to a single base in the reference sequence.
    The list of reads are represented as PileupRead objects in the PileupColumn.pileups property
    """
    # Each iteration returns a PileupColumn which represents all the reads in the SAM file that map to a single base in 
    # the reference sequence. 
    if baq == False:
        pileupcolumns = bam_file.pileup(chrom, start, end, truncate=True,
                                        min_base_quality=min_base_quality,
                                        adjust_capq_threshold=adjustment_threshold,
                                        min_mapping_quality=min_mapping_quality)
    elif baq == True:
        pileupcolumns = bam_file.pileup(chrom, start, end, truncate=True,
                                        stepper="samtools",
                                        fastafile=fastafile, compute_baq=True,
                                        min_base_quality=min_base_quality,
                                        adjust_capq_threshold=adjustment_threshold,
                                        min_mapping_quality=min_mapping_quality)
    else:
        raise ValueError('baq parameter selected neither True nor False')

    return pileupcolumns


def extr_rbases_bam(bamvsref_file: pysam.libcalignmentfile.AlignmentFile,
                    chrom: str, coordinate: int, ref: str, alt: str, baq: bool,
                    fasta_ref: pysam.libcfaidx.FastaFile,
                    adjustment_threshold: int, length_threshold: int,
                    min_base_quality: int = 0,
                    min_mapping_quality: int = 0) -> Tuple[
    list, list, list, list]:
    """ 
    Extract the read bases, i.e. the bases of the reads that map to a specific
    position in the bam file. bq and mq are base quality and mapping quality.
    Change them if I want to filter.

    :param bamvsref_file: bam file 
    :param chrom: chromosome number of the position to extract
    :param coordinate: coordinate of the position to extract
    :param ref: reference allele at the position
    :param alt: alternate allele at the position
    :param baq: whether to perform BAQ (base alignment quality) calculation
    :param fasta_ref: Fasta File of the genome with deletion
    :param adjustment_threshold: adjust mapping quality.
    :param length_threshold: value to keep only reads with read length < length_threshold. Not used in the final script, i.e. it is set to 1000 so no filtering is performed
    :param min_base_quality: mapping quality
    :param min_mapping_quality: base alignment quality filtering
    :return: reads_list, read bases that at the position correspond to
        either the Reference or Alternate allele 
    :return: other_list, read bases that at the position do not correspond to
        neither the Reference nor Alternate allele 
    :return: ref_list, read bases that at the position do not correspond to
        the Reference allele 
    :return: alt_list: read bases that at the position do not correspond to
        Alternate allele
    """
    # reads_list contains the reads overlapping to the position when the
    # base found corresponds to either the Ref or the Alt If the base found
    # is different from either the Ref or the Alt, I'll put the read in:pysam.libcalignmentfile.AlignmentFile
    # other_list

    reads_list, other_list, ref_list, alt_list = [], [], [], []

    pileupcolumns = get_pilecolumns(bamvsref_file, baq, chrom,
                                    min_base_quality, adjustment_threshold,
                                    min_mapping_quality, fastafile=fasta_ref,
                                    start=coordinate - 1, end=coordinate)

    for pileupcolumn in pileupcolumns:
        reads_list, other_list, ref_list, alt_list = extract_lists(
            pileupcolumn, ref, alt, length_threshold, reads_list, other_list,
            ref_list, alt_list)
    return reads_list, other_list, ref_list, alt_list


def extract_lists(pileupcolumn: pysam.libcalignedsegment.PileupColumn,
                  ref: str, alt: str, length_threshold: int, reads_list: list,
                  other_list: list, ref_list: list, alt_list: list) -> Tuple[
    list, list, list, list]:
    """
    Function to extract the reads mapping to a certain position from the
    PileupColumn object. If the read base at the position corresponds to
    either the Reference or Alternate allele of the SNP as written in the
    file, the read will be put in the object "reads_list". If not, the read
    will be put in the "other_list"

    :param pileupcolumn: pysam PileupColumn object which represents all the reads in the SAM file that map to a single base from the reference sequence. 
    :param ref: reference allele at the position
    :param alt: alternate allele at the position
    :param length_threshold: value to keep only reads with read length < length_threshold. Not used in the final script, i.e. it is set to 1000 so no filtering is performed
    :param reads_list: read bases that at the position correspond to
        either the Reference or Alternate allele 
    :param other_list: read bases that at the position do not correspond to
        neither the Reference nor Alternate allele 
    :param ref_list: read bases that at the position do not correspond to
        the Reference allele 
    :param alt_list: read bases that at the position do not correspond to
        Alternate allele 
    :return: reads_list, read bases that at the position correspond to
        either the Reference or Alternate allele 
    :return: other_list, read bases that at the position do not correspond to
        neither the Reference nor Alternate allele 
    :return: ref_list, read bases that at the position do not correspond to
        the Reference allele 
    :return: alt_list: read bases that at the position do not correspond to
        Alternate allele
    """

    # For each read mapping to the PileupColumn position
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:

            # Uncomment if want to print the coverage at the position
            # print('\tbase in read %s = %s' %
            #       (pileupread.alignment.query_name,
            #        pileupread.alignment.query_sequence[pileupread.query_position]))

            # I save the read name, the read base that map to the position,
            # and the base quality
            read_name = pileupread.alignment.query_name
            base = pileupread.alignment.query_sequence[
                pileupread.query_position]
            quality = pileupread.alignment.query_alignment_qualities[
                pileupread.query_position]
            read_length = pileupread.alignment.query_length

            read_info_list = [base, quality, read_name, read_length]

            if read_length <= length_threshold:
                # If the read base corresponds to the Reference or to the
                # Alternate allele, put it in the reads_list
                if base == ref or base == alt:
                    reads_list.append(read_info_list)

                if base == ref:
                    ref_list.append(read_info_list)

                if base == alt:
                    alt_list.append(read_info_list)

                # If the read base is not like the reference nor like the
                # alternate allele, put it in the other_list
                if base != ref and base != alt:
                    other_list.append(read_info_list)

    return reads_list, other_list, ref_list, alt_list


def calc_snps_posteriors(snp_list: list,
                         bamvsref: pysam.libcalignmentfile.AlignmentFile,
                         chrom: str, fasta_ref: pysam.libcfaidx.FastaFile,
                         baq_snp: bool,
                         adjustment_threshold: int, length_threshold: int,
                         top4_snps_list: list) -> Tuple[
    pd.DataFrame, float, float, float, dict]:
    """ 
    Calculate the Posterior Probabilities for each SNP, storing the
    information in a dataframe and saving the statistics relative to the
    total coverage on the SNPs and to the coverage of each of the 4 TOP SNPs

    :param snp_list: list of the SNPs
    :param bamvsref: bam file 
    :param chrom: chromosome number of the position to extract
    :param fasta_ref: Fasta File of the genome 
    :param baq_snp: whether to perform BAQ (base alignment quality) calculation on the SNPs
    :param adjustment_threshold: adjust mapping quality.
    :param length_threshold: value to keep only reads with read length < length_threshold. Not used in the final script, i.e. it is set to 1000 so no filtering is performed
    :param top4_snps_list: list of the 4 top SNPs

    :return: prob_df, dataframe containing the probabilities
    :return: coverage_ref and coverage_alt, counter containing the total number
        of bases mapping to all the SNPs, a total coverage
    :return: dict_snps_cov, dictionary containing the coverage for each of the
        4 TOP SNPs
    """

    # I initialize an empty dataframe where I'll store the Posterior
    # probabilities of each genotype at each position

    prob_list = []

    # I initialize an empty dictionary and an empty counter to save the
    # coverage statistics of the SNPs
    coverage_ref, coverage_alt, coverage_other = 0, 0, 0
    dict_snps_cov = {}

    # Iterate through each SNP 

    for snp in snp_list:

        idx, coordinate, ref, alt, rsquared = snp[0:5]
        coordinate = int(coordinate)

        # 1 - Extract all the reads bases mapping to each snp
        reads_list, other_list, ref_list, alt_list = extr_rbases_bam(bamvsref,
                                                                     chrom,
                                                                     coordinate,
                                                                     ref, alt,
                                                                     baq_snp,
                                                                     fasta_ref,
                                                                     adjustment_threshold,
                                                                     length_threshold)

        # 2 - Update SNPs coverage statistics
        dict_snps_cov = coverage_dict(dict_snps_cov, idx, alt_list,
                                      top4_snps_list)

        coverage_ref += len(ref_list)
        coverage_alt += len(alt_list)
        coverage_other += len(other_list)

        # 3 - Calculate p(D|G) for each genotype, so for Ref/Ref (pD_RR),
        # Ref/Del (pD_RA), Del/Del (pD_AA)
        pD_RR, pD_RA, pD_AA = pD_G_(reads_list, ref, alt)

        # 4 - Calculate p(D), i.e. denominator of the final equation. 
        pD = pD_(pD_RR, pD_RA, pD_AA)

        # 5 - Calculate p(G|D), the posterior probability of each deletion
        # genotype.
        pRR_D, pRA_D, pAA_D = pG_D_(pD_RR, pD_RA, pD_AA, pD)

        # If reads_list is empty it means that there were no reads
        # overlapping the position, so I'll put 0.33
        if not reads_list:
            snp_result = {
                "id": idx,
                "chrom": chrom,
                "ref": ref,
                "alt": alt,
                "rsquared": float(rsquared),
                "other": 0,
                "n_total_reads": 0,
                "n_ref_reads": 0,
                "n_alt_reads": 0,
                "P(RR|D)": 0.33,
                "P(RA|D)": 0.33,
                "P(AA|D)": 0.33
            }
        else:
            snp_result = {
                "id": idx,
                "chrom": chrom,
                "ref": ref,
                "alt": alt,
                "rsquared": float(rsquared),
                "other": 0,
                "n_total_reads": len(reads_list),
                "n_ref_reads": len(ref_list),
                "n_alt_reads": len(alt_list),
                "P(RR|D)": pRR_D,
                "P(RA|D)": pRA_D,
                "P(AA|D)": pAA_D
            }

        # If there are reads with a base other than the REF or the ALT,
        # add them as a list in the column "other"
        if not other_list:
            snp_result["other"] = other_list

        prob_list.append(snp_result)

    prob_df = pd.DataFrame.from_records(prob_list)

    # I calculate the average coverage haplotype
    coverage_ref = coverage_ref / len(snp_list)
    coverage_alt = coverage_alt / len(snp_list)
    coverage_other = coverage_other / len(snp_list)

    return prob_df, coverage_ref, coverage_alt, coverage_other, dict_snps_cov


def coverage_dict(dict_snps_cov: dict, snp_id: str, alt_list: list,
                  top4_snps_list: List[str]) -> dict:
    """
    Updates the dictionary counter "dict_snps_cov" with the number of reads
    overlapping the 4 TOP SNPs

    :param dict_snps_cov: dict with the number of reads overlapping the 
        4 TOP SNPs
    :param snp_id: SNP name
    :param alt_list: read bases that at the position do not correspond to
        Alternate allele
    :param top4_snps_list: list of the 4 top SNPs
    :return: dict with the number of reads overlapping the 
        4 TOP SNPs

    """
    for top_snp in top4_snps_list:
        if snp_id == top_snp:
            dict_snps_cov[top_snp] = len(alt_list)

    return dict_snps_cov


############## Part B FUNCTIONS DECLARATION ##############

def minimum_overlap(bam_file: pysam.libcalignmentfile.AlignmentFile,
                    chrom: str, position_list: list, adjustment_threshold: int,
                    mapping_all, length_threshold: int, ol_threshold: int,
                    sample: str,
                    fasta_coll: pysam.libcfaidx.FastaFile,
                    fasta_ref: pysam.libcfaidx.FastaFile, baq: bool,
                    overlap_type: str,
                    min_base_quality: int = 30,
                    min_mapping_quality: int = 30) -> Tuple[
    Union[dict, defaultdict], dict, dict, List[dict]]:
    """
    Function to calculate the minimum length that a read overlaps the: -
    Starting and Ending position of the 32deletion in the bam aligned
    against the reference genome GRCh37 --> To find reads NOT HAVING the
    deletion, i.e. having the reference 32bp sequence - Position breakpoint
    of the 32deletion in the bam aligned against the coll 32del reference
    --> To find reads HAVING the deletion, i.e. not having the reference
    32bp sequence

    :param bam_file: bam file
    :param chrom: chromosome number
    :param position_list: list containing the 4 different Starting and Ending coordinates of the deletion in the Reference or the unique Starting and Ending coordinates in the Collapsed Reference
    :param adjust_threshold: adjust mapping quality.
    :param mapping_all: a list of dicts with stored reseults
    :param length_threshold: value to keep only reads with read length < length_threshold. Not used in the final script, i.e. it is set to 1000 so no filtering is performed
    :param number of nucleotides with which each read needs to overlap the deletion or reference coordinates in order to be kept
    :sample: a SNP
    :fasta_coll: Fasta File of the genome with deletion
    :fasta_ref: Fasta File of the genome 
    :param baq: base alignment quality filtering; string
    :param overlap_type: is it reference or deletion
    :param min_base_quality: base quality; 
    :param min_mapping_quality: mapping quality

    :return: reads_dict: the dictionary containing reads with their minimum overlaps
    :return: lengths_dict: dictionary containing each read with its minimum length with which it overlaps the deletion or reference coordinates
    :return: mapping_all: a list of dicts with stored reseults
    :return: nm_tags_dict: dictionary containing each read with its Samtools NM tag, i.e. the count of mismatches between the read and the reference 
             Needed for the function remove_overlaps, with which, if there are reads that overlap both the reference and the collapsed genome, I'll keep only the one that has the lowest number of mismatches and the highest overlapping length
    """
    lengths_dict = {}
    nm_tags_dict = {}
    # If the list contains 4 elements --> Bam aligned vs Reference,
    # to detect reads NOT having the deletion
    if overlap_type == "ref":
        # I initialize an empty dictionary in which the values of the keys
        # will be of type list
        reads_dict = defaultdict(list)

        # For each coordinates couples
        for start_end in position_list:

            # For each of the S and E
            for pos in start_end:
                pileupcolumns = get_pilecolumns(bam_file, baq, chrom,
                                                min_base_quality,
                                                adjustment_threshold,
                                                min_mapping_quality,
                                                fastafile=fasta_ref,
                                                start=pos - 1, end=pos)

                for pileupcolumn in pileupcolumns:
                    reads_dict, lengths_dict, mapping_all, nm_tags_dict = min_over_reference_or_32del(
                        pileupcolumn, reads_dict, lengths_dict, mapping_all,
                        nm_tags_dict, length_threshold, ol_threshold, sample,
                        overlap_type, position_list=start_end, pos=pos)

    # If the list contains 1 element, i.e. P --> Bam aligned vs Collapsed Reference, 
    # to detect reads HAVING the deletion
    elif overlap_type == "del":

        reads_dict = {}

        pileupcolumns = get_pilecolumns(bam_file, baq, chrom, min_base_quality,
                                        adjustment_threshold,
                                        min_mapping_quality,
                                        fastafile=fasta_coll,
                                        start=position_list[0] - 1,
                                        end=position_list[0])

        for pileupcolumn in pileupcolumns:
            reads_dict, lengths_dict, mapping_all, nm_tags_dict = min_over_reference_or_32del(
                pileupcolumn, reads_dict, lengths_dict, mapping_all,
                nm_tags_dict, length_threshold, ol_threshold, sample,
                overlap_type, position_list=position_list, pos=None)
    else:
        raise ValueError('selected overlap_type is not del or ref')

    return reads_dict, lengths_dict, mapping_all, nm_tags_dict


def min_over_reference_or_32del(
        pileupcolumn: pysam.libcalignedsegment.PileupColumn,
        reads_dict: Union[dict, defaultdict], lengths_dict: dict,
        mapping_all: List[dict], nm_tags_dict: dict, length_threshold: int,
        ol_threshold: int, sample: str, overlap_type: str,
        position_list: list, pos: Union[int, Any]) -> Tuple[
    Union[dict, defaultdict], dict, List[dict], dict]:
    """
    Function to extract each read from the pileupcolumn and calculate the the minimum overlapping length with which it maps to the Reference or to the
    Collapsed Reference, number of mismatches, read length, read sequence etc to save everything in a file
    :param pileupcolumn: pysam PileupColumn object which represents all the reads in the SAM file that map to a single base from the reference sequence. 
    :param reads_dict: the dictionary containing reads with their minimum overlaps
    :param lengths_dict: dictionary containing each read with its minimum length with which it overlaps the deletion or reference coordinates
    :param mapping_all: a list of dicts with stored reseults
    :param nm_tags_dict: dictionary containing each read with its Samtools NM tag, i.e. the count of mismatches between the read and the reference. See param spec of function minimum_overlap for more details
    :param length_threshold: value to keep only reads with read length < length_threshold. Not used in the final script, i.e. it is set to 1000 so no filtering is performed
    :param number of nucleotides with which each read needs to overlap the deletion or reference coordinates in order to be kept
    :param sample: a SNP
    :param overlap_type: is it reference or deletion
    :param position_list: list containing the 4 different Starting and Ending coordinates of the deletion in the Reference or the unique Starting and Ending coordinates in the Collapsed Reference
    :param pos: used to iterate over the positions list
    :return: reads_dict: the dictionary containing the reads with their
        minimum overlaps
    :return: lengths_dict: dictionary containing each read with its minimum length with which it overlaps the deletion or reference coordinates
    :return: mapping_all: a list of dicts with stored reseults
    :return: nm_tags_dict: dictionary containing each read with its Samtools NM tag, i.e. the count of mismatches between the read and the reference. See param spec of function minimum_overlap for more details
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

            # Filter for only the reads under the threshold that I set (1000 bp, i.e. not filtering cause no read will be longer than 1000)
            if read_length <= length_threshold:

                if overlap_type == 'ref':
                    S, E = position_list
                    if reference_start <= S and reference_end >= E:
                        # Minimum overlapping length is 32
                        min_over = 32
                    else:
                        # I calculate the left and right overlap of the read
                        left = query_position
                        right = reference_end - pos + 1

                        # The minimum overlapping length is the minimum
                        # between these two
                        min_over = min(left, right)

                # I calculate the left and right overlaps of the read
                elif overlap_type == 'del':
                    left = query_position
                    right = reference_end - position_list[0] + 1

                    # I assign the minimum overlapping length and I add it
                    # with the relative read name to the dictionary
                    min_over = min(left, right)

                if min_over >= ol_threshold:
                    # I append the minimum overlapping length to the reads
                    # dictionary in a list

                    if overlap_type == 'ref':
                        reads_dict[read_name].append(int(min_over))
                    elif overlap_type == 'del':
                        reads_dict[read_name] = int(min_over)

                    lengths_dict[read_name] = int(read_length)
                    nm_tags_dict[read_name] = int(nm_tag)

                    row_to_add = {
                        "sample": sample,
                        "read_name": read_name,
                        "reference_start": reference_start,
                        "reference_end": reference_end,
                        "read_sequence": read_sequence,
                        "read_length": read_length,
                        "min_over": min_over,
                        "n_mismatches": nm_tag,
                        "alignment": overlap_type
                    }
                    mapping_all.append(row_to_add)

    return reads_dict, lengths_dict, mapping_all, nm_tags_dict


def average_minimum_overlap(reads_dict: Union[dict, defaultdict]) -> dict:
    """
    Function to take the average of the minimum overlapping lengths over the
    4 different coordinates couples of the deletion.

    :param reads_dict: dictionary containing read names and the relative
    minimum overlapping lengths for each read for each coordinate couple
    :return: dictionary containing read name - average minimum overlap
    on the genome as key,values
    """

    # If the read overlaps across both breakpoints for at least one
    # coordinate couple, set min_over to 32 Else, min_over will be calculate
    # as the average of the minimum overlaps for each coordinate couple
    return dict(
        (k, 32) if 32 in v else (k, mean(v)) for k, v in reads_dict.items())


def perfect_match_filtering(bamfile: pysam.libcalignmentfile.AlignmentFile,
                  reads_dict: Union[dict, defaultdict], lengths_dict: int,
                  chrom: str, start: int, end: int) -> Union[
    dict, defaultdict]:
    """
    Function to keep only the reads with the perfect matching 
    :param bamfile: bam file
    :param reads_dict: the dictionary containing the reads with their
        minimum overlaps
    :param lengths_dict: dictionary containing each read with its minimum length with which it overlaps the deletion or reference coordinates
    :param chrom: chromosome number of the position to extract
    :param start: starting coordinate of the deletion (46414943) in the Collapsed reference
    :param end: ending coordinate of the deletion (46414980) in the Collapsed reference
    :param reads_dict: the dictionary containing the reads with their
        minimum overlaps
    """

    # For each read name in reads_dict:
    for key in reads_dict.keys():

        # For each read aligning in this region.
        # Each element of the list is a pysam.AlignmentSegment object
        for read in bamfile.fetch(chrom, start, end):

            # If the read from reads_dict is in the list AND does not have a
            # perfect match, I'll remove it from the dictionary
            if key == read.query_name and read.get_tag("NM") != 0:
                del reads_dict[key]
                del lengths_dict[key]

    return reads_dict


def snps_reporting(haplotype_list: List[str],
                  bamvsref: pysam.libcalignmentfile.AlignmentFile, chrom: str,
                  baq_snp: bool,
                  adjustment_threshold: int, length_threshold: int,
                  sample: str, fasta_ref: pysam.libcfaidx.FastaFile) -> Tuple[
    List[dict], List[dict]]:
    """ -> tuple[list[dict], list[dict]]
    Function to save the info about each SNP from the haplotype list. Not needed for HAPI, function created
    just to have the equivalent of what Kirstine was doing before I started to work on the project.

    :param haplotype_list: list of the haplotype SNPs
    :param bamvsref: bam file
    :param chrom: chromosome number of the position to extract
    :param baq_snp: whether to perform BAQ (base alignment quality) calculation on the SNPs
    :param adjustment_threshold: adjust mapping quality.
    :param length_threshold: value to keep only reads with read length < length_threshold. Not used in the final script, i.e. it is set to 1000 so no filtering is performed
    :param sample: SNP
    :param fasta_ref: Fasta File of the genome with deletion
    :return: haplo_results_list: list with each sample and its reference and alternate bases called for each snp
    :return: ref_haplo_count_list: list with each sample and its referencecount, purereferenceount, haplocount, purehaplocount, notavail
    """

    # I initialize dictionary where I'll store the reference and alternate
    # bases called for each SNP
    haplo_results_list, ref_haplo_count_list = [], []
    dict_snps = OrderedDict()
    referencecount, purereferencecount, haplocount, purehaplocount, notavail = 0, 0, 0, 0, 0

    # Iterate through each SNP 
    for snp in haplotype_list:

        idx, coordinate, ref, alt, rsquared = snp[0:5]
        coordinate = int(coordinate)

        # 1 - Extract all the reads bases mapping to each snp
        reads_list, other_list, ref_list, alt_list = extr_rbases_bam(bamvsref,
                                                                     chrom,
                                                                     coordinate,
                                                                     ref, alt,
                                                                     baq_snp,
                                                                     fasta_ref,
                                                                     adjustment_threshold,
                                                                     length_threshold)

        # 2 - Calculate the number of reference and alternate bases called
        # for each SNP
        ref_bases_n = len(ref_list)
        alt_bases_n = len(alt_list)

        if ref_bases_n > 0:
            referencecount += 1

        if ref_bases_n > 0 and alt_bases_n == 0:
            purereferencecount += 1

        if alt_bases_n > 0:
            haplocount += 1

        if alt_bases_n > 0 and ref_bases_n == 0:
            purehaplocount += 1

        # 3 - Save these in a dictionary
        dict_snps["Sample"] = sample

        if ref_bases_n == 0 and alt_bases_n == 0:
            notavail += 1

            dict_snps[idx + "_ref"] = [math.nan]
            dict_snps[idx + "_alt"] = [math.nan]

        else:

            dict_snps[idx + "_ref"] = [ref_bases_n]
            dict_snps[idx + "_alt"] = [alt_bases_n]

    dict_ref_haplo_count = {
        "Sample": [sample],
        "referencecount": [referencecount],
        "purereferencecount": [purereferencecount],
        "haplocount": [haplocount],
        "purehaplocount": [purehaplocount],
        "notavail": [notavail]
    }

    haplo_results_list.append(dict_snps)
    ref_haplo_count_list.append(dict_ref_haplo_count)

    return haplo_results_list, ref_haplo_count_list


def remove_overlaps(reads_dict_del: dict, reads_dict_ref: dict,
                    nm_tags_dict_del: dict,
                    nm_tags_dict_ref: dict, lengths_dict_ref: dict,
                    lengths_dict_del: dict) -> Tuple[
    dict, dict, dict, dict, dict, dict, int]:
    """ 
    In case there are reads that overlap both the reference and the Collapsed genome, I'll keep only the one that has the lowest
    number of mismatches and the highest overlapping length
    :param reads_dict_del: dictionary containing each read mapping to the Collapsed genome
    :param reads_dict_ref: dictionary containing each read mapping to the Reference genome
    :param nm_tags_dict_del: dictionary containing each read with its Samtools NM tag, i.e. the count of mismatches between the read and the Collapsed reference. See param spec of function minimum_overlap for more details
    :param nm_tags_dict_ref: dictionary containing each read with its Samtools NM tag, i.e. the count of mismatches between the read and the reference. See param spec of function minimum_overlap for more details
    :param lengths_dict_ref: dictionary containing each read mapping to the Reference genome with its length
    :param lengths_dict_del: dictionary containing each read mapping to the Collapsed genome with its length
    """
    n_reads_mapping_both = 0

    for key in list(reads_dict_del.keys()):

        if key in reads_dict_ref:
            n_reads_mapping_both += 1

            # If the read vs del has a lower number of mismatches than the
            # read vs ref
            if nm_tags_dict_del[key] < nm_tags_dict_ref[key]:

                # Assign the read to del, i.e. remove the read vs ref from
                # its dictionary
                del reads_dict_ref[key]
                del lengths_dict_ref[key]


            # If the read vs del has a higher number of mismatches than the
            # read vs ref
            elif nm_tags_dict_del[key] > nm_tags_dict_ref[key]:

                # Assign the read to ref, i.e. remove the read vs del from
                # its dictionary
                del reads_dict_del[key]
                del lengths_dict_del[key]

            # If the number of mismatches is the same among the two:
            else:

                # If the overlapping length of the read vs del is lower than
                # the one of the read vs ref:
                if reads_dict_del[key] < reads_dict_ref[key]:

                    # Assign the read to ref, i.e. remove the read vs del
                    # from its dictionary
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

    return (reads_dict_del, reads_dict_ref, nm_tags_dict_del, nm_tags_dict_ref,
            lengths_dict_ref, lengths_dict_del, n_reads_mapping_both)
