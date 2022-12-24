import pandas as pd
from collections import defaultdict
from collections import OrderedDict
from statistics import mean
from hapi.utils.probabilities_helpers import *


############## Part A FUNCTIONS DECLARATION ##############

def get_pilecolumns(bam_file, baq, chrom, min_base_quality, adjustment_threshold, min_mapping_quality, fastafile, start, end):
    
    if baq == False:
        pileupcolumns = bam_file.pileup(chrom, start, end, truncate=True,
                                            min_base_quality=min_base_quality,
                                            adjust_capq_threshold=adjustment_threshold,
                                            min_mapping_quality=min_mapping_quality)
    elif baq == True:
        pileupcolumns = bam_file.pileup(chrom, start, end, truncate=True,
                                            stepper="samtools", fastafile=fastafile, compute_baq=True,
                                            min_base_quality=min_base_quality,
                                            adjust_capq_threshold=adjustment_threshold,
                                            min_mapping_quality=min_mapping_quality)
    else:
        raise ValueError('baq parameter selected neither True nor False')
    
    return pileupcolumns

# def extr_rbases_bam(bamvsref_file, chrom, coordinate, ref, alt, baq, adjustment_threshold, min_base_quality=30,
#                     min_mapping_quality=30):

# Trying to do the same as Kirstine
def extr_rbases_bam(bamvsref_file, chrom, coordinate, ref, alt, baq, fasta_ref, adjustment_threshold, length_threshold, min_base_quality=0, min_mapping_quality=0):
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

    # reads_list contains the reads overlapping to the position when the base found corresponds to either the Ref or the Alt
    # If the base found is different from either the Ref or the Alt, I'll put the read in other_list

    reads_list, other_list, ref_list, alt_list = [], [], [], []
    
    pileupcolumns = get_pilecolumns(bamvsref_file, baq, chrom, min_base_quality, adjustment_threshold, min_mapping_quality, fastafile=fasta_ref, start=coordinate-1, end=coordinate)
    
    for pileupcolumn in pileupcolumns:
        reads_list, other_list, ref_list, alt_list = extract_lists(pileupcolumn, ref, alt, length_threshold)
    
    return reads_list, other_list, ref_list, alt_list



def extract_lists(pileupcolumn, ref, alt, length_threshold):
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
            
            read_info_list = [base, quality, read_name, read_length]
            
            if read_length <= length_threshold:
            # if True:
                # If the read base corresponds to the Reference or to the Alternate allele, put it in the reads_list
                if base == ref or base == alt:
                    reads_list.append(read_info_list)

                if base == ref:
                    ref_list.append(read_info_list)

                if base == alt:
                    alt_list.append(read_info_list)

                # If the read base is not like the reference nor like the alternate allele, put it in the other_list
                if base != ref and base != alt:
                    other_list.append(read_info_list)

    return reads_list, other_list, ref_list, alt_list

def calc_snps_posteriors(snp_list, bamvsref, fasta_ref, baq_snp, adjustment_threshold, length_threshold):
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
    
    for snp in snp_list:

        chrom = str(3)
        idx, coordinate, ref, alt, rsquared = snp[0:5]
        coordinate = int(coordinate)

        # 1 - Extract all the reads bases mapping to each snp
        reads_list, other_list, ref_list, alt_list = extr_rbases_bam(bamvsref, chrom, coordinate, ref, alt, baq_snp, fasta_ref, adjustment_threshold, length_threshold)

        # 2 - Update SNPs coverage statistics
        dict_snps_cov = coverage_dict(dict_snps_cov, reads_list, idx, other_list, ref_list, alt_list)

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
            row = pd.DataFrame([[idx, chrom, ref, alt, float(rsquared), 0, 0, 0, 0, 0.33, 0.33, 0.33]],
                               columns=columns_df,
                               index=[coordinate])
        else:
            row = pd.DataFrame([[idx, chrom, ref, alt, float(rsquared), 0, len(reads_list), len(ref_list), len(alt_list), pRR_D, pRA_D, pAA_D]],
                               columns=columns_df, index=[coordinate])

        # Whatever is the case, I'll append the row with this SNP to the probability dataframe
        prob_df = prob_df.append(row)

        # If there are reads with a base other than the REF or the ALT, add them as a list in the column "other"
        if not other_list == []:
            prob_df = prob_df.astype({"other": object})

            prob_df.at[coordinate, "other"] = other_list


    return prob_df, coverage_ref, coverage_alt, coverage_other, dict_snps_cov



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

############## Part B FUNCTIONS DECLARATION ##############

def minimum_overlap(bam_file, chrom, position_list, adjustment_threshold, df_mapping_all, length_threshold, ol_threshold, sample, fasta_fake, fasta_ref, baq,  min_base_quality=30, min_mapping_quality=30):
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

            # For each of the S and E
            for pos in start_end:
                pileupcolumns = get_pilecolumns(bam_file, baq, chrom, min_base_quality, adjustment_threshold,
                                                min_mapping_quality, fastafile=fasta_ref, start=pos-1, end=pos)
              
                for pileupcolumn in pileupcolumns:
                    
                    reads_dict, lengths_dict, df_mapping_all, nm_tags_dict = min_over_reference_or_32del(pileupcolumn, reads_dict, lengths_dict, df_mapping_all, nm_tags_dict, length_threshold, ol_threshold, sample, position_list=start_end, pos=pos, min_over_type='ref')

        return reads_dict, lengths_dict, df_mapping_all, nm_tags_dict

    # If the list contains 1 element, i.e. P --> Bam aligned vs fake Reference, to detect reads HAVING the deletion
    elif len(position_list) == 1:

        reads_dict = {}
        
        pileupcolumns = get_pilecolumns(bam_file, baq, chrom, min_base_quality, adjustment_threshold, min_mapping_quality,
                                        fastafile=fasta_fake, start=position_list[0]-1, end=position_list[0])
    
        for pileupcolumn in pileupcolumns:
            reads_dict, lengths_dict, df_mapping_all, nm_tags_dict = min_over_reference_or_32del(pileupcolumn, reads_dict, lengths_dict, df_mapping_all, nm_tags_dict, length_threshold, ol_threshold, sample, position_list=position_list, pos=None, min_over_type='del')
        
        return reads_dict, lengths_dict, df_mapping_all, nm_tags_dict

    
def min_over_reference_or_32del(pileupcolumn, reads_dict, lengths_dict, df_mapping_all, nm_tags_dict, length_threshold, ol_threshold, sample, min_over_type, position_list, pos):
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
                
                if min_over_type == 'ref':
                    S, E = position_list
                    if reference_start <= S and reference_end >= E:
                        # Minimum overlapping length is 32
                        min_over = 32
                    else:
                        # I calculate the left and right overlap of the read
                        left = query_position
                        right = reference_end - pos + 1 
                        
                        # The minimum overlapping length is the minimum between these two
                        min_over = min(left, right)
                
                # I calculate the left and right overlaps of the read
                elif min_over_type == 'del':
                    left = query_position
                    right = reference_end - position_list[0] + 1

                    # I assign the minimum overlapping length and I add it with the relative read name to the dictionary
                    min_over = min(left, right)
                
                # I save all in a row
                row_to_add = [sample, read_name, reference_start, reference_end, read_sequence, read_length, min_over, nm_tag, min_over_type]

                df_length = len(df_mapping_all)
                
                # print("df_length", df_length)
                # That I add to a dataframe so I can analyse it afterwards
                # df_mapping_all.loc[df_length] = row_to_add

                if min_over >= ol_threshold:
                    # I append the minimum overlapping length to the reads dictionary in a list

                    to_add = [int(min_over)]
                    
                    if min_over_type == 'ref':
                        reads_dict[read_name].append(int(min_over))
                    elif min_over_type == 'del':
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

    # I create a list containing all the reads aligning in this region.
    # Each element of the list is a pysam.AlignmentSegment object
    aligned_list = [read for read in bamfile.fetch("CCR5_del32_120b.fasta", 30, 90)]
        
    # For each read name in reads_dict:
    for key in list(reads_dict):

        # For each read in the above list
        for read in aligned_list:

            # If the read from reads_dict is in the list AND does not have a perfect match, I'll remove it from the dictionary
            if key == read.query_name and read.get_tag("NM") != 0:
                del reads_dict[key]
                del lengths_dict[key]

    return reads_dict


def some_function(haplotype_list, haplo_df, ref_haplo_count_df, bamvsref, baq_snp, adjustment_threshold, length_threshold, sample):

    # I initialize dictionary where I'll store the reference and alternate bases called for each SNP
    dict_snps = OrderedDict()

    dict_ref_haplo_count = OrderedDict()
    referencecount, purereferencecount, haplocount, purehaplocount, notavail = 0, 0, 0, 0, 0
    # Iterate through each SNP 
    for snp in haplotype_list:

        chrom = str(3)
        idx, coordinate, ref, alt, rsquared = snp[0:5]
        coordinate = int(coordinate)
    
        # 1 - Extract all the reads bases mapping to each snp
        reads_list, other_list, ref_list, alt_list = extr_rbases_bam(bamvsref, chrom, coordinate, ref, alt, baq_snp, adjustment_threshold, length_threshold)

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


            dict_snps[idx + "_ref"] = [math.nan]
            dict_snps[idx + "_alt"] = [math.nan]

        else:

            dict_snps[idx + "_ref"] = [ref_bases_n]
            dict_snps[idx + "_alt"] = [alt_bases_n]

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