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
import csv
from itertools import chain
from time import time

import pandas as pd

from hapi.conf.config import create_parser
from hapi.utils.data_utils import (
    open_files_args,
    snp_haplo_list,
    dict_to_list,
    write_probdf,
    write_results
)
from hapi.utils.probabilities_helpers import (
    prob_to_weighted,
    calc_prob_joint,
    pD_RR_b_,
    p_D_G_2,
    pD_2_,
    pG_D_2
)
from hapi.utils.mappings_helpers import (
    calc_snps_posteriors,
    minimum_overlap,
    average_minimum_overlap,
    tag_filtering,
    some_function,
    remove_overlaps
)

############## Execution #################

def main():
    # Start time initiation
    start = time()
    
    # Get script arguments
    parser = create_parser()
    args = parser.parse_args()
    
    # Output folder
    results_filepath = args.output_folder / "results.tsv"
    outdir = args.output_folder / "prob_dfs/"
    outdir.mkdir(exist_ok=True, parents=True)
    
    # Initialize empty list and header assignments
    mapping_all = []
    header, header_del, header_ref = True, True, True
    
    # Open the samples list
    samples_list = args.samples_file.read().splitlines()

    # If --haplotype option is activated --> write a table to file containing the reporting of all the 86 SNPs
    if args.haplotype_file:
         # Initialize empty list and header assignments
        header_hapl, header_hapl_c = True, True
        
        # e.g. [['rs58697594', '46275570', 'G', 'A', '0.8602'], ['rs73833032', '46276490', 'T', 'C', '0.8602']]
        # Need to convert this list of lists in another list of just the column names with the format
        # rs58697594_ref, rs58697594_alt, rs73833032_ref, rs73833032_alt
        # I contruct the rs_list as a list containing all the 82 SNPs _alt and _ref
        
        haplotype_list = snp_haplo_list(args.haplotype_file)
        rs_list = list(chain.from_iterable((rs_id[0] + "_alt", rs_id[0] + "_ref") for rs_id in haplotype_list))
        rs_list.sort()


    # For each sample to analyse
    for sample in samples_list:
        print(sample)
        
        # I parse the arguments given when executing the script
        bamvsref, bamvsdel, fasta_ref, fasta_fake = open_files_args(args, sample)

        ############## Part 0: If --haplotype option is activated --> write a table to file containing the reporting of all the 86 SNPs
        if args.haplotype_file:

            # I report all the SNPs called of the haplotype
            haplo_results_list, ref_haplo_count_list = some_function(haplotype_list, bamvsref, args.baq_snps, args.adjustment_threshold, args.length_threshold, sample, fasta_ref)


        ############## Part A: Prior Probability calculated as joint probability of top 4 SNPs' posterior probabilities - Execution #################

        # 1 - Extract SNPs from file
        # List of lists containing the 4 SNPs of the CEU rs333 haplotype with coordinates and R squared value to the rs333
        # e.g. [['rs58697594', '46275570', 'G', 'A', '0.8602'], ['rs73833032', '46276490', 'T', 'C', '0.8602']]
        snp_list = snp_haplo_list(args.snps_file)

        # 2 - Calculate Posterior probability of each SNP given each possible Genotype
        prob_df, coverage_ref, coverage_alt, coverage_other, dict_snps_cov = calc_snps_posteriors(snp_list, bamvsref, fasta_ref, args.baq_snps, args.adjustment_threshold, args.length_threshold)


        # 3 - Weighted probs calculation: multiply each SNP's probability by the relative Rsquared and store in new columns.
        # 3 - Plus, add a column with their logarithm10 conversion
        prob_df = prob_to_weighted(prob_df)

        # I write the prob_df as a tab separated file
        write_probdf(prob_df, outdir, sample)

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
        reads_dict_ref, lengths_dict_ref, mapping_all, nm_tags_dict_ref = minimum_overlap(bamvsref, "3", position_list_reference, args.adjustment_threshold, mapping_all, args.length_threshold, args.overlapping_length_threshold, sample, fasta_fake, fasta_ref, baq = args.baq_deletion)
        
        reads_dict_del, lengths_dict_del, mapping_all, nm_tags_dict_del = minimum_overlap(bamvsdel, "3", position_list_deletion, args.adjustment_threshold, mapping_all, args.length_threshold, args.overlapping_length_threshold, sample, fasta_fake, fasta_ref, baq = args.baq_deletion)

        # 3 - Average of the overlapping lengths of all the 4 coordinates couples in the bam vs GRCh37
        reads_dict_ref = average_minimum_overlap(reads_dict_ref)


        # In case there are reads that overlap both the reference and the collapsed genome, I'll keep only the one that
        # has the lowest number of mismatches and the highest overlapping length                            
        reads_dict_del, reads_dict_ref, nm_tags_dict_del, nm_tags_dict_ref, lengths_dict_ref, lengths_dict_del, N_reads_mapping_both = remove_overlaps(reads_dict_del, reads_dict_ref, nm_tags_dict_del, lengths_dict_ref, lengths_dict_del, nm_tags_dict_ref)

        # 4 - Filter by the XM:i:0 tag --> keep only the perfect matching reads

        if args.perfect_match:
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
        pD_RR_b, pD_RD_b, pD_DD_b = pD_RR_b_(pD_RR_g, pD_RR_d, pD_RD_g, pD_RD_d, pD_DD_g, pD_DD_d)


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

        
        # 11 - Make records
        
        records = [{
            "Sample": sample,
            "pRR_Data_n": pRR_D_2_norm,
            "pRD_Data_n": pRD_D_2_norm,
            "pDD_Data_n": pDD_D_2_norm,
            "N_reads_ref": len(reads_dict_ref),
            "N_reads_del": len(reads_dict_del),
            "Min_over_ref": reads_list_ref,
            "Min_over_del": reads_list_del,
            "Lengths_ref": lengths_list_ref,
            "Lengths_del": lengths_list_del,
            "Coverage_ref": coverage_ref,
            "Coverage_alt": coverage_alt, 
            "SNP_1_rs113341849": dict_snps_cov["rs113341849"],
            "SNP_2_rs113010081": dict_snps_cov["rs113010081"],
            "SNP_3_rs11574435": dict_snps_cov["rs11574435"],
            "SNP_4_rs79815064": dict_snps_cov["rs79815064"],
            "p_RR": pRR_D_joint_norm,
            "p_RA": pRA_D_joint_norm,
            "p_AA": pAA_D_joint_norm,
            "pData_RR": pD_RR_b,
            "pData_RD": pD_RD_b,
            "pData_DD": pD_DD_b,
            "pD_norm": pD_2_norm,
            "pRR_Data_r": pRR_D_2_r,
            "pRD_Data_r": pRD_D_2_r,
            "pDD_Data_r": pDD_D_2_r,
            "N_reads_mapping_both": N_reads_mapping_both
        }]
        
            
        records_ref, records_del = [], []
        ref_del_references = zip([reads_dict_ref, reads_dict_del], ["ref", "del"], [records_ref, records_del])
        for reads_dict, _class, records_class in ref_del_references:
            for read in reads_dict.keys():
                result_class = {
                    "sample": sample, 
                    "read_name": read,
                    "class": _class
                }
            
                records_class.append(result_class)  
        
        # 11 - Append the results to the output file
        
        header = write_results(results_filepath, records, header)
        header_ref = write_results(args.output_folder / "reads_assigned_ref.tsv", records_ref, header_ref)
        header_del = write_results(args.output_folder / "reads_assigned_del.tsv", records_del, header_del)
        if args.haplotype_file:
            header_hapl = write_results(args.output_folder / "SNPS_reporting.tsv", haplo_results_list, header_hapl) 
            header_hapl_c = write_results(args.output_folder / "ref_haplo_counts.tsv", ref_haplo_count_list, header_hapl_c)

    # I average the overlapping lengths 
    df_mapping_all = pd.DataFrame.from_records(mapping_all)
    df_mapping_all = (df_mapping_all
    .assign(average_min_over = 
        lambda x: 
            x.groupby(["sample","read_name","alignment"])
            ["min_over"].transform("mean")))

    df_mapping_all = df_mapping_all.drop(columns = ["min_over"])
    df_mapping_all = df_mapping_all.drop_duplicates()

    # I need to average the overlapping lengths of the ref
    df_mapping_all.to_csv(args.output_folder / "all_reads_mapping.tsv", sep="\t", quoting=csv.QUOTE_NONE, index = False)

        
    end = time()
    length = end - start
    print("Time:", length)

    
if __name__ == "__main__":
    main()