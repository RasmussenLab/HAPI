import pysam
from collections import OrderedDict
from pathlib import Path

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
    reads_list = []

    # I convert the reads_dict from a normal dict to an Ordered Dict
    reads_dict = OrderedDict(reads_dict)
    if reads_dict != {}:
        
        for key, value in reads_dict.items():
            reads_list.append(value)

    return reads_list

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


def write_results(results_filepath, coverage_ref, coverage_alt, dict_snps_cov,
                  sample, N_reads_mapping_both, pRR_D_2_norm, pRD_D_2_norm,
                  pDD_D_2_norm, reads_dict_ref, reads_dict_del, reads_list_ref,
                  reads_list_del, lengths_list_ref, lengths_list_del,
                  pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm,
                  pD_RR_b, pD_RD_b, pD_DD_b, pD_2_norm, pRR_D_2_r, pRD_D_2_r,
                  pDD_D_2_r):
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