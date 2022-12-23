import argparse

def create_parser():

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
    
    return parser