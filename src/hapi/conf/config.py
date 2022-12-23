import argparse
import pathlib

def parse_bool(s: str) -> bool:
    try:
        return {"true": True, "false": False}[s.lower()]
    except:
        raise argparse.ArgumentTypeError(f"expected true/false, got: {s}")

def create_parser():

    # ArgParse to parse the arguments given in the command-line
    parser = argparse.ArgumentParser()

    parser.add_argument("--samples-file", required=True, type=argparse.FileType('r'), help = "File containing list of samples to process")
    # parser.add_argument("--bam-ref", required=True,
    #                     help="write the input bam file aligned against the GRCH37 reference genome")
    # parser.add_argument("--bam-del", required=True, help="write the input bam file aligned against the 32del fake genome")

    parser.add_argument("--files-extension", required=True, type=str, help = "String describing the extension of the file, e.g. .rmdup.realign.md.cram or .rmdup.realign.bam")
    parser.add_argument("--folder-ref", required=True, type=pathlib.Path,
                        help="Folder path of the bam files aligned against the GRCH37 reference genome")
    parser.add_argument("--folder-fake", required=True, type=pathlib.Path,
                        help="Folder path of the bam files aligned against the 32del fake genome")


    parser.add_argument("--fasta-ref-file", required=True, type=pathlib.Path, help="fasta file containing the reference genome")
    parser.add_argument("--fasta-fake-file", required=True, type=pathlib.Path, help="fasta file containing the fake reference genome")

    parser.add_argument("--snps-file", required=True, type=pathlib.Path, help="text file containing list of the 4 top SNPs")
    parser.add_argument("--haplotype-file", required=False, type=pathlib.Path, help="text file containing list of the 86 SNPs")
    parser.add_argument("--output-folder", required=False, type=pathlib.Path, default="./results/",
                        help="write the output folder in which to append the results of the probability calculations")
    parser.add_argument("--baq-snps", type=parse_bool,  default="False", required=False)
    parser.add_argument("--baq-deletion",  type=parse_bool,  default="False", required=False)

    parser.add_argument("--length-threshold", type=int, required=True)
    parser.add_argument("--overlapping_length_threshold", type=int, default=4, help="overlapping length threshold")
    parser.add_argument("--perfect-match",  type=parse_bool, required=True)
    parser.add_argument("--adjustment-threshold", type=int, required=True)
    # baq_snp = "no"
    # baq_deletion = "no"
    
    return parser