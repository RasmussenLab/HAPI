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

    parser.add_argument("--samples-file", required=True,
                        type=argparse.FileType('r'),
                        help="File containing list of samples to process")
    parser.add_argument("--files-extension", required=True, type=str,
                        help="String describing the extension of the file, e.g. .rmdup.realign.md.cram or .rmdup.realign.bam")
    parser.add_argument("--folder-ref", required=True, type=pathlib.Path,
                        help="Folder path of the bam files aligned against the GRCH37 reference genome")
    parser.add_argument("--folder-coll", required=True, type=pathlib.Path,
                        help="Folder path of the bam files aligned against the Collapsed reference genome")

    parser.add_argument("--fasta-ref-file", required=True, type=pathlib.Path,
                        help="fasta file containing the reference genome")
    parser.add_argument("--fasta-coll-file", required=True, type=pathlib.Path,
                        help="fasta file containing the coll reference genome")

    parser.add_argument("--snps-file", required=True, type=pathlib.Path,
                        help="text file containing list of the 4 top SNPs")
    parser.add_argument("--haplotype-file", required=False, type=pathlib.Path,
                        help="text file containing list of the 86 SNPs")
    parser.add_argument("--output-folder", required=False, type=pathlib.Path,
                        default="./results/",
                        help="Output folder in which to append the results of the probability calculations")
    parser.add_argument("--baq-snps", type=parse_bool, default="False",required=False,
                        help="option to perform Base alignment quality (BAQ) in samtools for the reads mapping the 4 tag SNPs. Default is False")
    parser.add_argument("--baq-deletion", type=parse_bool, default="False", required=False,
                        help="option to perform Base alignment quality (BAQ) in samtools for the reads mapping to the deletion. Default is False")

    parser.add_argument("--length-threshold", type=int, required=True, default=1000,
                        help="value to keep only reads with read length < length_threshold. Not used in the final script, i.e. it is set to 1000 so no filtering is performed")
    parser.add_argument("--overlapping-length-threshold", type=int, required=False, default=6,
                        help="number of nucleotides with which each read needs to overlap the deletion or reference coordinates in order to be kept")
    parser.add_argument("--perfect-match", type=parse_bool, required=False, default=False,
                        help="option to keep only reads with perfect match to the  genome. Not used in the final script") 
    parser.add_argument("--adjustment-threshold", type=int, required=False, default = 50,
                        help="option adjust_capq_threshold of pysam to adjust mapping quality. Recommended value is 50 and is the default")

    return parser
