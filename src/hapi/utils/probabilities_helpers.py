import math
from typing import Tuple
import numpy as np
import pandas as pd


def phred2prob(x: int):
    """
    Convert Phred score to probability

    :param x: Phred score
    :return: probability
    """
    """"""
    return 10 ** (-x / 10)


def prob2phred(x: float) -> float:
    """Convert probability to Phred score"""
    return -10 * math.log10(x)


def pD_G_(reads_list: list, ref: str, alt: str) -> Tuple[float, float, float]:
    """
    Calculate step 1, 2, and 3 of genotype likelihood for each SNP
    :param reads_list: list of read bases mapping each of the CEU rs33
        haplotype SNPs; list
    :param ref: reference allele
    :param alt: alternate allele
    :return: pD_RR, pD_RA, pD_AA: genotypes likelihoods, i.e. probability of
        the data given each genotype
    """

    # Here with RR, RA, and AA I actually mean RR,RD,DD, so the deletion
    # genotypes p(D|RR), p(D|RA), p(D|AA)
    pD_RR, pD_RA, pD_AA = 0, 0, 0

    # If the list is not empty
    if reads_list != []:
        # Iterate over each read base of the list
        for base in reads_list:

            # Following Simon's slides sequence: Step 1: calculate prob of
            # each base given the observed allele, given that the true base
            # is the observed p_b_A = p(b|A)
            p_b_A_main = (1 - phred2prob(base[1]))
            p_b_A_rest = (phred2prob(base[1])) / 3

            # Step 2: calculate prob of each observed base given EACH
            # POSSIBLE genotype, so Ref/Ref, Ref/Alt, Alt/Alt. Step 3:
            # multiply (in this case sum, since they are logarithms),
            # over all bases to calculate the Likelihood of each genotype p(
            # D|G)
            if base[0] == ref:
                pD_RR += math.log10(p_b_A_main)  # Would be (p_b_A_main +
                # p_b_A_main)/2, so it simplifies to p_b_A_main
                pD_RA += math.log10((p_b_A_main + p_b_A_rest) / 2)
                pD_AA += math.log10(p_b_A_rest)

            elif base[0] == alt:
                pD_RR += math.log10(p_b_A_rest)
                pD_RA += math.log10((p_b_A_main + p_b_A_rest) / 2)
                pD_AA += math.log10(p_b_A_main)

    else:
        # If the reads_list is empty, i.e. there are no reads mapping that
        # position, set the probability of the data given the genotype (the
        # genotype likelihood) as 0.33, i.e. as random
        pD_RR = math.log10(0.33)
        pD_RA = math.log10(0.33)
        pD_AA = math.log10(0.33)

    return 10 ** pD_RR, 10 ** pD_RA, 10 ** pD_AA


# SNPs genotype calculation
# pD_G = p(D|G) likelihood
# pG_D = p(G|D)
# pD_ = p(D)
# p(G|D) = p(G) * p(D|G) / p(D)

# p(D|G)


# p(D)
def pD_(pD_RR: float, pD_RA: float, pD_AA: float) -> float:
    """
    p(D) = SUM p(Gi) p(D|Gi), for each genotype RR, RA, AA. Here we assume a
    uniform prior distribution, giving it a value of 0.33
    :param pD_RR: likelihood of getting the observed Data given the Genotype
        Ref Ref, calculated in the function pD_G_
    :param pD_RA: likelihood of getting the observed Data given the Genotype
        Ref Alt, calculated in the function pD_G_
    :param pD_AA: likelihood of getting the observed Data given the Genotype
        Alt Alt, calculated in the function pD_G_
    :return: pD: p(D)= SUM p(Gi) p(D|Gi)
    """

    pD = (0.33 * pD_RR) + (0.33 * pD_RA) + (0.33 * pD_AA)
    return pD


# p(G|D)
def pG_D_(pD_RR: float, pD_RA: float, pD_AA: float, pD: float,
          pG: float = 0.33) -> Tuple[float, float, float]:
    """
    Posterior probability of the specific base
    :param pD_RR: likelihood of getting the observed Data given the Genotype
        Ref Ref, calculated in the function pD_G_
    :param pD_RA: likelihood of getting the observed Data given the Genotype
        Ref Alt, calculated in the function pD_G_
    :param pD_AA: likelihood of getting the observed Data given the Genotype
        Alt Alt, calculated in the function pD_G_
    :param pD: p(D)
    :param pG: prior probability. Here we assume a uniform prior, so 0.33
    :return: pRR_D, pRA_D, pAA_D
    """
    # These are the probabilities of each genotype given the data FOR THIS BASE
    pRR_D = (pD_RR * pG) / pD
    pRA_D = (pD_RA * pG) / pD
    pAA_D = (pD_AA * pG) / pD

    return pRR_D, pRA_D, pAA_D


def prob_to_weighted(prob_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add, to the dataframe prob_df, columns containing the probabilities of each
        SNP multiplied by each relative R squared value
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


def calc_prob_joint(prob_df: pd.DataFrame) -> Tuple[
    float, float, float]:
    """
    Calculate Joint Posterior Probabilities of all the SNPs and normalize them
    :param prob_df:
    :return: pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm
    """

    log_list = np.array([prob_df["logP(RR|D)w"].sum(),
                         prob_df["logP(RA|D)w"].sum(),
                         prob_df["logP(AA|D)w"].sum()])

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

    # These will be used as the prior probability in the formula to
    # calculate the Posterior of the deletion and the reference

    pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm = normalized

    return pRR_D_joint_norm, pRA_D_joint_norm, pAA_D_joint_norm


def pD_RR_b_(pD_RR_g: float, pD_RR_d: float, pD_RD_g: float, pD_RD_d: float,
             pD_DD_g: float, pD_DD_d: float) -> Tuple[float, float, float]:
    """
    Joint likelihoods from both the bam files, vs GRCH37 and vs coll
    :return: pD_RR_b, pD_RD_b, pD_DD_b
    """
    pD_RR_b = pD_RR_g * pD_RR_d
    pD_RD_b = pD_RD_g * pD_RD_d
    pD_DD_b = pD_DD_g * pD_DD_d

    return pD_RR_b, pD_RD_b, pD_DD_b


# 32bp sequence posterior probabilities genotype calculation
# pD_G_2 = p(D|G) likelihood
# pD_2 = p(D)
# pG_D_2 = p(G|D)
# p(G|D) = p(G) * p(D|G) / p(D)

# p(D|G)
def p_D_G_2(reads_dict: dict, which_bam: str) -> Tuple[float, float, float]:
    """
    Calculate p(D|G), i.e. the probability of the data (the reads) given the
    Genotype, i.e. that the sample has RR, RD, or DD genotype
    :param reads_dict:
    :param which_bam: to specify which bam I want to analyze
    :return:
    """
    # p_ref_r = p(ref|r), i.e. the probability of having the reference
    # sequence given the observed read p_del_r = p(del|r), i.e. the
    # probability of having the deleted sequence given the observed read

    pD_RR_list, pD_RD_list, pD_DD_list = [], [], []

    # If the dict is not empty, i.e. if there is at least one read
    # overlapping the region
    if reads_dict != {}:
        for key, value in reads_dict.items():

            # Step 1: calculate prob of reference or deleted sequence given
            # the observed read
            if which_bam == "ref":
                p_ref_r = 1 - (1 / value) ** 2
                p_del_r = (1 - p_ref_r) / 2
            elif which_bam == "del":
                p_del_r = 1 - (1 / value) ** 2
                p_ref_r = (1 - p_del_r) / 2

            # Step 2: calculate prob of Data, of reads, given each possible
            # genotype, so Ref/Ref, Ref/Del, Del/Del
            # Step 3: multiply (in this case sum, since they are logarithms),
            # over all reads to calculate the Likelihood of each genotype
            # p(D|G). This would be like doing p_ref_r/2 + p_ref_r/2, and it
            # will just give p_ref_r so I simplified
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


# p(D)
def pD_2_(pRR_D_joint_norm: float, pRA_D_joint_norm: float,
          pAA_D_joint_norm: float,
          pD_RR_b: float, pD_RD_b: float, pD_DD_b: float) -> Tuple[
    float, float]:
    """
    p(D) = SUM p(Gi) p(D|Gi), for each genotype RR, RD, DD. As a prior I'll
    use the posterior probability for the rs333 haplotype that I calculated
    earlier """

    # Calculate with normalized
    pD_2_norm = (pRR_D_joint_norm * pD_RR_b) + \
                (pRA_D_joint_norm * pD_RD_b) + \
                (pAA_D_joint_norm * pD_DD_b)

    # Calculate it when it's random
    pD_2_r = (0.33 * pD_RR_b) + (0.33 * pD_RD_b) + (0.33 * pD_DD_b)

    return pD_2_norm, pD_2_r


# p(G|D)
def pG_D_2(pG: float, pD_G: float, pD: float) -> float:
    pG_D_2 = (pG * pD_G) / pD

    return pG_D_2
