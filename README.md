# HAPI

Instructions to run HAPI (Haplotype-Aware Probabilistic model for Indels) to identify the CCR5delta32 deletion in ancient low coverage DNA samples, as published in the pre-print:

```
Tracing the evolutionary path of the CCR5delta32 deletion via ancient and modern genomes
Kirstine Ravn, Leonardo Cobuccio, Rasa Audange Muktupavela, Jonas Meisner, Michael Eriksen Benros, Thorfinn Sand Korneliussen, Martin Sikora, Eske Willerslev, Morten E. Allentoft, Evan K. Irving-Pease, Fernando Racimo, Simon Rasmussen
medRxiv 2023.06.15.23290026; doi: https://doi.org/10.1101/2023.06.15.23290026
```

The 144 ancient simulated DNA samples, together with the folder containing the results ran by HAPI, are available at [this link](https://doi.org/10.17894/ucph.a31d9052-546d-4f8f-8e16-e5bd896df67b).

After unzipping the file, HAPI can be installed and run with the following commands:

```
pip install hapi-pyth

mkdir results

hapi-pyth \
--samples-file list_samples.txt \
--files-extension .cram \
--folder-ref GRCh37 \
--folder-coll Collapsed \
--fasta-ref-file references/hs.build37.1.fa \
--fasta-coll-file references/ceuhaplo_collapsed.hs.build37.1.fa \
--snps-file top_4_snps.txt \
--length-threshold 1000 \
--output-folder results
```

HAPI will output several files in the results folder. The most important file is `results.tsv`, which is a table containing the prediction for each sample run.
```
* Sample: sample ID
* pRR_Data_n: Posterior probability of a sample being homozygous for the reference sequence, given the Data
* pRD_Data_n: Posterior probability of a sample being heterozygous for CCR5delta32, given the Data
* pDD_Data_n: Posterior probability of a sample being homozygous for CCR5delta32, given the Data
* N_reads_ref: Number of reads mapping to the reference sequence in the canonical reference
* N_reads_del: Number of reads mapping to the CCR5delta32 sequence in the collapsed reference
* Min_over_ref: List containing the minimum overlapping length of each read mapping to the canonical reference
* Min_over_del: List containing the minimum overlapping length of each read mapping to the collapsed reference
* Lengths_ref: List containing the length of each reads mapping to the canonical reference
* Lengths_del: List containing the length of each reads mapping to the collapsed reference
* Coverage_ref: Average coverage of reference SNPs for the four variants in the canonical reference
* Coverage_alt: Average coverage of alternate SNPs for the four variants in the canonical reference
* p_RR: Posterior probability of a sample having each of the top 4 variants in the SNP genotype ref|ref
* p_RA: Posterior probability of a sample having each of the top 4 variants in the SNP genotype ref|alt
* p_AA: Posterior probability of a sample having each of the top 4 variants in the SNP genotype alt|alt
* pData_RR: likelihood of the Data, given that the sample is homozygous for the reference sequence. Calculated as the joint likelihood from the alignments to the two references
* pData_RD: likelihood of the Data, given that the sample is heterozygous for CCR5delta32. Calculated as the joint likelihood from the alignments to the two references
* pData_DD: likelihood of the Data, given that the sample is homozygous for CCR5delta32. Calculated as the joint likelihood from the alignments to the two references
* pD_norm: marginal likelihood (denominator of the equation)
* pRR_Data_r: Posterior probability of a sample being homozygous for the reference sequence, given the random (uniform) haplotype
* pRD_Data_r: Posterior probability of a sample being heterozygous for CCR5delta32, given the random (uniform) haplotype
* pDD_Data_r: Posterior probability of a sample being homozygous for CCR5delta32, given the random (uniform) haplotype
* N_reads_mapping_both: Number of reads mapping both the canonical and the collapsed references
* SNP_1_rs113341849: Number of ALT alleles called for SNP 1
* SNP_2_rs113010081: Number of ALT alleles called for SNP 2
* SNP_3_rs11574435: Number of ALT alleles called for SNP 3
* SNP_4_rs79815064: Number of ALT alleles called for SNP 4
```

The file `settings.tsv` contains the options used when the HAPI analysis got run, and is useful for reproducibility purposes

The file `all_reads_mapping.tsv` contain all the reads mapping to either the canonical and the collapsed reference. 

The file `reads_assigned_ref.tsv` and `reads_assigned_del.tsv` contain the list of reads assigned to the canonical and collapsed reference, respectively, according to the minimum overlapping length option. This means that these two files will contain less reads than the previous one, i.e. only those that were assigned to the references, not merely all the reads mapping.

The folder `prob_dfs` contains, for each sample, the calls of each of the top 4 variants, together with how many calls, and the values of the probabilities for each variant

Please note that in its current state, HAPI can identify only the CCR5delta32 deletion. This is for two reasons:
1. The CCR5delta32 deletion has 4 equivalent representations, each with its own coordinates ([click here for details](https://varsome.com/variant/hg19/rs333?annotation-mode=germline)). HAPI was developed with these 4 different sets of coordinates in mind. Another deletion of interest might have only one representation or a set of different ones
2. HAPI uses the information from the top 4 tag variants in high LD with the CCR5delta32 as Prior. Another deletion of interest might not have known tag variants
  
Therefore, one could potentially extend HAPI to identify other deletions if they have information about tag variants in high LD with them, and adapting the code regarding the deletion coordinates.

For more details about HAPI, please refer to the pre-print references above.
