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


Please note that in its current state, HAPI can identify only the CCR5delta32 deletion. This is for two reasons:
1. The CCR5delta32 deletion has 4 equivalent representations, each with its own coordinates ([click here for details](https://varsome.com/variant/hg19/rs333?annotation-mode=germline)). HAPI was developed with these 4 different sets of coordinates in mind. Another deletion of interest might have only one representation or a set of different ones
2. HAPI uses the information from the top 4 tag variants in high LD with the CCR5delta32 as Prior. Another deletion of interest might not have known tag variants
  
Therefore, one could potentially extend HAPI to identify other deletions if they have information about tag variants in high LD with them, and adapting the code regarding the deletion coordinates.

For more details about HAPI, please refer to the pre-print references above.
