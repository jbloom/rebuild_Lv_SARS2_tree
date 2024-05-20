# Interactive phylogeny of early SARS-CoV-2 that attempts to reproduce Fig 3a of 2024 study by group of Yong-Zhen Zhang

## Overview
This repository builds a phylogeny designed to reproduce that in Fig. 3a of [Lv et al, Virus Evolution (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252), which reports a phylogeny of SARS-CoV-2 sequences in previously in public databases and new sequences collected by Zhang's group from what Zhang calls "stages 0 and 1" of the outbreak (before March-1-2020).

The phylogeny here was built by trying to exactly reproduce the methods of [Lv et al, Virus Evolution (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252) used to make Fig. 3a of that paper, so see that paper for details about the sequences included, their lineage designations, rooting of the tree, etc.
However, the methods were not sufficient to completely reproduce every detail, so the tree here is very close but may not be exactly identical to that from the paper.
In addition to the metadata from [Lv et al, Virus Evolution (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252), sequences noted as being from 2019 in Tables 6 and 7 of the [joint WHO-China report](https://www.who.int/publications/i/item/who-convened-global-study-of-origins-of-sars-cov-2-china-part) are also annotated on this phylogeny.

Note that the best way to root the SARS-CoV-2 phylogeny is still an open topic with substantial uncertainty (see [Pipes et al (2021)](https://academic.oup.com/mbe/article/38/4/1537/6028993)).
There are also open questions about the best way to subsample, de-duplicate, and quality control early SARS-CoV-2 sequences. Further, the high similarity of the sequences mean there is also just limited statistical support to reliably draw conclusions about some aspects of the phylogeny (see the [paper "Phylogenetic analysis of SARS-CoV-2 data is difficult" by Morel et al (2020)](https://academic.oup.com/mbe/article/38/5/1777/6030946)).

This phylogeny attempts to reproduce the one from the study by Zhang's group ([Lv et al, Virus Evolution (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252)) since that is the most recent major study on the topic and includes additional sequence data not in prior analyses.
However, further work on rooting and sequence curation are needed (this phylogeny should help with that), and recall the conclusion of [Pipes et al (2021)](https://academic.oup.com/mbe/article/38/4/1537/6028993) that _"These results suggest that phylogenetic evidence alone is unlikely to identify the origin of the SARS-CoV-2 virus and we caution against strong inferences regarding the early spread of the virus based solely on such evidence."_

## Details

The goal of this repository is to re-build the portion of the phylogenetic tree covering what the paper refers to as Stage 0 and Stage 1 of the outbreak, which is everything prior to March-1-2020 so that it can be rendered interactively using the Nextstrain [auspice.us](https://auspice.us/) platform.
This is the tree shown in Figure 3A of the paper.

Because that tree is not provided in raw form in the paper, this repository tries to repeat the analysis in [Lv et al (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252) to produce it.

The `conda` environment to run the analysis is in [environment.yml](environment.yml).

The input data are in [./data/](data):
  - [data/Supplementary Tables 1-3.xlsx](data/Supplementary Tables 1-3.xlsx): supplementary tables from [Lv et al (2024)](https://academic.oup.com/ve/article/10/1/veae020/7619252)
  - [data/extract_accessions.ipynb](data/extract_accessions.ipynb): Jupyter notebook that extracts lists of accessions from  Supplementary Table 3. It should be run manually to create the following three files:
    - [data/all_accessions.csv](data/all_accessions.csv): information on all accessions in Supplementary Table 3
    - [data/genbank_accessions.csv](data/genbank_accessions.csv): all Genbank accessions in Supplementary Table 3
    - [data/gisaid_accessions.csv](data/gisaid_accessions.csv): all GISAID accessions in Supplementary Table 3
  - [data/genbank_sequences.fa](data/genbank_sequences.fa): all sequences specified by [data/genbank_accessions.csv](data/genbank_accessions.csv) as manually downloaded from Genbank.
  - [data/gisaid_sequences.fa](data/gisaid_sequences.fa): all sequences specified by [data/gisaid_accessions.csv](data/gisaid_accessions.csv) as manually downloaded from GISAID. **Due to GISAID data sharing restrictions, this sequence file is not actually tracked in the repository. You will need to manually generate it yourself from the list of accessions.**
  - [data/reference_seq.gb](data/reference_seq.gff): SARS-CoV-2 Genbank file for calling amino-acid mutations from [https://github.com/nextstrain/ncov/blob/master/defaults/reference_seq.gb](https://github.com/nextstrain/ncov/blob/master/defaults/reference_seq.gb).
  - [data/jointWHO_sequences.csv](data/jointWHO_sequences.csv): sequences from the Tables 6 and 7 of the [joint WHO-China report](https://www.who.int/publications/i/item/who-convened-global-study-of-origins-of-sars-cov-2-china-part), which are stated to be the only sequences from patients with onset date before Dec-31-2019.

The pipeline to build the tree mostly uses the [augur](https://docs.nextstrain.org/projects/augur) pipeline, and is run by the `snakemake` file in [Snakefile](Snakefile).

All the files created by the pipeline are placed in `./results/` subdirectory which is not tracked in this GitHub repo due to GISAID data-sharing restrictions.
The final tree is in [auspice/rebuild_Lv_SARS2_tree.json](auspice/rebuild_Lv_SARS2_tree.json), and can be viewed with Nextstrain using the [instructions here](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html) at the link [https://nextstrain.org/community/jbloom/rebuild_Lv_SARS2_tree](https://nextstrain.org/community/jbloom/rebuild_Lv_SARS2_tree).
