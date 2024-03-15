# Phascolarctobacterium_CRC

This repository contains code for data analysis used to generate results in the paper *Species-level verification of Phascolarctobacterium association to colorectal cancer*.

## Data generation and analyses
The script used for preparation of datasets is the. The script used for analyses of the data is the [/scripts/data_generation.R](https://github.com/Rounge-lab/Phascolarctobacterium_CRC/blob/main/scripts/data_analyses.R). Required R packages include [tidyverse](https://www.tidyverse.org/packages/), [rstatix](https://cran.r-project.org/web/packages/rstatix/index.html), [vegan](https://cran.r-project.org/web/packages/vegan/index.html), [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html) and [MicrobiomeProfiler](https://bioconductor.org/packages/release/bioc/html/MicrobiomeProfiler.html).

In this study we used four different cohorts to verify a previously reported association between Phascolarctobacterium species and colorectal cancer found in [Bucher-Johannessen et al. (2023)](https://www.ncbi.nlm.nih.gov/pubmed/37182146) and [Senthakumaran et al. (2023)](https://www.ncbi.nlm.nih.gov/pubmed/36703031). One of these cohorts was the publicly available [CuratedmetagenomeData](https://waldronlab.io/curatedMetagenomicData/articles/curatedMetagenomicData.html). 

## Pangenome analyses
Pangenome analyses of the CRCbiome samples was generated using a snakemake pipeline with the script [pangenome.smk](https://github.com/Rounge-lab/Phascolarctobacterium_CRC/blob/main/pangenome.smk). Average nucleotide identity was estimated using the python script [genomes_for_ANI.py](https://github.com/Rounge-lab/Phascolarctobacterium_CRC/blob/main/scripts/genomes_for_ANI.py).

See our webpage for more information about the [CRCbiome](https://www.mn.uio.no/sbi/english/groups/rounge-group/crcbiome/) and the [NORCCAP](https://www.kreftregisteret.no/Forskning/Prosjekter/NORCCAP/) study 
