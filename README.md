# Blend-Seq Paper

This repository contains the code used in the paper:

* Blended Length Genome Sequencing (blend-seq):
  Combining Short Reads with Low-Coverage Long Reads to Maximize Variant Discovery

A copy of the preprint can be found [here](https://www.biorxiv.org/content/10.1101/2024.11.01.621515v2). Please refer to the paper for more details on the methods and results.

## Data

The GATK-SV VCF used in the paper for a short-read control can be found under `data`. It comes from the following paper:
```
Zan Koenig, Mary T Yohannes, Lethukuthula L Nkambule, and others. A harmonized public 
resource of deeply sequenced diverse human genomes. Genome Research, 2024
```

## Code Overview

The code is split into `scripts` and `pipelines`. For documentation on what they do and how they were used, please see the Methods section of the manuscript.

## External Code

The following external code was used in the paper:
* [BenchmarkSVs](https://github.com/broadinstitute/palantir-workflows/blob/main/BenchmarkSVs/BenchmarkSVs.wdl) for benchmarking structural variants
* [BenchmarkVCFs.wdl](https://github.com/broadinstitute/palantir-workflows/blob/main/BenchmarkVCFs/BenchmarkVCFs.wdl) for benchmarking small variants
* [chop_reads](https://github.com/rickymagner/chop_reads) for chopping long read alignments into synthetic short read alignments
* [CleanSVs.wdl](https://github.com/broadinstitute/palantir-workflows/blob/main/BenchmarkSVs/CleanSVs.wdl) for cleaning SV VCFs to normalize for downstream tools