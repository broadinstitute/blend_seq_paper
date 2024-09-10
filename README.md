# Blend-Seq Paper

This repository contains the code used in the paper:

- Blended Length Genome Sequencing (blend-seq): blending short and long reads to maximize variant discovery.

A copy of the manuscript can be found at **(insert link)**

## Code Overview

The code is split into `scripts` and `pipelines`.

### Scripts

The scripts included are:
- `convert_large_indels.py`: Takes a BAM as input (first argument), and writes a new BAM (second argument) with deletions/insertions of length above a threshold (50) converted to `N`s.
- `TrivialCaller.java`: Takes a BAM and produces a VCF by taking all deletions/insertions above a threshold (50) in each alignment and writing them as variants.

### Pipelines

The pipelines included are:
- `VgCallSVs.wdl`: Takes in an `augmentation_bam` file to use to augment linear reference into graph reference using variants produced by running the `TrivialCaller.java`. Optionally takes in a `custom_augmentation_calls` VCF instead for specific variants to use to add paths to the reference graph. After it takes a fastq (or two for paired reads) and aligns them to the graph reference. Then it calls structural variants using `vg call` and outputs a VCF.