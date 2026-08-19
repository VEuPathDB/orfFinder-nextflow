# orfFinder-nextflow

A Nextflow pipeline that finds open reading frames (ORFs) in nucleotide sequences using a BioPerl-based ORF finder.

## Overview

This pipeline scans DNA sequences for open reading frames — searching for start/stop codon boundaries and translating candidate frames using a configurable genetic code table — and reports the qualifying ORFs (those meeting a minimum peptide length) as a sorted, tabix-indexed GFF file. It is used within VEuPathDB's genomic data workflows to generate ORF annotations, for example over genomic sequence for which gene models are not yet available. The input FASTA is split into subsets and processed in parallel, and the per-subset GFF results are merged, sorted, and indexed.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Docker or Singularity/Apptainer — processes run in the `bioperl/bioperl:stable` and `biocontainers/tabix` container images (select the engine via the `docker` or `singularity` profile/config in `conf/`)

## Usage

```
nextflow run VEuPathDB/orfFinder-nextflow -r main \
  --inputFilePath /path/to/sequences.fasta \
  --minPepLength 50 \
  --outputDir /path/to/output \
  -C conf/docker.config \
  -resume
```

The pipeline has a single (default) entry point:

1. `orfFinder` splits `params.inputFilePath` into subsets of `params.fastaSubsetSize` sequences and runs the `bin/orfFinder` BioPerl script on each subset, emitting a per-subset GFF of ORFs at least `params.minPepLength` amino acids long.
2. `indexResults` collects, sorts, and `bgzip`-compresses the merged GFF, then indexes it with `tabix`.

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `inputFilePath` | `data/input.fa` | FASTA file of nucleotide sequences to scan for ORFs |
| `minPepLength` | `50` | Minimum translated peptide length (amino acids) for an ORF to be reported |
| `fastaSubsetSize` | `1` | Number of sequences per chunk when splitting the input FASTA for parallel processing |
| `outputFileName` | `Orf50.gff` | Base filename for the sorted, indexed ORF GFF |
| `outputDir` | `$launchDir/output` | Directory the final indexed GFF is published to |

`bin/orfFinder` also supports `--startCodon`, `--stopCodon` (default `taa|tga|tag`), and `--translTable` (default `1`, the standard genetic code), though these are not currently exposed as pipeline-level params and would need to be added to the `orfFinder` process invocation to override.

## Output

- `<outputFileName>.gz` — sorted, `bgzip`-compressed GFF of predicted ORFs, published to `outputDir`
- `<outputFileName>.gz.tbi` — `tabix` index for the compressed GFF
