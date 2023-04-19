# HIV-64148 Pipeline

THIS PIPELINE IS UNDER DEVELOPMENT, NO FEATURES ARE GUARANTEE TO WORK AT THIS MOMENT.

## Purpose

- Dependency management and installation of tools is a time consuming process.
- Allows scienctist to focus on science rather than software.
- This container contains software used in HIV quasispecies analysis pipeline
- The pipeline can be used for long-read only or hybrid analysis
- The pipeline produces quansispecie sequences, consensus sequence, and variant calling files

This is a software package incoperated with pipeline for an analysis of viral genome. It can performs reference based consensus assembly, variant calling and reference-free quasispecies assembly and identification on long-read data particularily from Oxford Nanopore technology. This software package aims to consolidate analytical tools needed for the analysis into a single, easy to use, portable Docker container which can be rapidly deployed on any machine regardless of its environment.

Wrapped around a software package is an analysis pipeline uses viral haplotypes assembler pipeline, [Strainline](https://github.com/HaploKit/Strainline), from [HaploKit](https://github.com/HaploKit) to produce quansispecie sequences witout the need of reference genome, this pipeline produces corrected read using local de bruijn graph construction method. These corrected reads then go into mulitple iteration of seed and cluster producing haplotigs which are then finalized into full-lenght haplotypes, To identify a subtype of each haplotype, the sequences are compared against a collection of selected 32 HIV strains consists of basic clades, CRFs and laboratory clones using BLAST, subtype of each haplotype is determined by a subtype of a sequence in the collection which returns a highest bitscore when aligned to the haplotype sequence. The result then used for selecting a reference genome from a list of 12 representative genome for each subtypes in order to perform variant calling using [Snippy](https://github.com/tseemann/snippy) from [tseemann](https://github.com/tseemann).
Finally, to perform sequence analysis the haplotype sequences are submitted to [Sierra Web Service 2](https://hivdb.stanford.edu/page/webservice/), a GraphQL-based web service for accessing [Stanford HIV database](https://hivdb.stanford.edu/) which returns mutation and drug resistant profiles.

## Installation

This pipeline can be download directly from GitHub in either ZIP file format or with `git clone https://github.com/minaminii/HIV-64148.git`. It is recommended that you use this pipeline within the provided Docker environment to ensure its functionality. Docker daemon is required for building the environment, it can be downloaded from [Docker website](https://docs.docker.com/get-docker/).

```shell
# On host environment
# Clone the repository
git clone https://github.com/minaminii/HIV-64148.git
cd HIV-64148.git
# Build a Docker image
docker build -t hiv64148:latest -f ./docker/Dockerfile .
# Mount directory and run Docker image in interactive mode
docker run -it \
-v $(pwd)/workspace:/workspace \
hiv64148:latest bash
```

After the image is built, the required dataset representative sequences of each HIV-1 subtypes and a set of 32 HIV-1 sequences for building BLAST database, must be obtained before running the analysis. This process can be done by running `python3 /hiv64148/scripts/main.py setup` in the Docker environment. After the setup is finished, the analysis pipeline can be run with command `python3 /hiv64148/scripts/main.py run`, with argument `-i` as raw reads input in FASTA or FASTQ format; `-o` as a path to output directory. Additionally the pipeline can be run with command `python3 /hiv64148/scripts/main.py run_cli` to start a command-line input dialogue.

```shell
# Inside Docker environment
# Obtain selected sequences for running the analysis
python3 /hiv64148/scripts/main.py setup

python3 /hiv64148/scripts/main.py run \
-i ${input_reads} \
-o ${output_directory}
```

OR run the Command-line input dialogue version.

``` shell
python3 /hiv64148/scripts/main.py run_cli
# > Path to FASTQ: ${input_reads}
# > Output directory: ${output_directory}
```

## The pipeline

We only require user to specified input and output directory. The configuration file comes with a preconfigured settings located in `settings/settings.json` changing a configuration is optional, beware that changing keys and , both dependencies and installation of the analysis softwares are managed by Python provided with the pipeline.

### Configuration

The pipeline configuration file is in a JSON format, located in a config directory, users can customize parameters as needed for each run.

#### Input

You can place your input into `./input` directory. On this current version, custom file name is not supported, please name your input as …..

#### Output

Your output will be present in `./output` directory after the pipline finished. Important output will be structured in directories as follow.

```ascii
${output_directory}/
├─ NanoPlot_QC/
│  ├─ NanoPlot-report.html
│  └─ NanoStats.txt
├─ Variants/
│  └─ hapN_vs_subtype_X/
│     ├─ snps.aligned.fasta
│     ├─ snps.bam
│     ├─ snps.bam.bai
│     ├─ snps.bed
│     ├─ snps.consensus.fasta
│     ├─ snps.vcf.gz
│     └─ snps.vcf.gz.csi
├─ haplotype.blast.csv
├─ haplotypes.final.fasta
├─ haplotype_divergence.txt
└─ hiv-64148_report.html
```

#### BLAST database

The database was created using selected laboratory clones, subtypes, and CRFs of 32 HIV strains listed in the supplementary ___ . This database is used to identified quasispecies produced by Strainline via BLASTN software.

### Tools

- Strainline
  - d'accord v0.0.10
  - spoa v4.0.9
  - jgi_summarize_bam_contig_depths from metabat2
- Minimap2 v2.24
- Samtools v.1.10
- Snippy v4.6.0
- BLASTN v2.13.0+

### Pipeline flowchart



## To optimize

- BLAST+ produces a report with cross-alignment data of all quasispecies vs database which makes the report cluttered and difficult to interprete.
- Interactive command-line based software to control custom input/output and configuration (currently I/O require user to change their name into specified name)

## Acknowledgement

This project is supported by The Health System Research Institute, Thailand (Grant HSRI-64-148)
