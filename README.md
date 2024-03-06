# HIV-64148 Pipeline

This is a software package pipeline for an analysis of HIV-1 genome. It can performs de novo or reference-based quasispecies assembly and identification on long-read data particularily from Oxford Nanopore technology. This package aims to consolidate analytical tools needed for such analysis into a single, easy to use, portable Docker container, capable of functioning across different computational environments with minimal setup.

Within this pipeline, we've configured the installation of these following long-read assembers, (1) de novo long-read assemblers, which include [Canu](https://github.com/marbl/canu), [Goldrush](https://github.com/bcgsc/goldrush), [Flye](https://github.com/fenderglass/Flye), and [Strainline](https://github.com/HaploKit/Strainline), and (2) reference-based long-read assemblers, which comprise [HaploDMF](https://github.com/dhcai21/HaploDMF), [iGDA](https://github.com/zhixingfeng/iGDA), and [RVHaplo](https://github.com/dhcai21/RVHaplo).

The Docker file comes with [Strainline](https://github.com/HaploKit/Strainline) and [MetaFlye](https://github.com/fenderglass/Flye) installed. Other assemblers mentioned in the article are installed as needed to minimize the container size. \(the container may take up to 15Gb for the full installation\).


## Usage

### Minimum system requirements

Based on our benchmarking study, the following optimal computational requirements are suggested:

- All assemblers should be executed in a Linux-based operating system.
- A quad-core (4 cores) or higher CPU with x86 microarchitecture is recommended,
- An 8 GB memory is sufficient but at least 16Gb is recommended.
    - Note: Strainline requires at least 16Gb baseline to run on a single threaded mode (`-t 1`) and scaled up roughly 4Gb per core.


### Installation

This pipeline can be download directly from GitHub with `git clone https://github.com/minaminii/HIV-64148.git`. It is recommended that you use this pipeline within the provided Docker environment to ensure its functionality. Docker daemon is required for building the environment, it can be downloaded from [Docker website](https://docs.docker.com/get-docker/).

```shell
# On host environment
# Clone the repository
git clone https://github.com/minaminii/HIV-64148.git
cd HIV-64148.git
# Build a Docker image
docker build -t hiv64148:latest -f ./docker/Dockerfile .
# Print help message
docker run --rm hiv64148:latest hiv64148 -h
```

> For using with Singulararity user have to build the Docker image and convert it to a Singularity image `singularity build hiv64148.sif docker-daemon://local/hiv64148:latest`

### Basic
After the setup is finished, you can attach a local directory to the container and run the pipeline with command `python3 /hiv64148/scripts/main.py run`, with argument `-i` or `--input` as raw reads input in FASTA or FASTQ format; `-o` or `--output_dir` as a path to output directory. 

```shell
docker run \
    -v ${YOUR_WORK_DIR}:/workspace \
    --rm hiv64148:latest \
    hiv64148 run \
        -i /workspace/${YOUR_FASTQ} \
        -o /workspace/${YOUR_OUTPUT}
```
### Switching between assemblers
To swith between differnet assemblers, user have to specify an `-a` or `--assember` argument with the desired assembler of choice `{canu,strainline,goldrush,metaflye,rvhaplo,haplodmf,igda} [default: strainline]`.

```shell
docker run \
    -v ${YOUR_WORK_DIR}:/workspace \
    hiv64148:latest \
    hiv64148 run \
        -i /workspace/${YOUR_FASTQ} \
        -o /workspace/${YOUR_OUTPUT} \
        --assember canu
```

> The created environment is not persistance, if you need to reuse the assembler installed in the container pleae use [`docke exec`](https://docs.docker.com/reference/cli/docker/container/exec/) to execute command on the container.

### Customize assembler parameters

Customized parameters for each assembler can be pass through with `-ag=` or `--assember-args=` argument, for the supported assembler argument, please refer to the documentation of your selected assembler.

```shell
docker run \
    -v ${YOUR_WORK_DIR}:/workspace \
    hiv64148:latest \
    hiv64148 run \
        -i /workspace/${YOUR_FASTQ} \
        -o /workspace/${YOUR_OUTPUT} \
        --assember strainline \
        --assember-args="--minTrimmedLen 500 --minOvlpLen 1000 -t 1"
```

## The pipeline

An implementation of HIV-64148 pipeline for HIV-1 genomic surveillance. The pipeline comprises three main stages: a read quality control analysis of long-read FASTQ files, followed by assembly using either de novo or reference-based assemblers, and concluding with the identification of HIV-1 subtype and drug resistance analysis.

To identify a subtype of each haplotype, the sequences are compared against a collection of selected 32 HIV strains consists of basic clades, CRFs and laboratory clones using BLAST, subtype of each haplotype is determined by a subtype of a sequence in the collection which returns a highest bitscore when aligned to the haplotype sequence.

Finally, to perform sequence analysis the haplotype sequences are submitted to [Sierra Web Service 2](https://hivdb.stanford.edu/page/webservice/), a GraphQL-based web service for accessing [Stanford HIV database](https://hivdb.stanford.edu/) which returns mutation and drug resistant profiles.

#### BLAST database

A selection of 12,000 HIV-1 genomes from the Los Alamos HIV-1 sequence database was included within the container to be used as a custom database for a local NCBI BLASTN used for HIV-1 subtype identification of the assembled quasispecies.

#### Output

Your output will be present in the directory specified with `-o` argument. The output will be structured in directories as follow.

```ascii
${output_directory}/
├─ NanoPlot_QC/
│  ├─ NanoPlot-report.html
│  └─ NanoStats.txt
├─ haplotype.blast.csv
├─ haplotypes.final.fasta
└─ hiv-64148_report.html
```

## Acknowledgement

The authors would like to thank The ERAWAN HPC Service, Chiang Mai University, Thailand, and The Center of Multidisciplinary Technology for Advanced Medicine (CMUTEAM), Faculty of Medicine, Chiang Mai University, Thailand for their supports on computational resources and server. We would like to thank Support the Children Foundation, Chiang Mai, Thailand for financial support. We also would like to thank Yuphin Chromwinya and Somporn Sankonkit of The Immunology Lab, Department of Microbiology, Faculty of Medicine, Chiang Mai University, Thailand for their administrative supports.

### Funding

This work was supported by the following funding bodies:

- The Health Systems Research Institute (Grant No. 64-148)
- The Faculty of Medicine Research Fund, Chiang Mai University (Grant No. 099-2563)
- Support the Children Foundation, Chiang Mai, Thailand.

## Citation

____________________________________
