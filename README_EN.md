# Metagen phage pipeline

The goal of this repository is to demonstrate how a pipeline for recovering Bacteriophages from metagenome samples can be made.

This pipeline was built using Python3 to build a command line interface (CLI) and the Workflow Description Language (WDL) to build the pipeline.

## Pipeline description

Two _inputs_  are required to run the pipeline:

- __Reads:__ _Reads_ used to assemble the contigs of interest. This _reads_ needs to be _paired-end_ and named with the suffixes  "_1.fastq" and "_2.fastq";
- __Multifasta__: Multifasta file containing all the contigs of interest that have been assembled from the __reads__ provided.

The steps below are developed throughout the pipeline.

### Dereplication

Dereplication is the step of removing the same _contigs_ built during assembly. For this step we used the [dedupe from bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/). This tool removes all _contigs_ that are identical (also considering the complementary reverse).

### Mapping

After dereplication, the next stage is mapping the _contigs_ against the _reads_. This step is required in order to calculate the coverage of each _contig_. This information is crucial for the _binning_ process (See [Basics of metagenomic](https://www.youtube.com/watch?v=MqD4aN1p1qA)). This process was carried out using [BWA](https://bio-bwa.sourceforge.net/).

### Sorting

The SAM archives created in the mapping stage were transformed into BAM files (binary form of SAM) and then _sorted_. This whole process was carried out using the [samtools](http://www.htslib.org/).

### Binning

In this stage we will assemble the _bins_, also called _Metagenome Assembled Genomes_ (MAGs). This assembly was done using the [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) with some parameters changed to improve the _binning_ of viral MAGs (-m 1500 -s 10000).

### Prediction

The final step is the prediction of bacteriophages in the assembled MAGs. The prediction tool [VirSorter2](https://github.com/jiarong/VirSorter2) was used for this task. Before prediction, the pipeline concatenates MAGs contigs, placing Ns between them. The MAGs has its _header_ renamed to the name of the file and then is submitted for prediction. Since the interest of this pipeline is the recovery of MAGs from phages, only the prediction group `dsDNAphage` was used.

After prediction, the results can be analyzed and interpreted to identify and characterize the bacteriophages present in the metagenomes. This information is valuable for understanding microbial ecology and the interaction between bacteria and bacteriophages in different environments.

## Installation

O primeiro passo é realizar o download do repositório:

```console
$ git clone https://github.com/visflores/metagen-phage-pipeline.git
```

Navigate to the downloaded directory and create a python virtual environment (v3.10 or higher):

```console
$ cd metagen-phage-pipeline
$ python3 -m venv env
```

With the virtual environment created, we can start installing the packages needed to run the pipeline:

```console
$ source env/bin/activate
$ (env) pip install .
```

To check that the pipeline is installed correctly, run the test script:

```console
$ (env) cd tests/
$ (env) python3 -m unittest -v
```

If the execution is finished without any errors, the pipeline is ready for use.

## Running the pipeline

Use the command below to check all CLI options:

```console
$ (env) metagen-phage -h

Pipeline para a recuperação de MAGs de bacteriófagos a partir de metagenomas.

options:
  -h, --help            show this help message and exit
  -file SCAFFOLD_FILE, --scaffold-file SCAFFOLD_FILE
                        Arquivo multifasta com todos os contigs montados a partir dos reads passados a flag '--reads-
                        dir'.
  -reads READS_DIR, --reads-dir READS_DIR
                        "Diretório contendo os reads usados para a montagem dos contigos do arquivo multifasta
                        passado para a flag '--scaffold-file'.
  -t THREADS, --threads THREADS
                        Quantidade de threads que será usada para a execução da pipeline.
  -pro, --prophage-find
                        Ativa a busca por profagos nos MAGs montados. Ativar essa opção torna o processo mais
                        demorado.

```

Two parameters are essential to run the pipeline:

- `-file`: Multifasta file containing all _contigs_ of interest and which have been assembled using the _reads_ pointed in the `-reads` parameter;
- `-reads`: Directory containing the _reads_ used to assemble the _contigs_ in the multifasta file. Only _reads paired end_ are accepted by the pipeline. The _reads_ must be named according to the [NAME]_1.fastq and [NAME]_2.fastq rules. The _reads_ must also be decompressed before use in the pipeline.

Use the command below to run the pipeline:

```console
$ (env) metagen-phage -file multifasta.fasta -reads reads/ -t 8
```

The results of the pipeline executions will be placed in a directory called `results` (if the directory doesn't exist, it will be created). 

At each run, a new directory is created within `results` containing the _outputs_ generated (these directories are named according to the time of execution), and the MAGs generated can be found in the `out` folder of these directories. 

## Case Study

In order to demonstrate the efficiency of this pipeline, a brief case study was carried out using real data. This data comes from the metagenome of newborn babies and was explored in the article by [SHARON, Itai et al](https://pubmed.ncbi.nlm.nih.gov/22936250/).

### Recovering _reads_ and assembling _contigs_

The initial step in this study was to assemble the _reads_ so that we can get the _contigs_ that will be used in the pipeline. To do this, we will collect the [SRA _reads_](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR492065&display=download) from NCBI and then use [metaSPADes](https://github.com/ablab/spades) to assemble.

To download the _reads_ from SRA, we used the `fastq-dump` program from [__SRA ToolKit__](https://github.com/ncbi/sra-tools). The download was done using the command below, changing only the identifier of the reads in the SRA:

```console
$ fastq-dump --split-files --gzip [IDENTIFICADOR]
```
The R1 and R2 _reads_ (_forward_ and _reverse_) are then downloaded as separate files for the next steps. These _reads_ have already been processed (adapters removed and quality filtered).

The _contigs_ were assembled using the program [MetaSPADes (v3.15.5)](https://github.com/ablab/spades). The program was installed as shown in the tool's repository.

The command below was used to assemble the _contigs_:
```console
$ python3 spades.py --meta -k 21,33,55,77,99,113,117,121,127 -1 [READ_R1.fastq.gz] -2 [READ_R2.fastq.gz] -t 12 -m 12 -o [OUTPUT]
```
All the _contigs_ assembled were grouped into a single fasta file (multifasta) to be used in the pipeline.

###  Running the pipeline

O arquivo multifasta e todos os _reads_ coletados foram utilizados como input para a pipeline com o seguinte comando:

```console
$ (env) metagen-phage -file multifasta.fasta -reads reads/ -t 12
```

Optei por deixar a busca por profagos desligadas (-pro option) para reduzir o tempo de processamento. A execução foi feita em um desktop com __16 GB de memória RAM e um processador AMD Ryzen 5 4600H__. O tempo total utilizado para executar a pipeline foi de __cerca de 5h38min__.

At the end of the run, __34 possible phage MAGs__ were retrieved. Analyzing these MAGs with the CheckV tool, we found that 15 were of high quality, indicating a high degree of genome completeness. Only 1 was identified as a possible prophage. The results found with the CheckV tool can be seen in the 'quality_summary.tsv' file deposited in the 'analysis' directory.

---
