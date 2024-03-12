# Metagen phage pipeline

O objetivo desse repositório é o de demonstrar como uma pipeline para a recuperação de Bacteriófagos a partir de amostras de metagenomas pode ser feita.

Essa pipeline foi construída utilizando Python3 para a construção de uma interface de linha de comando (CLI) e a Workflow Description Language (WDL) para a construção da pipeline.

## Descrição da Pipeline

Dois _inputs_ são necessários para o funcionamento adequado da pipeline:

- __Reads:__ _Reads_ dos _contigs_ de interesse para a busca de bacteriófagos. Esses devem ser _reads_ _paired-end_ e nomeados com _1.fastq e _2.fastq para identificação de correta;
- __Multifasta:__ Arquivo _multifasta_ contendo todos os contigs de interesse e que foram montados a partir dos __reads__ fornecidos.

As etapas abaixo são desenvolvidas ao longo da pipeline.

### Deduplicação

A deduplicação é a etapa de remoção de _contigs_ iguais construídos durante a montagem. Para essa etapa foi utilizado a ferramenta [dedupe do bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/). Essa ferramenta remove todos os _contigs_ que são idênticos (considerando também o reverso complementar).

### Mapeamento

Após a deduplicação, entramos na etapa de mapeamento dos _contigs_ contra os _reads_. Essa etapa é necessária para que seja calculada a cobertura de cada um dos _contigs_. Essa informação é importante para o processo de _binagem_ (Assista [Basics of metagenomic](https://www.youtube.com/watch?v=MqD4aN1p1qA)  para uma boa explicação!). Esse processo foi feito usando o [BWA](https://bio-bwa.sourceforge.net/).

### Sorting

Nessa etapa os arquivos SAM criados na etapa de mapeamento foram transformados em arquivos BAM (forma binária do SAM) e, então, foi feito o _sorting_ dos arquivos BAM. Todo esse processo foi feito utilizando as ferramentas do [samtools](http://www.htslib.org/).

### Binning

Nessa etapa iremos fazer a montagem dos _bins_, também chamados de _Metagenome Assembled Genomes_ (MAGs). Essa montagem foi feita por meio da ferramenta [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) com alguns parâmetros alterados para melhorar a _binagem_ de MAGs virais (-m 1500 -s 10000).

### Predição

A etapa final da pipeline é a predição de bacteriófagos nos MAGs montados. Para essa tarefa foi utilizada a ferramenta de predição [VirSorter2](https://github.com/jiarong/VirSorter2). Antes de realizar a predição a pipeline faz a concatenação dos contigs dos MAGs, colocando Ns entre eles. Esse MAG tem o seu _header_ renomeado para o nome do arquivo e, então, é submetido para predição pela ferramenta. Dado que o interesse dessa pipeline é a recuperação de MAGs de fagos, foi utilizado apenas o grupo de predição `dsDNAphage` do VirSorter2.

Após a predição, os resultados podem ser analisados e interpretados para identificar e caracterizar os bacteriófagos presentes nos metagenomas. Essa informação é valiosa para entender a ecologia microbiana e a interação entre bactérias e bacteriófagos em diferentes ambientes.

## Instalação

## Estudo de Caso

Buscando demontrar a eficiência dessa pipeline, foi feito um estudo de caso utilizando dados reais. Esses dados provêm do metagenoma de recém nascidos e foram explorados no artigo de [SHARON, Itai et al](https://pubmed.ncbi.nlm.nih.gov/22936250/).

### Recuperação dos _reads_ e montagem dos _contigs_

O passo inicial desse projeto foi a montagem das _reads_ para que consigamos os _contigs_ que serão usados na pipeline. Para tal iremos coletar as [_reads_ do SRA](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR492065&display=download) do NCBI e então usar o [metaSPADes](https://github.com/ablab/spades) para poder fazer a montagem.

Para realizar o download das _reads_ a partir do SRA foi utilizado o programa `fastq-dump` do [__SRA ToolKit__](https://github.com/ncbi/sra-tools). O download foi feito usando o comando abaixo, alterando-se apenas o identificador dos reads no SRA:

```console
$ fastq-dump --split-files --gzip [IDENTIFICADOR]
```

Dessa forma é feito o download das _reads_ R1 e R2 (_forward_ e _reverse_) já em arquivos separados para as próximas etapas. Essas _reads_ já tinham sido processadas (tiveram adaptadores removidos e qualidade avaliada).

A montagem dos _contigs_ foi feita utilizando o programa [MetaSPADes (v3.15.5)](https://github.com/ablab/spades). A instalação do programa foi feita conforme o apresentado no repositório da ferramenta.
O comando abaixo foi usado para fazer a montagem dos _contigs_:
```console
$ python3 spades.py --meta -k 21,33,55,77,99,113,117,121,127 -1 [READ_R1.fastq.gz] -2 [READ_R2.fastq.gz] -t 12 -m 12 -o [OUTPUT]
```

Todos os _contigs_ montados foram agrupados em um único arquivo fasta (multifasta) para ser usado na pipeline de recuperação de bacteriófagos.

