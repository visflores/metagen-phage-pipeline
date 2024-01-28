#!/bin/bash

# ======================================================================================================================================
# ================================================== CABEÇALHO DE INFORMAÇÕES ==========================================================
# ======================================================================================================================================

# O programa desenhado nesse código é a automatização do pipeline disponível no repositório do MARVEL. Teremos como resultado final de
#sua execução a criação dos bins e sua triagem pelo MARVEL dizendo quais são fagos. Tenha em mente que é necessário que o sistema tenha
#instalado o MetaBAT2, SamToolse o BWA. O Marvel será baixado direto do repositório do github e usado para a triagem dos bins.

# Todos os resultados da execução desse programa irão ser depositados na pasta "execucao/", qual será criada por ele no diretório raiz
#de execução do algoritmo.

# Para que seja feita a execução do programa é necessário passar o seu arquivo de montagem (assembly.fa ou assembly.fasta) no formato
#fasta e a pasta contendo os reads de sequenciamento. Esse pipeline pode funcionar tanto com paired reads, single reads e uma combina
#ção de ambos.

# ====================================================================================================================================
# ==================================================== CHECANDO PARAMETROS ===========================================================
# ====================================================================================================================================

# Vamos checar se os parametros necessários foram passados

# $1 --> assembly.fasta $2 --> reads/ $3 --> Threads $4 --> RAM $5 --> Seed $6 --> tipo de execução

# O parametro $6 indica qual será o tipo de execução: paired, single ou ambos.

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ] | [ -z $5 ] | [ -z $6 ]
	then
		printf "Estão faltando parametros para serem passados!\n"
		exit -1
fi

# =====================================================================================================================================
# ========================================================= VALIDAÇÕES ================================================================
# =====================================================================================================================================

# Seção do código dedicada a validação das informações que devem ser passadas para o pipeline.

# Verificar se os arquivos estão no formato correto.

# Irei ver isso com base nas extensões dos arquivos.

# Verificando se o arquivo de montagem tem a extensão .fa ou .fasta

ls -1 $1 > $1.tmp

assembly=$(sed 's/\./,/' $1.tmp| grep -e 'fasta' -e 'fa')

if [ -z $assembly ]
	then
		printf "O arquivo $1 precisa possuir a extensão .fasta ou .fa \n\n"
		printf "Abortando execução do pipeline!"
		rm $1.tmp
		exit 1
	else
		printf "Tudo certo com o arquivo $1 \n\n"
fi

rm $1.tmp

# Vamos verificar agora se os reads estão no formato correto. Para essa etapa o usuário deve passar uma pasta contendo apenas os reads.
#Iremos verificar se eles possuem a extensão fastq e quais são os reads. Caso a execução seja por ambos, devemos deixar todos os reads
#no mesmo diretório. Leia a seção MAPEAMENTO para revisar a regra de nomenclatura.

if [ -d $2 ]
	then
		printf "O diretório $2 apresenta `ls -1 $2 | wc -l` arquivos. \n\n"
		
		printf "Avaliando os arquivos... \n\n"

		ls -1 $2 > $2.tmp

		grep -r ".*\.fastq$" $2.tmp > $2.tmp.teste

		if [ -s $2.tmp.teste ]
			then
				printf "Os arquivos estão no formato fastq! \n\n"
			else
				printf "Os reads não estão no formato fastq. Por favor verificar seus reads."
				exit 1
		fi

		rm $2.tmp.teste
	else
		printf "Por favor, submeta um diretório contendo os reads"
		exit 1
fi

# =======================================================================================================================================
# ================================================ DEDUPLICAÇÃO E MAPEAMENTO ============================================================
# =======================================================================================================================================


# Após avaliarmos se os arquivos estão no formato correto, vamos fazer a deduplicação do fasta com os contigs e em seguida realizar o
#mapeamento usando o BWA.

# ---
# Retirando os contigs duplicados da montagem. Vamos colocar o novo arquivo da deuplicação no diretório "deduped_assembly/" dentro da
# pasta "execucao".

# Verificando se a pasta execucao foi criada. Caso não, ele irá a criar. Caso sim, ela irá ser apagada e novamente criada.

if [ -d "./execucao" ]
	then
		rm -r execucao/
		mkdir execucao
		mkdir -p execucao/deduped_assembly execucao/mapeamento execucao/metabat execucao/phage_bins execucao/MARVEL
	else
		mkdir execucao
		mkdir -p execucao/deduped_assembly execucao/mapeamento execucao/metabat execucao/phage_bins execucao/MARVEL
fi

# Retirando duplicações com o dedupe do Samtools

dedupe in=$1 out=execucao/deduped_assembly/deduped_assembly.fasta

if [ ${?} -eq 0 ]
	then
		printf "O arquivo foi deduplicado com sucesso!! \n\n"
	else
		printf "Ocorreu um problema durante a execução do dedupe. \n\n"
		printf "Abortando execução do pipeline!!!"
		exit -1
fi


# Após a deduplicação, iremos realizar o mapeamento dos contigs diante dos reads. Inicialmente temos que preparar o ambiente
#com o comando BWA index no diretório em que o arquivo deduplicado da montagem irá ser depositado.

printf "Preparando o ambiente para o mapeamento!!\n\n"

bwa index execucao/deduped_assembly/deduped_assembly.fasta

if [ ${?} -eq 0 ]
	then
		printf "Ambiente preparado... Vamos iniciar o mapeamento!\n\n"
		printf "Esse processo é demorado, tenha paciência!!\n\n"
	else
		printf "Houve algum problema com a preparação do ambiente\n"
		printf "Abortando o pipeline!"
		exit -1
fi

# REGRAS DE NOMENCLATURA E DISPOSIÇÃO DOS READS
# PAIRED READS
# Caso opte por executar o pipeline apenas com reads paired-ends, você deve nomear os fastq da seguinte forma:
	# readR_1.fastq e readR_2.fastq
# A particula R_[num].fastq será importante para que os reads sejam encontrados e mapeados corretamente. É necessário que  o nome
#que vier antes da particula R_[num].fastq seja igual em ambos os reads, por exemplo:
	# ZC4DAY01R_1.fastq e ZC4DAY01R_2.fastq --> Assim não ocorrerá erro durante a execução.

# SINGLE READS
# Caso opte por executar o pipeline com reads single-end, tenham em mente apenas que os reads devem possuir a extensão fastq. Dessa
#forma não ocorrerá problemas com a validação inicial. Aqui nenhuma nomenclatura específica no read é necessária, já que o pipeline
#irá percorrer a pasta fastq/ e executar o BWA em single reads em cada um dos arquivos dela.

# AMBOS OS TIPOS
# Aqui a nomenclatura utilizada para paired-reads também é necessária, portanto deve ser seguida. Agora temos também uma nomenclatura
#específica para os reads single-end. Devemos nomear esses reads com a seguinte regra:
	#readS.fastq
# A particula S.fastq será usada para diferenciar quais são os reads que devem ser alinhados com o programa paired-end e quais devem
#ser alinhados com o programa single-end do BWA MEM. Respeitando essas regras, o programa irá executar sem nenhum problema.

# Antes de iniciarmos o mapeamento, tenho que montar os pares dos readsR1 e R2. Para isso irei fazer um match com a porção dos nomes
#dos reads que devem ser iguais, caso sejam irei armazenar esse par em um arquivo e executar o BWA.

# Esse bloco de código executa o pipeline apenas para single-reads
if [ $6 == 'single' ]
	then
		printf "Executando pipeline para Single-reads!\n"

		grep -r ".*\.fastq" $2.tmp > $2.tmp.single

		sed 's/.fastq//' $2.tmp.single > $2.tmp.single.edited

		for name in `cat $2.tmp.single.edited`
			do
				printf "Mapeando o read $name!\n"	

				bwa mem -t $3 execucao/deduped_assembly/deduped_assembly.fasta $2$name.fastq > execucao/mapeamento/$name.aln-se.sam
				if [ ${?} -eq 0 ]
					then
						printf "Mapeando....\n\n"
					else
						rm $2.tmp*
					
						printf "Ocorreu um problema durante o mapeamento!\n Abortando o Pipeline!\n"
						
						exit -1

				fi
				 
			done
		
		printf "Mapeamento finalizado!!\n"

		rm $2.tmp
# Bloco para executar o pipeline para apenas Paired-reads
elif [ $6 == 'paired' ]
	then
		printf "Executando pipeline para Paired-reads!\n"

		# Separando os reads R1 e R2 em arquivos separados

		grep -r ".*R_1\.fastq" $2.tmp > $2.tmp.r1

		grep -r ".*R_2\.fastq" $2.tmp > $2.tmp.r2

		# Vamos agora combinar os reads R1 e R2

		sed 's/R_1.fastq//' $2.tmp.r1 > $2.tmp.r1.edited

		for name in `cat $2.tmp.r1.edited`
			do
				search=$(grep -w "$name\R_2.fastq" $2.tmp.r2)

				if [ -z $search ]
					then
						printf "Não foi possível econtrar o read R2 de $name.\n"
				
						printf "Abortando o pipeline!\n"
				
						rm $2.tmp*

						exit -1
				else
					erreum=$(grep -w "$name\R_1.fastq" $2.tmp.r1)
					erredois=$(grep -w "$name\R_2.fastq" $2.tmp.r2)
				
					printf "Irei mapear os reads $erreum e $erredois agora!!\n\n"

					bwa mem -t $3 execucao/deduped_assembly/deduped_assembly.fasta $2$erreum $2$erredois > execucao/mapeamento/$name.aln-se.sam

					if [ ${?} -eq 0 ]
						then
							printf "Mapeando....\n\n"
						else
							rm $2.tmp*
					
							printf "Ocorreu um problema durante o mapeamento!\n Abortando o Pipeline!\n"
						
							exit -1

					fi

				fi

		done

		printf "Mapeamento finalizado!!\n"

		rm $2.tmp*
# Bloco para executar o pipeline tanto para Single quanto Paired-reads
elif [ $6 == 'ambos' ]
	then
		printf "Executando pipeline para Single-reads e Paired-reads!\n"
		
		printf "Vamos mapear os Single-reads primeiro!\n"

		grep -r '.*S\.fastq' $2.tmp > $2.tmp.single
		
		sed 's/S.fastq//' $2.tmp.single > $2.tmp.single.edited
		
		counter=1

		for name in `cat $2.tmp.single.edited`
			do
				printf "Mapenado o read $name!\n"

				bwa mem -t $3 execucao/deduped_assembly/deduped_assembly.fasta $2$name\S.fastq > execucao/mapeamento/single_$counter\_$name.aln-se.sam
				if [ ${?} -eq 0 ]
					then
						printf "Mapeando Single-reads....\n\n"
					else
						rm $2.tmp*
					
						printf "Ocorreu um problema durante o mapeamento dos Single-reads!\n Abortando o Pipeline!\n"
						
						exit -1

				fi
				
				((counter=counter+1))
			done

		printf "Vamos agora mapear os Paired-reads!\n"

		# Separando os reads R1 e R2 em arquivos separados

		grep -r ".*R_1\.fastq" $2.tmp > $2.tmp.r1

		grep -r ".*R_2\.fastq" $2.tmp > $2.tmp.r2

		# Vamos agora combinar os reads R1 e R2

		sed 's/R_1.fastq//' $2.tmp.r1 > $2.tmp.r1.edited

		for name in `cat $2.tmp.r1.edited`
			do
				search=$(grep -w "$name\R_2.fastq" $2.tmp.r2)

				if [ -z $search ]
					then
						printf "Não foi possível econtrar o read R2 de $name.\n"
				
						printf "Abortando o pipeline!\n"
				
						rm $2.tmp*

						exit -1
				else
					erreum=$(grep -w "$name\R_1.fastq" $2.tmp.r1)
					erredois=$(grep -w "$name\R_2.fastq" $2.tmp.r2)
				
					printf "Irei mapear os reads $erreum e $erredois agora!!\n\n"

					bwa mem -t $3 execucao/deduped_assembly/deduped_assembly.fasta $2$erreum $2$erredois > execucao/mapeamento/$name.aln-se.sam

					if [ ${?} -eq 0 ]
						then
							printf "Mapeando Paired-reads....\n\n"
						else
							rm $2.tmp*
					
							printf "Ocorreu um problema durante o mapeamento!\n Abortando o Pipeline!\n"
						
							exit -1

					fi

				fi

		done

		printf "Mapeamento finalizado!!\n"

		rm $2.tmp
else
	printf "\nAlgum problema ocorreu durante o mapeamento, abortanto pipeline!!\n"

	exit -1
fi

# ======================================================================================================================================
# ==================================================  SAMs -> BAMs -> Sorted BAMs ======================================================
# ======================================================================================================================================

# Criando estrutura de diretórios para armazenar os BAMs e BAMs sorteados.

if [ -d "./execucao/mapeamento/BAM" ] || [ -d "./execucao/mapeamento/BAM/SORTED" ]
	then
		printf "Removendo diretórios BAMs antigos. \n\n"
		rm -r execucao/mapeamento/BAM
		printf "Criando novo diretório BAM.\n\n"
		mkdir -p execucao/mapeamento/BAM execucao/mapeamento/BAM/SORTED
	else
		printf "Criando diretórios para armazenar BAMs e BAMs sorteados.\n\n"
		mkdir -p execucao/mapeamento/BAM execucao/mapeamento/BAM/SORTED
fi

# Transformando arquivos SAM em BAM e depois fazendo o sorting dos arquivos BAM.

ls -1 execucao/mapeamento/*aln-se.sam > sorting.tmp

sed -i "s#execucao/mapeamento/##" sorting.tmp
sed -i "s/\.aln-se\.sam//" sorting.tmp

# Loop para converter para BAM e depois sortear.
for sam in `cat sorting.tmp`
	do
		printf "Convertendo o arquivo $sam\n"
		
		samtools view -b execucao/mapeamento/$sam.aln-se.sam > execucao/mapeamento/BAM/$sam.raw.bam

		if [ ${?} -eq 0 ]
			then
				printf "Ok, tudo certo com a conversão para BAM. Realizando o Sorting!\n"

				samtools sort execucao/mapeamento/BAM/$sam.raw.bam -o execucao/mapeamento/BAM/SORTED/$sam.bam

				if [ ${?} -eq 0 ]
					then
						printf "O arquivo BAM $sam foi sorteado com sucesso!!\n\n"
					else
						printf "Houve algum problema ao realizar o sorteamento do arquivo BAM $sam.\n\n"
						printf "Abortando o pipeline!\n"
						exit -1
				fi

			else
				printf "Houve algum erro ao realizar a conversão do arquivo $sam.\n\n"
				printf "Abortando o pipeline!\n"
				exit -1
		fi

	done

# Eliminando arquivo temporário.
rm sorting.tmp

# ======================================================================================================================================
# =========================================================  MetaBAT2 ==================================================================
# ======================================================================================================================================

# Antes de rodar o MetaBAT2 em si, devemos construir um arquivo contendo as informações de cobertura para a montagem. Para isso iremos
#usar o programa jgi_summarize_bam_contig_depths.

printf "Criando arquivo depths.txt em execucao/metabat/\n\n"

jgi_summarize_bam_contig_depths --outputDepth execucao/metabat/depths.txt execucao/mapeamento/BAM/SORTED/*.bam

if [ ${?} -eq 0 ]
	then
		printf "Arquivo depths.txt criado com sucesso!\n\n"
	else
		printf "Houve algum problema ao rodar o programa jgi_summarize_bam_contig_depths.\n\n"
		printf "Abortando o pipeline\n"
		exit -1
fi

# Vamos agora rodar o MetaBat2. Para que ele seja rodado é necessário que um seed seja passado para fazer a execução. Caso queira que
#esse seed seja aleatório (escolhido pelo MetaBAT2), coloque o valor 0.

printf "Construíndo os bins com o MetaBAT2!\nTodos os bins serão depositados no diretório execucao/metabat/bins/\n\n"

mkdir -p execucao/metabat/bins

metabat2 -i execucao/deduped_assembly/deduped_assembly.fasta -a execucao/metabat/depths.txt -o execucao/metabat/bins/bin -m 1500 -s 10000 --seed $5 -t $3

if [ ${?} -eq 0 ]
	then
		printf "Os bins foram criados com sucesso!!\n\n"
		printf "Foram construídos `ls -1 execucao/metabat/bins/bin.[0-9]*.fa | wc -l`.\n"
		#printf "Contigs não binados --> `grep -c '^>' execucao/metabat/bins/bin.unbinned*`.\n\n"

	else
		printf "Algum problema ocorreu com a execução do MetaBAT2!\n\n"
		printf "Abortando o pipeline\n"
		exit -1
fi

# ======================================================================================================================================
# =============================================================  MARVEL ================================================================
# ======================================================================================================================================

# Abaixo iremos baixar, treinar o modelo e rodar o MARVEL para buscar os bins que são fagos e foram produzidos pelo MetaBAT2.

# Fazendo o Download do MARVEL

git clone https://github.com/LaboratorioBioinformatica/MARVEL execucao/MARVEL

# Verificando se houve algum problema com a execução do comando anterior.
if [ ${?} -eq 0 ]
	then
		printf "Download do MARVEL feito com sucesso!\n\n"
	else
		printf "Houve algum problema ao baixar o MARVEL. Por favor verificar sua conexão.\n\n"
		printf "Abortando o pipeline\n"
		exit -1
fi

# Carregando os modelos do MARVEl. Realizando o treinamento do algoritmo de Random Forest.
# Aqui é necessário rodar um subprocesso para executar o script. Ele apenas roda do diretório corrente.
(cd execucao/MARVEL/ && python3 download_and_set_models.py)


# Verificando se houve algum problema com a excução do comando anterior.
if [ ${?} -eq 0 ]
	then
		printf "MARVEL treinado com sucesso, vamos agora executar o script principal!\n\n"
	else
		printf "Houve algum problema ao treinar o modelo.\n\n"
		printf "Abortando pipeline!\n"
		exit -1
fi

(cd execucao/MARVEL/ && python3 marvel_bins.py -i ../metabat/bins -t $3 > marvel.out)

if [ ${?} -eq 0 ]
	then
		mv execucao/metabat/bins/results execucao/phage_bins
		printf "MARVEL executado com sucesso!!\n\n"
		printf "O MARVEL predisse `ls -1 execucao/phage_bins/results/phage_bins/ | wc -l` bins como sendo fagos!\n\n"
		printf "Finalizando o Pipeline!\n"
	else
		printf "Houve algum problema ao rodar o MARVEL.\n\n"
		printf "Abortando o Pipeline!"
		exit -1
fi

