version 1.1

task ExecuteDedupe {
	input {
		File multiFasta
	}

	command <<< 
		dedupe.sh in="~{multiFasta}" out=deduped_multifasta.fasta
	>>>

	output {
		File dedupedFasta = "deduped_multifasta.fasta"
	}

	runtime {
		container: "staphb/bbtools:latest"
	}
}
