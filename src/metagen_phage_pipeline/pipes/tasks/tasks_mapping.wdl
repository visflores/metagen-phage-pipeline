version 1.1

task ExecuteMapping {
	input {
		File dedupedFasta
		Pair[File, File] reads
	}

	String outName = basename(reads.left, ".fastq")

	command <<< 
		bwa index "~{dedupedFasta}" 
		bwa mem -t 2 "~{dedupedFasta}" "~{reads.left}" "~{reads.right}" > "~{outName}.aln-pe.sam" 
	>>>

	output {
		File mappingFile = "~{outName}.aln-pe.sam"
	}

	runtime {
		container: "staphb/bwa:latest"
	}
}

