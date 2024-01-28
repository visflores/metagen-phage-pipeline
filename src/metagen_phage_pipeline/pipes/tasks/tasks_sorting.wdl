version 1.1

task ConvertToBam {
	input {
		File samFile
	}

	String outName = basename(samFile, ".aln-pe.sam")

	command <<<
		samtools view -b "~{samFile}" > "~{outName}.bam"
	>>>

	output {
		File bamFile = "~{outName}.bam"
	}

	runtime {
		container: "staphb/samtools:latest"
	}
}

task ExecuteSorting {
	input {
		File bamFile
	}

	String outName = basename(bamFile, ".bam")

	command <<<
		samtools sort "~{bamFile}" -o "~{outName}.sorted.bam"
	>>>

	output {
		File sortedBam = "~{outName}.sorted.bam"
	}

	runtime {
		container: "staphb/samtools:latest"
	}
}
