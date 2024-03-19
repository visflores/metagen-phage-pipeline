version 1.1

task SummarizeDepth {
	input {
		Array[File] sortedBams
	}

	command <<< 
		jgi_summarize_bam_contig_depths --outputDepth depths.txt "~{sep(' ', sortedBams)}"
	>>>

	output {
		File depths = "depths.txt"
	}

	runtime {
		container: "metabat/metabat:latest"
	}
}

task ExecuteBinning {
	input {
		File depths
		File dedupedFasta
		Int seed = 1
	}

	command <<<
		metabat2 -i "~{dedupedFasta}" -a "~{depths}" -o bin \
			-m 1500 -s 10000 --seed "~{seed}" --unbinned -t 4
	>>>

	output {
		Array[File] bined = glob("*[0-9]*.fa")
		File unbinned = "bin.unbinned.fa"
	}

	runtime {
		container: "metabat/metabat:latest"
	}
}
