version 1.1

task AssembleInput {
	input {
		File bin
	}

	String outName = basename(bin)

	command <<< 
		header=$(head -n 1 "~{bin}")
		N=$(printf "N%.0s" {1..50})
		sed "s/^>.*/$N/g" "~{bin}" > concat."~{outName}"
		sed -i "1i >~{outName}" concat."~{outName}"
	>>>

	output {
		File concatBins = "concat.~{outName}"
	}

	runtime {
		container: "ubuntu:latest"
	}
}

task ExecutePrediction {
	input {
		Array[File] bins
		Int threads
		Boolean findProvirus = false
	}

	command <<<
		cat ~{sep(' ', bins)} > tmp.fa
		if ~{findProvirus}; then
			virsorter run -w results -i tmp.fa -j ~{threads} --include-groups dsDNAphage
		else
			virsorter run -w results -i tmp.fa -j ~{threads} --include-groups dsDNAphage --provirus-off
		fi
	>>>

	output {
		File phageBins = "results/final-viral-score.tsv"
	}

	runtime {
		container: "vsflores/metagen-virsorter:latest"
	}
}

task PhageBins {
	input {
		File phageBins
		Array[File] bins
	}

	command <<<
		mkdir phage_mags no_phage_mags
		cut -f 1 -d $'\t' "~{phageBins}" > 1.tmp
		awk -F '|' '{print $1}' 1.tmp | sed '/seqname/d' | sort -u > 2.tmp
		for bin in ~{sep(' ', bins)}; do
			name=$(basename $bin)
			if [ $(grep -ic "$name" 2.tmp) -eq 1 ]; then
				cp $bin phage_mags/
			else
				cp $bin no_phage_mags/
			fi
		done
	>>>

	output {
		Array[File] phageMags = glob("phage_mags/*.fa")
		Array[File] noPhageMags = glob("no_phage_mags/*.fa")
	}

	runtime {
		container: "ubuntu:latest"
	}
}
