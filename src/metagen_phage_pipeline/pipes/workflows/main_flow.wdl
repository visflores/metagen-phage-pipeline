version 1.1

import "../tasks/tasks_dedup.wdl" as DedupTasks
import "../tasks/tasks_mapping.wdl" as MappingTasks
import "../tasks/tasks_sorting.wdl" as SortingTasks
import "../tasks/tasks_binning.wdl" as BinningTasks

workflow MainFlow {
	input {
		File multiFasta
		Array[Pair[File, File]] readsToMap
	}

	call DedupTasks.ExecuteDedupe {
		input:
			multiFasta = multiFasta
	}

	scatter (reads in readsToMap) {
		call MappingTasks.ExecuteMapping {
			input:
				dedupedFasta = ExecuteDedupe.dedupedFasta,
				reads = reads
		}
	}

	scatter (sam in ExecuteMapping.mappingFile) {
		call SortingTasks.ConvertToBam {
			input:
				samFile = sam
		}

		call SortingTasks.ExecuteSorting {
			input:
				bamFile = ConvertToBam.bamFile
		}
	}

	call BinningTasks.SummarizeDepth {
		input:
			sortedBams = ExecuteSorting.sortedBam
	}

	call BinningTasks.ExecuteBinning {
		input:
			depths = SummarizeDepth.depths,
			dedupedFasta = ExecuteDedupe.dedupedFasta 
	}	 	
}
