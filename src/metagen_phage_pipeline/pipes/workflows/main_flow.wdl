version 1.1

import "../tasks/tasks_dedup.wdl" as DedupTasks
import "../tasks/tasks_mapping.wdl" as MappingTasks
import "../tasks/tasks_sorting.wdl" as SortingTasks
import "../tasks/tasks_binning.wdl" as BinningTasks
import "../tasks/tasks_prediction.wdl" as PredictionTasks

workflow MainFlow {
	input {
		File multiFasta
		Array[Pair[File, File]] readsToMap
		Int threads = 8
		Boolean findProvirus = true
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

	scatter (bin in ExecuteBinning.bined) {
		call PredictionTasks.AssembleInput {
			input:
				bin = bin
		}

		call PredictionTasks.ExecutePrediction {
			input:
				bin = AssembleInput.concatBins,
				threads = threads,
				findProvirus = findProvirus
		}
	}

	call PredictionTasks.MergeResults {
		input:
			phageBins = ExecutePrediction.phageBins
	}

	call PredictionTasks.PhageBins {
		input:
			phageBins = MergeResults.phageBinsMerged,
			bins = ExecuteBinning.bined,
			errorCheck = ExecutePrediction.errorCheck
	}

	output {
		File phageInfos = MergeResults.phageBinsMerged
		Array[File] phageMags = PhageBins.phageMags
		Array[File] noPhageMags = PhageBins.noPhageMags
		File errorFile = PhageBins.errorFile
	}
}
