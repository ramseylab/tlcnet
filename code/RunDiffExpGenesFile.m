function RunDiffExpGenesFile

load 'DiffExpGenes';
load 'GeneOntology';
load 'ClusterAnalysis';

N=size(means, 1);
SaveAnnotatedGeneFile([1:N], ...
		  psGeneSymbols, ...
		  psTitles, ...
		  psEntrezIDs, ...
		  psNames, ...
		  'DiffExpGenes.tsv', ...
		  psGOProcess, ...
		  psGOComponent, ...
		  psGOFunction, ...
		  goTerms, ...
		  means, ...
		  stdDevs, ...
		  uniqueExptNames, ...
		  numReplicates, ...
		  [], ...
		  {}, ...
		  bestKidx);

		  