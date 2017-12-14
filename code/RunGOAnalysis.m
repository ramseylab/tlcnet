function clustGO=RunGOAnalysis(GO)

load 'AllGenes' 'masterInds';
geneListSuper = masterInds;
clear 'masterInds';

load MasterAnnotations 'psEntrezIDs';
load DiffExpGenes 'masterInds' 'psGeneSymbols';

geneListSub = masterInds;

if nargin == 0
  GO = geneont('File', 'gene_ontology.obo');
end

load ClusterAnalysis;

K = size(bestKC, 1);

threshLPV = -log10(0.1);

diffExpGO = AnalyzeGO(geneListSuper, ...
		      geneListSub, ...
		      threshLPV, ...
 		      GO);

clustGO = cell(1,K);
for k=1:K
  k
  clustInd = find(bestKidx == k);
  psGeneSymbols{clustInd}
  threshLPVClust = threshLPV;
  clustGO{k} = AnalyzeGO(geneListSuper, ...
			 geneListSub(clustInd), ...
			 threshLPVClust, ...
			 GO);
end

save 'GOAnalysis' 'diffExpGO' 'clustGO';


