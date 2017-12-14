function RunGetExpGenes

load 'NormalizedReplicates';
load 'MasterAnnotations';

indProbes = find(strncmp(psNames,'AFFX-',5)==0);

[subselectedIndices, ...
 means, ...
 stdDevs] = ...
    CombineDataAndSelectProbesets(normDat(indProbes,:), ...
				  psNames(indProbes), ...
				  psEntrezIDs(indProbes), ...
				  7, ...
				  exptNames);

means=means(subselectedIndices,:);
stdDevs=stdDevs(subselectedIndices,:);

masterInds = indProbes(subselectedIndices);
psNames=psNames(masterInds);
psEntrezIDs=psEntrezIDs(masterInds);
psGeneSymbols=psGeneSymbols(masterInds);
psGOProcess=psGOProcess(masterInds);
psGOComponent=psGOComponent(masterInds);
psGOFunction=psGOFunction(masterInds);
psChromosome=psChromosome(masterInds);
psStartCoord=psStartCoord(masterInds);
psEndCoord=psEndCoord(masterInds);
psTitles=psTitles(masterInds);

length(subselectedIndices)

uniqueExptNames = unique(exptNames);

clear 'normDat' 'indProbes' 'subselectedIndices';

save 'ExpGenes';


