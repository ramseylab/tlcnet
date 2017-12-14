function RunGetDiffExpGenes

stims = {'LPS', 'PIC', 'PAM3', 'PAM2', 'CPG', 'R848', 'PAM3PIC'};
strains = {'WT'};

load 'NormalizedReplicates';
load 'MasterAnnotations';

indProbes = find(strncmp(psNames,'AFFX-',5)==0);

[subselectedIndices, ...
 means, ...
 stdDevs, ...
 pValues, ...
 pValueCutoffs, ...
 tcExptNames] = ...
    CombineDataAndSelectProbesets(normDat(indProbes,:), ...
				  psNames(indProbes), ...
				  psEntrezIDs(indProbes), ...
				  7, ...
				  exptNames, ...
				  arrayNames, ...
				  stims, ...
				  strains, ...
				  1000, ...
				  0.4, ...
				  50, ...
				  0.0001);

means=means(subselectedIndices,:);
stdDevs=stdDevs(subselectedIndices,:);
pValues=pValues(subselectedIndices,:);

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

save 'DiffExpGenes';

