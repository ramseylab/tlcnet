function RunGetDiffExpTFGenes

stims = {'LPS', 'PIC', 'PAM3', 'PAM2', 'CPG', 'R848', 'PAM3PIC'};
strains = {'WT'};

load 'NormalizedReplicates';
load 'MasterAnnotations';

indProbes = find(strncmp(psNames,'AFFX-',5)==0);
geneSymbolsRealProbes = psGeneSymbols(indProbes);

N = length(indProbes);
geneSymbolsRealProbesShort = cell(1,N);
for i=1:N
  geneSymbol = geneSymbolsRealProbes{i};
  geneSymbolsRealProbesShort{i} = geneSymbol;
  k = strfind(geneSymbol, ' /// ');
  if length(k) > 0
    geneSymbolsRealProbesShort{i} = geneSymbol(1:(k(1)-1));
  end
end

[MatName, VMatName, TNum, tfEntrezGeneIDs, tfGeneSymbols]= ...
    textread('ThorssonTFsWithMatrices.tsv', '%s %s %s %d %s');

tfEntrezGeneIDs = unique(tfEntrezGeneIDs);
G = length(tfEntrezGeneIDs);

psEntrezGeneIDsRealProbes = psEntrezIDs(indProbes);
indProbesTFs = [];
P = length(indProbes);
for g=1:G
  tfEntrezGeneID = tfEntrezGeneIDs(g);
  for p=1:P
    if find(psEntrezGeneIDsRealProbes{p}==tfEntrezGeneID)
      indProbesTFs = [indProbesTFs indProbes(p)];
    end
  end
end

tfGeneSymbols = unique(tfGeneSymbols);
G = length(tfGeneSymbols);
for g=1:G
  tfGeneSymbol = tfGeneSymbols{g};
  ind = find(strcmp(tfGeneSymbol, geneSymbolsRealProbesShort));
  if length(ind > 0)
    indProbesTFs = [indProbesTFs indProbes(ind)];
  end
end

indProbesTFs = unique(indProbesTFs);

[subselectedIndices, ...
 means, ...
 stdDevs] = ...
    CombineDataAndSelectProbesets(normDat(indProbesTFs,:), ...
				  psNames(indProbesTFs), ...
				  psEntrezIDs(indProbesTFs), ...
				  7, ...
				  exptNames);

means=means(subselectedIndices,:);
stdDevs=stdDevs(subselectedIndices,:);
masterInds = indProbesTFs(subselectedIndices);
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

clear 'normDat' 'MatName' 'VMatName' 'P' 'G' 'g' 'TNum' ...
    'geneSymbolsRealProbes' 'arrayDataOnly' 'tfGeneSymbols' ...
    'tfEntrezGeneIDs' 'tfEntrezGeneID' 'tfGeneSymbol' 'tfGeneSymbols' ...
    'p' 'ind' 'subselectedIndices' 'indProbes' 'psEntrezGeneIDsRealProbes' ...
    'indProbesTFs';

uniqueExptNames = unique(exptNames);

save 'ExpTFGenes';


