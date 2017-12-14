function RunClusterAnalysis

load 'DiffExpGenes' 'means' 'refCol';

[N,M]=size(means);

absInt = means;

refMeans = repmat(absInt(:,refCol),1,M);

relMeans = absInt - refMeans;

absRelMeans = abs(relMeans);

maxAbsRelMeans = repmat(max(absRelMeans, [], 2), ...
			1, M);

normRelMeans = relMeans ./ maxAbsRelMeans;

minK = 18;
maxK = 50;
varK = 3;
numRep = 1000;
krange = [minK:1:maxK];
[bestKidx, bestKC, kvals,kBIC]=ClustBIC(normRelMeans, krange, numRep, varK);
max(bestKidx)
K=size(bestKC,1);
clear 'krange' 'absRelMeans' 'relMeans' 'refMeans' 'maxAbsRelMeans' ...
    'means' 'absInt';
save 'ClusterAnalysis';



