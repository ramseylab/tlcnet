function clust=GetClusterForGenes(symbols)

load 'DiffExpTFGenes';
load 'ClusterAnalysis' 'K' 'bestKC';

G = length(symbols);
M=size(means,2);
meanDiff = means - repmat(means(:,refCol),1,M);
normRelMeansMore = meanDiff ./ repmat(max(abs(meanDiff),[],2),1,M);

clust = zeros(1,G);
for g=1:G
  geneSymbol = symbols{g};
  psind = find(1==strcmp(psGeneSymbols, geneSymbol));
  if length(psind) > 0
    nrlgene = normRelMeansMore(psind, :);
    bestClustDiff = inf;
    bestClustNum = -1;
    for k=1:K
      theClustDiff=sum(abs(bestKC(k,:)-nrlgene));
      if theClustDiff < bestClustDiff
	bestClustDiff = theClustDiff;
	bestClustNum=k;
      end
    end
    clust(g)=bestClustNum;
  end
end

