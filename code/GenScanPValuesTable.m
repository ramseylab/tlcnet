function GenScanPValuesTable

'loading clusteranalysis'
load 'ClusterAnalysis' 'K';

'loading ClusterSizes'
clusterSizes = load('ClusterSizes.tsv');

'loading ExpressedGeneHits'
refHits = load('ExpressedGeneHits.tsv');

'loading ClusterHits'
data = load('ClusterHits.tsv');

T = size(data,1);

%[status, expNum] = system(sprintf('grep ''>'' ../PromScan/ExpGenes.fa | wc -l'));
%expNum=str2num(expNum)
expNum = load('ExpGenesNumScanned.txt')

pvs = zeros(T, K);

apv = zeros(1,T);

for m=1:T
  for c=1:K
    fref = refHits(m)/expNum;
    fclu = data(m,c)/clusterSizes(c);
	
      pvs(m,c) = BetterHygeccdf(data(m,c), ...
				expNum, ...
				refHits(m), ...
				clusterSizes(c));
  end

  apv(m) = BetterHygeccdf(sum(data(m,:)), ...
			  expNum, ...
			  refHits(m), ...
			  sum(clusterSizes));
end

save '../PromScan/ScanPValues2.tsv' 'pvs' '-ascii' '-tabs';
return;

