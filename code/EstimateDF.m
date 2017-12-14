function EstimateDF

load ClusterAnalysis
load TimeLaggedCorr

tlcColsAll = [];
L = length(tlcCols);
for l=1:L
  tlcColsAll = [tlcColsAll tlcCols{l}'];
end
tlcColsAll = unique(tlcColsAll);

dfk = zeros(1,K);
for k=1:K
  kind = find(k==bestKidx);
  kexp = normRelMeans(kind, tlcColsAll);
  
  [bk, bkc, kv, kbic]=ClustBIC(kexp, ...
			       [4:min(max(round(0.2*length(kind)),4),20)], ...
			       100, ...
			       5);
  dfk(k) = size(bkc, 1);
end

dfk

clustSizes = zeros(1,32);
for k=1:K
  clustSizes(k) = length(find(k==bestKidx));
end

save 'EstimateDF' 'dfk' 'clustSizes';
