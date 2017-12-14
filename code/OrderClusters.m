function OrderClusters

load 'ClusterAnalysis';

inputMat = bestKC';

initDist = distance(inputMat)

S=simulatedannealing(inputMat, ...
		   40, ...
		   0.1, ...
		   50000, ...
		   3);

nonCycleDist = zeros(1,K);
Sshift=zeros(K,K);
for k=1:K
  for kp=1:K
    Sshift(k,kp)= rem(k + kp, K) + 1;
  end
  nonCycleDist(k) = distanceNonCyclic(inputMat(:,S(Sshift(k,:))));
end

[minDist, minDistInd] = min(nonCycleDist);

clusterOrder = S(Sshift(minDistInd,:));

figure;
imagesc(inputMat(:,clusterOrder));
load 'RedGreenColormap';
colormap(rgcMap);

save 'ClusterOrder' 'clusterOrder';
