function GetGenesForCluster(cluster, fileName)

load 'ClusterAnalysis';
load 'DiffExpGenes';

ind = find(bestKidx == cluster);
L = length(ind);

if nargin > 1
  fid=fopen(fileName, 'w+');
else
  fid=1;
end
for l=1:L
  fprintf(fid, '%s\n', psGeneSymbols{ind(l)});
end
if nargin > 1
  fclose(fid);
end


