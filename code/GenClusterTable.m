function GenClusterTable

load 'ClusterAnalysis';
load 'DiffExpGenes' 'psGeneSymbols';

cytochemos=textread('CytoChemokines.txt','%s');

[lpsInds,lpsTimes] = GetColsets('WT', 'LPS', 1, 2880);
M=size(normRelMeans,2);
normRelMeansMeds = zeros(K,M);

% create array to hold extreme change in each cluster
extChanges = [];
for k=1:K
  kind = find(k==bestKidx);
  normRelMeansMeds(k,:)=bestKC(k,:);
  extChanges(k) = Extreme(normRelMeansMeds(k,lpsInds));
end
% convert the extreme changes to plus one or minus one
% (just retaining the sign information about the values)
extChanges = sign(extChanges);
extChanges(33)=1;

ilpsTimes = [0 lpsTimes']
ilpsInds = [refCol lpsInds'];

indTimes = [];
for k=1:K
  thresh = 0.25;
  try
    indTimes(k)=GetInductionTime(ilpsTimes, ...
				 normRelMeansMeds(k,ilpsInds), ...
				 thresh);
  catch
    'unable to get induction time: ' 
    k
  end
end

load 'InductionTimes' 'indTimes';

load 'DiffExpGenesTFStatusStrict';

fid=fopen('ClusterTable.tsv','w+');
for k=1:K
  kind = find(bestKidx == k);
  [sortnames, kindsi]=sort(psGeneSymbols(kind));
  clustTFInds = kind(find(psIsTF(kind)));
  C = length(kind);
  cytochemosInCluster = '';
  for j=1:C
    geneName = psGeneSymbols{kind(kindsi(j))}
    if length(find(strcmp(geneName, cytochemos))) > 0
      if length(cytochemosInCluster) > 0
	cytochemosInCluster = [cytochemosInCluster ', '];
      end
      cytochemosInCluster = [cytochemosInCluster geneName];
    end
  end
  
  tfsInCluster = '';
  I = length(clustTFInds);
  [sctf, sctfi]=sort(psGeneSymbols(clustTFInds));
  for i=1:I
    tfName = psGeneSymbols{clustTFInds(sctfi(i))};
    tfNameInd = strfind(tfName, ' /// ');
    if strcmp(tfName, 'LOC669007')
      tfName = 'BAP135';
    end
    if length(tfNameInd) > 0
      tfName = tfName(1:(tfNameInd(1)-1));
    end
    tfsInCluster = [tfsInCluster tfName];
    if i < I
      tfsInCluster = [tfsInCluster ', '];
    end
  end
  
  fprintf(fid, ...
	  'C%d\t%d\t\t%d\t\t%s\t%s\n', ...
	  k, ...
	  C, ...
	  round(indTimes(k)), ...
	  tfsInCluster, ...
	  cytochemosInCluster);
end

fclose(fid);

