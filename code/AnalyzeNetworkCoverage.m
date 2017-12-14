function AnalyzeNetworkCoverage

% in this file, we want to go through the list of all 1,960
% differentially expressed genes, and for each gene, determine 
% if it possesses a motif hit for any TF for which a TF gene was 
% associated with the cluster with which the gene is a member

load 'DiffExpGenes' 'psGeneSymbols';

N = length(psGeneSymbols);

[tfgs, tfps, tggs, tgps, tgcs, corrs, lags, mscs, tlcpvs] = ...
    textread('TfLagToMotifTargets2.tsv', ...
	     '%s %s %s %s %d %f %d %f %f', ...
	     'delimiter', '\t');
H = length(tfgs);

[netcs, nettfs, nettfcs, netlpvs, netmats, netpctbs, nethits, ...
 netlpvexps, nettoptpvs, netmeanss, netlpvcs, netavgcs] = ...
    textread('Network2.tsv', ...
	     '%d %s %d %f %s %f %d %f %f %f %f %f', ...
	     'delimiter', '\t', ...
	     'headerlines', 1);

hitTargets = zeros(1,N);
E = length(netcs);
for i=1:E
  tfGene = nettfs{i};
  clust = netcs(i);
  
  foundInds = find(strcmp(tfgs, tfGene) & (tgcs==clust));
  targetGenes = tggs(foundInds);
  T = length(targetGenes);
  for j=1:T
    hitTargets(find(strcmp(psGeneSymbols, targetGenes{j})))=1;
  end
end

length(find(hitTargets))
