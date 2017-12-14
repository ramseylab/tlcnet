function ChIPChipGlobalAnalysis(tfList)

% get the number of clusters
load 'ClusterAnalysis' 'K' 'bestKidx';
load 'DiffExpGenes' 'psGeneSymbols';
load 'ClusterScans3' 'tfNames' 'pcsPV';
load 'TimeLaggedCorr2' 'tlcPVC' 'toptPVC' 'tlcShifts';

fid=fopen('ChIPChipGlobalAnalysis.tsv','w+');

for j=1:length(tfList)
  tfName = tfList{j};
  
  tfInd = find(strcmp(tfName, tfNames));
  if length(tfInd) == 0
    error 'could not find TF gene in list of scanned TFs';
  end

  tlcPVComb = tlcPVC .* toptPVC;

  combPV = gammainc(-log(tlcPVComb .* pcsPV), 2, 'upper')/gamma(2);

  % get list of all genes represented on the ChIP-chip array
  chipGenes = unique(upper(textread('../ChIPChip/NewListOfAllGenesTiled.txt', '%s')));

  % get the list of genes identified through ChIP-chip for the
  % specific TF target
  tfChIPFileName = sprintf('../ChIPChip/LPS_%s_ChIPChip_Close5.txt', ...
			   upper(tfName));
  chipTFGenes = unique(upper(textread(tfChIPFileName, '%s')));

  % loop through the clusters, one at a time
  for k=1:K
    % get the list of gene symbols for genes in this cluster, as a
    % cell array
    kind = find(k == bestKidx);
    clustGenes = psGeneSymbols(kind);
    
    % get the number of genes in the cluster
    G = length(kind);
    
    % for any gene name that contains the separator " /// ", just
    % take the first gene symbol before the separator
    for g=1:G
      geneName = clustGenes{g};
      gk = strfind(geneName, ' /// ');
      if length(gk) > 0
	clustGenes{g} = geneName(1:(gk(1) - 1));
      end
    end
    
    % upper-case all the gene symbols
    clustGenesUC = upper(clustGenes);
    
    % count the number of genes on the array
    M = length(chipGenes);
    
    % count the number of genes on the array that were found to have the TF bound
    Y = length(chipTFGenes);
    
    % verify that every gene listed in the TF-specific ChIP-chip data
    % file is also listed in the master list of genes for the
    % ChIP-chip array
    if Y ~= length(intersect(chipGenes, chipTFGenes))
      error 'there are genes in the TF-specific ChIP-Chip file that are not in the master gene list for the ChIP-Chip array';
    end
    
    % count the number of genes in the cluster that are represented
    % on the ChIP-Chip array
    N = length(intersect(clustGenesUC, chipGenes));
    
    % count the number of genes in the cluster that have the TF bound
    X = length(intersect(clustGenesUC, chipTFGenes));
    
    pv = 1 - hygecdf(X - 1, M, Y, N);
    
    inNet = '';
    
    meanShift = mean(tlcShifts(tfInd, kind));
    
    if (pcsPV(tfInd, k) < 0.05) && (combPV(tfInd, k) < 0.0248) && ...
	  (meanShift > 10)
      inNet = 'yes';
    end
    
    fracOnChip = N / length(kind);
    
    if fracOnChip >= 0.30
      
      line = sprintf('%s\t%d\t%d\t%d\t%d\t%0.2f\t%0.3e\t%0.3e\t%0.3e\t%0.1f\t%0.3e\t%s', ...
		   tfName, ...
		   k, Y, N, X, ...
		   fracOnChip, ...
		   pv, ...
		   pcsPV(tfInd, k), ...
		   tlcPVComb(tfInd, k), ...
		   meanShift, ...
		   combPV(tfInd, k), ...
		   inNet);
      
      disp(line)
      fprintf(fid, '%s\n', line);
    end
  end
end

fclose(fid);
