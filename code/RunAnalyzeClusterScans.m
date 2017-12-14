function RunAnalyzeClusterScans2

clusterHits = load('ClusterHits3.tsv');
clusterSizes = load('../PromScan/ClusterSizes2.tsv');
expGeneHits = load('ExpressedGeneHits3.tsv');
tfMatrices = textread('../PromScan/ScannedMatrices.tsv', '%s');
expNum = textread('../PromScan/ExpGenesNumScanned.txt','%d')

% get number of matrices
Q = length(tfMatrices);

% get number of clusters
K = length(clusterSizes);

%[status, expNum] = system('grep ''>'' ../PromScan/ExpGenes.fa | wc -l');
%expNum=str2num(expNum);

clusterMatrixPVs = ones(Q,K);
for m=1:Q
  freqBig = expGeneHits(m)/expNum;
  for k=1:K
    freqSmall = clusterHits(m,k)/clusterSizes(k);
      clusterMatrixPVs(m,k)=BetterHygeccdf(clusterHits(m,k), ...
					   expNum, ...
					   expGeneHits(m), ...
					   clusterSizes(k));
  end
end

tfMatricesV = cell(size(tfMatrices));
for m=1:Q
  tfMatricesV{m}=sprintf('V$%s', tfMatrices{m});
end

[f1, allMatrixNames, f3, tfEntrezIDs, tfGeneNames]=textread('ThorssonTFsWithMatrices.tsv', ...
						  '%s %s %s %d %s', ...
						  'delimiter', '\t');

uniqueMatrixNames = unique(allMatrixNames);
M = length(uniqueMatrixNames);

% get a list of unique transcription factor genes from Vesteinns file
tfNames = unique(tfGeneNames);
T = length(tfNames);

tfInds = [];

% get the set of differentially expressed genes, that are
% associated with transcription factors for which we have TRANSFAC
% matrices
load 'TimeLaggedCorr2' 'tfInds' 'tfNames';
%load 'DiffExpGenesMore' 'psGeneSymbols';
%for t=1:T
%  geneSym = tfNames{t};
%  ind = find(1==strcmp(psGeneSymbols, geneSym));
%  if length(ind) > 0
%    tfInds = [tfInds, ind];
%  end
%end
%tfInds = unique(tfInds);
%tfNames = psGeneSymbols(tfInds);
%[tfNames, tfi]=sort(tfNames);
%tfInds = tfInds(tfi);
T = length(tfNames);

pcsPV = ones(T,K);
pctBind = zeros(T,K);
matNames = cell(T,K);
hits = zeros(T,K);

for m=1:M
  % get the matrix name (from the unique list of matrices in
  % Vesteinns file)
  matrixName = uniqueMatrixNames{m};

  %  Did we scan this matrix?
  % find this matrix in the list of matrices that I scanned, and
  % get the index if it appears in the list
  scanMatInd = find(1==strcmp(matrixName, tfMatricesV));
  if length(scanMatInd)==1

    % get the list of genes associated with this TF matrix, from
    % among the list of genes in Vesteinns file that are 
    % differentially expressed
    geneInd = find(1==strcmp(allMatrixNames, matrixName));
    G = length(geneInd);
    for g=1:G
      % get the gene symbol
      geneName = tfGeneNames{geneInd(g)};

      % get the index of this gene in the final list of genes
      tfi = find(1==strcmp(tfNames, geneName));
      if length(tfi)==1
	for k=1:K
	  if clusterMatrixPVs(scanMatInd, k) < pcsPV(tfi,k)
	    pcsPV(tfi,k)=clusterMatrixPVs(scanMatInd, k);
	    matNames{tfi,k} = matrixName;
	    pctBind(tfi,k) = clusterHits(scanMatInd,k)/ ...
		clusterSizes(k);
	    hits(tfi,k) = clusterHits(scanMatInd,k);
	  end
	end
      end
    end

    
  else
    matrixName
    warning 'unable to find matrix';
  end
end

save 'ClusterScans3' 'pcsPV' 'tfNames' 'tfInds' 'matNames' ...
    'pctBind' 'hits';

    

          
