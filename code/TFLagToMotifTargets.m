function TFLagToMotifTargets

% get the number of clusters
load 'ClusterAnalysis' 'K';
load 'TimeLaggedCorr' 'tfNames' 'tfProbesets';
load 'DiffExpGenes' 'psNames' 'psGeneSymbols';

% load the list of scanned matrices
matrixNames = textread('ScannedMatrices.tsv','%s');

M = length(matrixNames);

T = length(tfNames);

tfTargets = cell(1,T);
tfTargetScores = cell(1,T);
for i=1:T
  tfTargets{i} = [];
  tfTargetScores{i} = [];
end

% load thorsson map
[foo1, vmatrixNames, foo2, foo3, tfGenes] = textread('ThorssonTFsWithMatrices.tsv',...
						  '%s %s %s %s %s');

% build a mapping between RefSeq IDs and Gene indices

% go through each matrix, one at a time
for m = 1:M
  matrixName = matrixNames{m};
  disp(sprintf('Analyzing matrix: %s\n', matrixName))
  vMatrixName = ['V$' matrixName];

  mind = find(strcmp(vmatrixNames, vMatrixName));
  if length(mind)==0
    matrixName
    error 'unable to find indices for matrix name';
  end
  tfInds = [];
  % go through each gene that is mapped to this matrix
  for i=1:length(mind)
    tfGene = tfGenes{mind(i)};
    tfGeneInd = find(strcmp(tfNames, tfGene));
    if length(tfGeneInd) == 1
      tfInds = [tfInds tfGeneInd];
    else
      if length(tfGeneInd) > 1
	error 'more than one index found';
      end
    end
  end
  
  if length(tfInds > 0)
    % we found at least one TF for this index, so it is worthwhile
    % to check the scanning results
    % go through each cluster, one at a time
    for k=1:K
      fileName = sprintf('../PromScan/GFF/%s_C%d.fa.gff.processed.tsv', matrixName, k);
      [gffProbesets, startCoords, endCoords, strands, scores, refSeqIDs] ...
	  = textread(fileName, '%s %d %d %s %f %s');
      
      F = length(gffProbesets);
      % go through each feature in the file
      for f=1:F
	% get the probeset
	probeset = gffProbesets{f};
	
	score = scores(f);

	% find the index of this probeset in the 
	psInd = find(strcmp(psNames, probeset));
	if length(psInd) == 0
	  probeset
	  error 'unable to find probeset index';
	end

	for i = tfInds
	  if length(intersect(tfTargets{i}, psInd))==0
	    tfTargets{i} = [tfTargets{i} psInd];
	    tfTargetScores{i} = [tfTargetScores{i} score];
	  end
	end
      end
    end
  end
end
 
load 'TimeLaggedCorr' 'tlcOpt' 'tlcShifts' 'tlcPV';
load 'ClusterAnalysis' 'bestKidx';

fid=fopen('TFLagToMotifTargets2.tsv','w+');
for i=1:T
  targets = tfTargets{i};
  targetScores = tfTargetScores{i};
  H = length(targets);
  for j=1:H
    tgi = targets(j);
    fprintf(fid, '%s\t%s\t%s\t%s\t%d\t%f\t%d\t%f\t%f\n', ...
	    tfNames{i}, ...
	    tfProbesets{i}, ...
	    psGeneSymbols{tgi}, ...
	    psNames{tgi}, ...
	    bestKidx(tgi), ...
	    tlcOpt(i,tgi), ...
	    tlcShifts(i,tgi), ...
	    targetScores(j), ...
	    tlcPV(i,tgi));
  end
end
fclose(fid);
