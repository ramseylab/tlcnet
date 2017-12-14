function [keptIndices, ...
	  rcMeans, ...
	  rcStdDevs, ...
	  pValues, ...
	  pValueCutoffs, ...
	  tcExptNames]=...
                    CombineDataAndSelectProbesets(arrayDataLogScale, ...
						  probeSetNames, ...
						  entrezGeneIDAnnotations, ...
						  minLogIntensity, ...
						  exptNames, ...
						  arrayNames, ...
						  stims, ...
						  strains, ...
						  numIterBootstrap, ...
						  smoothingLengthGOFScore, ...
						  numBinsGOFScore, ...
						  falseDiscoveryRate)

N = length(probeSetNames);

% encode the probeset type into an integer score
probeSetScores = zeros(1,N);
for i=1:N
  score = 4;
  psName = probeSetNames{i};
  if length(strfind(psName, '_x_'))>0
    score = 1;
  else
    if length(strfind(psName, '_s_'))>0
      score = 2;
    else
      if length(strfind(psName, '_a_'))>0
	score = 3;
      end
    end
  end
  probeSetScores(i)=score;
end
% only keep probesets with a score > 2
keepIndScore = find(probeSetScores > 2);

geneAnnotationsStr = {};
for k=1:N
  geneAnnotationsStr{k} = num2str(entrezGeneIDAnnotations{k});
end

% get the list of probesets that are annotated
keepIndAnnot = find(strcmp(geneAnnotationsStr, '')==0);

% replicate-combine the array data
retDat = CombineAllReplicates(arrayDataLogScale, ...
			      exptNames);
uniqueExptNames = retDat{1};
uniqueExptReplicates = retDat{2};
rcMeans = retDat{3};
rcStdDevs = retDat{4};

% get the probesets that have a log2 intensity greater than or
% equal to the minimum
keepIndLogInt = find(max(rcMeans,[],2)>= ...
		     minLogIntensity);


% get the list of probesets that passed the first set of filters
% (acceptable probeset type, minimum log intensity, and possessing
% an Entrez GeneID annotation)
keepInd = intersect(intersect(keepIndScore, keepIndAnnot), ...
		    keepIndLogInt);

sprintf('number of genes passing first filter: %d', length(keepInd))

if nargin > 5
  [diffExpGenes, pValuesGenes, pValueCutoffs, tcExptNames]= ...
      GetDiffExpGenes(arrayDataLogScale(keepInd, :), ...
		      arrayNames, ...
		      exptNames, ...
		      stims, ...
		      strains, ...
		      numIterBootstrap, ...
		      smoothingLengthGOFScore, ...
		      numBinsGOFScore, ...
		      falseDiscoveryRate);
  pValues = ones(N,length(tcExptNames));
  pValues(keepInd,:)=pValuesGenes;
  keepInd = keepInd(find(diffExpGenes));
end

keepIndNR = [];

% remove redundant probesets
keepGeneAnnot = geneAnnotationsStr(keepInd);
keepGeneAnnotUnique = unique(keepGeneAnnot);
Q = length(keepGeneAnnotUnique);
for q=1:Q
  uniqueGeneAnnot = keepGeneAnnotUnique{q};
  geneInd = find(strcmp(keepGeneAnnot, uniqueGeneAnnot)==1);
  origInd = keepInd(geneInd);
  if length(origInd > 1)
    if nargin > 5
      [mpv,bestInd] = min(min(pValues(origInd,:),[],2));
    else
      [mlg,bestInd] = max(max(rcMeans(origInd,:),[], ...
			      2));
    end
    keepIndNR = [keepIndNR origInd(bestInd)];
  else
    keepIndNR = [keepIndNR origInd];
  end
end

keptIndices = sort(keepIndNR);

