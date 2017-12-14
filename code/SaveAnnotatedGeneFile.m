% SaveAnnotatedGeneFile
%
% From a list of probesets, print out a spreadsheet of a sorted,
% annotated list of the genes that corresponds to this probeset.
%
% "psInd" is the list of probeset indices that identifies which
% probesets are to be included in the output.
%
% "means" is an optional matrix of all microarray data, that
% can be included in the output 
%
% "doCluster" is optional, allowing you to additionally specify
% that the probesets should be clustered based on fold-change
% relative to wild-type
function SaveAnnotatedGeneFile(psInd, ...
			   geneSymbols, ...
			   geneTitles, ...
			   entrezIDs, ...
			   probeSets, ...
			   outFile, ...
			   goProcess, ...
			   goComponent, ...
			   goFunction, ...
			   goTerms, ...
			   means, ...
			   stdDevs, ...
			   uniqueExptNames, ....
			   numReplicates, ...
			   colSets, ...
			   colSetNames, ...			   
			   psClust)


[geneNames, geneInd] = sort(geneSymbols(psInd));

fid=fopen(outFile, 'w+');
fprintf(fid, 'Gene symbol\tEntrez ID\tProbeset(s)');
fprintf(fid, '\tRepresentative probeset\tGene title');
if nargin > 16
  fprintf(fid, '\tCluster');
end

fprintf(fid, '\tGO Process\tGO Component\tGO Function');

selmeans = means(psInd,:);
lmeans = selmeans;
N=size(selmeans,1);
M=length(uniqueExptNames);
for j=1:M
  fprintf(fid, '\t%s', uniqueExptNames{j});
end

fprintf(fid, '\tMax Exp');

if nargin > 14
  % user wants colSet summaries
  S=length(colSets);
  colSetVals = zeros(N,S);
  for j=1:S
    cj=colSets{j};
    colSetVals(:,j)=max(lmeans(:,cj),[],2);
    fprintf(fid, '\t%s', colSetNames{j});
  end
end

fprintf(fid, '\n');
N=length(geneNames);

% go through the list of unique gene names
for i=1:N
  geneName = geneNames{i};

  % the list of indices (indices within psInd) of probesets
  % corresponding to this gene
  psForGeneRel = geneInd(i);
  
  % the corresponding indices within "psNames", of these probesets
  probeSetsForGene = psInd(psForGeneRel);
  
  psNamesStr = '';
  numPS = length(probeSetsForGene);
  entrezIDList = [];

  gopAr = [];
  gocAr = [];
  gofAr = [];

  % go through the list of probesets corresponding to this gene
  for j=1:numPS
    psi = probeSetsForGene(j);
    
    entrezIDList = cat(2, entrezIDList, entrezIDs{psi});
    
    % append this probeset to the probeset names string
    psNamesStr = strcat(psNamesStr, ...
			probeSets{psi});
    if j < numPS
      psNamesStr = strcat(psNamesStr, ', ');
    end

    gopAr = [gopAr, goProcess{psi}];
    gocAr = [gocAr, goComponent{psi}];
    gofAr = [gofAr, goFunction{psi}];
  end

  gopAr = unique(gopAr);
  gocAr = unique(gocAr);
  gofAr = unique(gofAr);

  gopstr = '';
  for k=1:length(gopAr)
    goid = gopAr(k);
    gopstr = strcat(gopstr, sprintf('GO:%07d (%s)', ...
				    goid, ...
				    goTerms{goid}));
    if k < length(gopAr)
      gopstr = strcat(gopstr, ', ');
    end
  end

  gocstr = '';
  for k=1:length(gocAr)
    goid = gocAr(k);
    gocstr = strcat(gocstr, sprintf('GO:%07d (%s)', ...
				    goid, ...
				    goTerms{goid}));
    if k < length(gocAr)
      gocstr = strcat(gocstr, ', ');
    end
  end

  gofstr = '';
  for k=1:length(gofAr)
    goid = gofAr(k);
    gofstr = strcat(gofstr, sprintf('GO:%07d (%s)', ...
				    goid, ...
				    goTerms{goid}));
    if k < length(gofAr)
      gofstr = strcat(gofstr, ', ');
    end
  end
  
  entrezIDStr = '';
  entrezIDList = unique(entrezIDList);
  Q = length(entrezIDList);
  for k=1:Q
    entrezIDStr = strcat(entrezIDStr, sprintf('%d', ...
					      entrezIDList(k)));
    if k < Q
      entrezIDStr = strcat(entrezIDStr, ', ');
    end
  end
  
  geneTitle = '';
  
  % we are printing expression data into this spreadsheet
  genePSExp = selmeans(psForGeneRel,:);
  chosenPSI = psForGeneRel;
  
  fprintf(fid, ...
	  '%s\t%s\t%s', geneName, entrezIDStr, psNamesStr);
  
  % print out the representative probeset
  fprintf(fid, ...
	  '\t%s', probeSets{psInd(chosenPSI)});

  geneTitle = geneTitles{psInd(chosenPSI)};

  % print out the gene title
  fprintf(fid, '\t%s', geneTitle);
  
  if nargin > 16
    % print out the cluster of this probeset
    fprintf(fid, '\t%d', psClust(psInd(chosenPSI)));
  end

  % print GO function information
  fprintf(fid, ...
	  '\t%s\t%s\t%s', ...
	  gopstr, ...
	  gocstr, ...
	  gofstr);

  % print out the gene expression values
  geneExp = lmeans(chosenPSI, :);
  for j=1:M
    fprintf(fid, '\t%f', geneExp(j));
  end

  fprintf(fid, '\t%f', max(geneExp));
  
  % print out the column-set summaries
  if nargin > 14
    for k=1:S
      fprintf(fid, '\t%f', colSetVals(chosenPSI, k));
    end
  end
  
  fprintf(fid, '\n');


end

fclose(fid);


