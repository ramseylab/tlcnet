function outData=AnalyzeGO(geneListSuper, ...
			   geneListSub, ...
			   threshLPV, ...
			   GO)

outData = cell(1,1);

load 'MasterAnnotations';
load 'GeneOntology';

N = length(psGeneSymbols);

if nargin < 3
  GO = geneont('File', '../Informatics/gene_ontology.obo');
end

bigN = length(geneListSuper);
smallN = length(geneListSub);

[allStats, allCount, allFreq]=ComputeGeneOntFreq(geneListSuper, ...
						 psGOProcess, ...
						 psGOComponent, ...
						 psGOFunction, ...
						 ontCode, ...
						 GO);
[smallStats, smallCount, smallFreq]=ComputeGeneOntFreq(geneListSub, ...
						  psGOProcess, ...
						  psGOComponent, ...
						  psGOFunction, ...
						  ontCode, ...
						  GO);

goids = find(ontCode);

lpvs = GeneOntLPV(smallCount(goids), ...
		  allStats(ontCode(goids)), ...
		  allCount(goids), ...
		  smallStats(ontCode(goids)));

smallCount(3677)
allCount(3677)
allStats(ontCode(3677))
smallStats(ontCode(3677))

lpvs(6355)

[slpvs, slpvi]=sort(lpvs);
numT = length(find(lpvs > threshLPV));
outStr = '';
M = length(goids);
ctr = 1;
for j=M:-1:(M-numT)
  goid = goids(slpvi(j));
  ontl = ontLevel(goid);
  smallCt = smallCount(goid);
  superCt = allCount(goid);
  ont = ontCode(goid);
  if ontl > 2
    superFreq = superCt/bigN;
    smallFreq = smallCt/smallN;
    if smallFreq > superFreq
      if smallCount(goid) > 1
	if smallFreq > 0.05
	  lineStr = sprintf('%s\t%f\t%d\tGO:%d\t%d', ...
			    goTerms{goid}, ...
			    slpvs(j), ...
			    ontl, ...
			    goid, ...
			    smallCt);
	  outStr = strvcat(outStr, lineStr);
	  outStr
	  outData{ctr} = {goid;
			 goTerms{goid}; 
			 slpvs(j);
			 ontl;
			 smallCt;
			 superFreq;
			 smallFreq};
	  ctr = ctr + 1;
	end
      end
    end
  end
end

