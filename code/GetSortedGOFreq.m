% Take the p-values for gene ontology frequencies computed by
% "AnalyzeGeneOntFreq.m", and print out a sorted list of
% high-significance gene ontology terms.
%
% Stephen Ramsey
% Institute for Systems Biology
% June 2006
%
% Let M be the max GOID number.  Let N be the number of GO terms.
%
% goCountAll:  length M vector of counts of GO terms, 
%              in the larger set of genes
% goCountSub:  length M vector of counts of GO terms,
%              in the smaller set of genes
% goFreqAll:   length M vector of frequencies of GO terms,
%              in the larger set of genes
% goFreqSub:   length M vector of frequencies of GO terms,
%              in the smaller set of genes
% ontCode:     length M vector of numbers indicating the
%              ontology of each GO term (1 = process, 
%              2 = component, and 3 = function)
% goCountLPV:  length N vector of the log10 p-values of
%              the change in GO count from the larger set
%              to the smaller set of genes
% goIDs:       length N vector of GO ID numbers for each
%              GO term
% goTerms:     length N cell arary of GO terms
% goFreqFile:  the name of the TSV file in which the output
%              will be saved

function GetSortedGOFreq(goCountAll, ...
			      goCountSub, ...
			      goFreqAll, ...
			      goFreqSub, ...
			      ontCode, ...
			      goCountLPV, ...
			      goIDs, ...
			      goTerms, ...
			      goFreqFile)

N=length(goCountLPV);

fid=fopen(goFreqFile, 'w+');


[featSort, featSortInd]=sort(goCountLPV);

featSort = featSort([N:-1:1]);
featSortInd = featSortInd([N:-1:1]);

ont = {'process', 'component', 'function'};

for i=1:100
  goid = goIDs(featSortInd(i));
  fprintf(fid, '%f\tGO:%07d\t%f\t%f\t%d\t%d\t%s\t%s\n', featSort(i), ...
	  goid, ...
	  goCountAll(goid), ...
	  goCountSub(goid), ...
	  goFreqAll(goid), ...
	  goFreqSub(goid), ...
	  ont{ontCode(goid)}, ...
	  goTerms{featSortInd(i)});
end

fclose(fid);
