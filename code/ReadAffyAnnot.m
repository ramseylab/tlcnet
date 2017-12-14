% Reads the annotations data file that is created by the script 
% "GetAnnotations.pl", from the Affymetrix Mouse GeneChip
% master annotation file.
% 
% Stephen Ramsey
% Institute for Systems Biology
% Aug 2006
function retData=ReadAffyAnnot(fileName)

[probeSets, geneSymbols, geneIDs, ...
 goProcess, goComponent, goFunction, ...
 chromosomeS, startCoordS, endCoordS titles]=textread(fileName, ...
					      '%s %s %s %s %s %s %s %s %s %s', ...
					      'delimiter','\t',...
					      'whitespace','',...
					      'commentstyle','matlab',...
					      'bufsize',65535);

retData = {};
retData{1}=probeSets';
retData{2}=geneSymbols';
N=length(probeSets);
geneIDsCell = {};
goProcessCell = {};
goComponentCell = {};
goFunctionCell = {};
chromosome = {};
startCoord = {};
endCoord = {};
for i=1:N
geneIDsCell{i}=(eval(sprintf('[ %s ]',geneIDs{i})));;
goProcessCell{i}=(eval(sprintf('[ %s ]',goProcess{i})))';
goComponentCell{i}=(eval(sprintf('[ %s ]',goComponent{i})))';
goFunctionCell{i}=(eval(sprintf('[ %s ]',goFunction{i})))';
chromosome{i} = (eval(sprintf('[ %s ]',chromosomeS{i})))';
startCoord{i} = (eval(sprintf('[ %s ]',startCoordS{i})))';
endCoord{i} = (eval(sprintf('[ %s ]',endCoordS{i})))';
end
retData{3}=geneIDsCell;
retData{4}=goProcessCell;
retData{5}=goComponentCell;
retData{6}=goFunctionCell;
retData{7}=chromosome;
retData{8}=startCoord;
retData{9}=endCoord;
retData{10}=titles';

