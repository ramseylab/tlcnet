function RunCombineReplicates
% Preliminaries:
% (1) Get the Affymetrix mouse GeneChip annotation file, in CSV format
% (2) Open the file in Excel, and save it as "text" (tab delimited)
%     - this produces a "TSV"file 
% (3) Run the TSV file through "GetAnnotations.pl" (in the Scripts
%     directory)
% (4) The "ReadAffyAnnot.m" MATLAB function is used to load in
%     the annotations (see below)

load 'MasterAnnotations.mat';

% read the collection of Affymetrix probeset names
indProbesToUse=find(strncmp(psNames,'AFFX-',5)==0);
sprintf('number of probes in data set: %d\n', length(indProbesToUse))

% read the microarray data.  This data comes from R, so it is a TSV
% file and not a MAT file
load 'NormalizedReplicates.mat';

psNames = psNames(indProbesToUse);

% Read the experiment names data file; the first column is the
% experiment name in SBEAMS; the second column is the simplified
% experiment name that I typed in (replicates are given the same name
% in the simplified column).  In the simplified experiment names, two
% CEL data files with the same strain/condition/time-point will have
% the same "simplified experiment name" (e.g., MYD88_LPS_060)
[dataFileNames, experimentNames] = textread('ExperimentNames.tsv', '%s %s');
uniqueExptNames = unique(experimentNames);
refCol = find(strcmp(uniqueExptNames, 'WT_UNSTIM_0000'))

retData=CombineAllReplicates(normDat(indProbesToUse,:),...
			     experimentNames);

numReplicates=retData{2};
means=retData{3};
stdDevs=retData{4};

clear('retData', 'normDat');

nu=length(uniqueExptNames);
fid=fopen('UniqueExptNames.txt', 'w+');
for i=1:nu
  fprintf(fid, '%s\n', uniqueExptNames{i});
end
fclose(fid);

save 'CombinedReplicates' ...
     'means' ...
     'numReplicates' ...
     'stdDevs' ...
     'uniqueExptNames' ...
     'psNames' ...
     'refCol' ...
     'indProbesToUse';
return;

