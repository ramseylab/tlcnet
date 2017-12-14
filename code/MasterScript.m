function MasterScript

% 
% system('cat Mouse430_2.na22.annot.txt | perl GetAnnotation.pl > MasterAnnotations.tsv');
% system('cat Mouse430_2.na22.annot.txt | perl GetAnnotationStrict.pl > MasterAnnotationsStrict.tsv');
%

% 
%---------------------------------------
% preliminary stuff:
%cd '../NormDat';
%system('NormalizeData.sh');
%---------------------------------------

cd '../Analysis';
% generate MAT files for expression data, GO terms, and annotations
RunSetup;

% generate a MAT file of all expressed genes 
RunGetAllGenes

% generate a TSV file of all expressed genes
RunGetExpGenes;
RunExpGenesFile;

% select differentially expressed genes, at 0.0001 P-value cutoff 
% (used for selecting target genes that are strongly differentially expressed)
RunGetDiffExpGenes;
RunDiffExpGenesFile;
GO=geneont('File', 'gene_ontology.obo');
GetListOfDiffExpTransRegs(GO, 1);
system('/usr/bin/perl ../Scripts/AddTFStatusToDiffExpGenes.pl > DiffExpGenesWithTFStatus.tsv');

% select differentially expressed genes, at 0.001 P-value cutoff
% (used for selecting transcription factors that are differentially expressed)
RunGetDiffExpTFGenes;
RunDiffExpTFGenesFile;
system('/usr/bin/perl ../Scripts/AddTRANSFACInfoToTFGenes.pl > DiffExpTFGenesWithTFInfo.tsv');

% generate a TSV file of the unique experiment names, and num
% replicates
GenUniqueExptNamesFile;
GenTableExperiments;

% cluster based on expression pattern and save results
RunClusterAnalysis;

% use simulated annealing to find the best order for the clusters
OrderClusters;

% perform GO enrichment analysis
RunGOAnalysis;
RunGOAnalysisFile;

% (1) Generate a file "TFs.txt" containing a list of all genes
% within DiffExpGenes.tsv that are transcription factors
% (2) Need to incorporate this information as a column in DiffExpGenes.tsv
% (3) Compute number of TF genes within each cluster

GenClusterTable;

GetListOfDiffExpTransRegs(GO, 0);
GenProbableTransRegs;
% generate time-lagged correlations
RunTimeLaggedCorr;
GenTimeCourseTable;

% ==============================================================
% Scan the promoters of all genes, and store the results
% in a directory "PromScan", in one GFF file for each cluster.
% ==============================================================
% perl GetBackgroundListOfGenes.pl
% perl GenClusterFastaFiles.pl
% <generate file TRANSFACMatricesDiffExp.txt>
% perl GetThreshScore.pl {0,1,2} 0.2
% 
% analyze the promoter scanning results
%RunAnalyzeClusterScans;

% build the combined network
%BuildNetwork;
