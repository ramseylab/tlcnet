function RunSetup(arrayDataOnly)

psNames=textread('ProbeSetNames.txt','%s');
normDat=load('NormalizedReplicates.tsv');
[arrayNames,exptNames]=textread('ExperimentNames.tsv','%s %s');
uniqueExptNames = unique(exptNames);
U = length(uniqueExptNames);
numReplicates = zeros(1,U);
for i=1:U
  numReplicates(i) = length(find(strcmp(exptNames, uniqueExptNames{i})));
end
refCol = find(strcmp(uniqueExptNames, 'WT_UNSTIM_0000'));
save 'NormalizedReplicates' 'psNames' 'normDat' 'arrayNames' ...
    'exptNames' 'uniqueExptNames' 'numReplicates' 'refCol';
clear 'normDat' 'psNames' 'exptNames' 'arrayNames' 'araryDataOnly' ...
    'U' 'numReplicates' 'uniqueExptNames' 'i' 'refCol';
  
if nargin > 0 && arrayDataOnly
  return;
end

% get list of GO terms
GO=geneont('File', 'gene_ontology.obo');
[goTerms, ontCode, ontLevel, ontNames]=ReadGeneOntology(GO);
clear 'GO';
save 'GeneOntology';
clear 'goTerms' 'ontCode' 'ontLevel' 'ontNames';


% read Affymetrix annotations
annotData=ReadAffyAnnot('MasterAnnotationsStrict.tsv');
psNames=annotData{1};
psGeneSymbols=annotData{2};
psEntrezIDs=annotData{3};
psGOProcess=annotData{4};
psGOComponent=annotData{5};
psGOFunction=annotData{6};
psChromosome=annotData{7};
psStartCoord=annotData{8};
psEndCoord=annotData{9};
psTitles=annotData{10};
clear 'annotData';

% Add Jackson Lab gene ontology annot
strict = 1;
[psGOProcess, ...
 psGOComponent, ...
 psGOFunction]=ReadMGIAnnotGOStrict('gene_association.mgi', ...
				    psGeneSymbols, ...
				    psGOProcess, ...
				    psGOComponent, ...
				    psGOFunction, ...
				    strict);

save 'MasterAnnotationsStrict';



% read Affymetrix annotations
annotData=ReadAffyAnnot('MasterAnnotations.tsv');
psNames=annotData{1};
psGeneSymbols=annotData{2};
psEntrezIDs=annotData{3};
psGOProcess=annotData{4};
psGOComponent=annotData{5};
psGOFunction=annotData{6};
psChromosome=annotData{7};
psStartCoord=annotData{8};
psEndCoord=annotData{9};
psTitles=annotData{10};
clear 'annotData';

% Add Jackson Lab gene ontology annot
strict = 0;
[psGOProcess, ...
 psGOComponent, ...
 psGOFunction]=ReadMGIAnnotGO('gene_association.mgi', ...
			      psGeneSymbols, ...
			      psGOProcess, ...
			      psGOComponent, ...
			      psGOFunction, ...
			      strict);

save 'MasterAnnotations';