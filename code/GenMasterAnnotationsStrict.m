function GenMasterAnnotationsStrict

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

strict = 1;
% Add Jackson Lab gene ontology annot
[psGOProcess, ...
 psGOComponent, ...
 psGOFunction]=ReadMGIAnnotGO('gene_association.mgi', ...
				    psGeneSymbols, ...
				    psGOProcess, ...
				    psGOComponent, ...
				    psGOFunction, ...
				    strict);

save 'MasterAnnotationsStrict';

