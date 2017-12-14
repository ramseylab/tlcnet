function GetListOfExpTransRegs(GO, strict)


load 'ExpGenes';
if nargin > 1 && strict
load 'MasterAnnotationsStrict' 'psGOFunction' 'psGOProcess';
else
load 'MasterAnnotations' 'psGOFunction' 'psGOProcess';
end

psGOFunction = psGOFunction(masterInds);
psGOProcess = psGOProcess(masterInds);

[mmat, vmat, tmat, entGeneID, geneSymbols]= ...
    textread('ThorssonTFsWithMatrices.tsv', ...
	     '%s %s %s %d %s');
thorssonUniqueGeneSymbols = unique(geneSymbols);
thorssonUniqueGeneIDs = unique(entGeneID);

[geneSymbols, geneDescs, geneAliases, junk, roachAnnot, junk2]= ...
    textread('RoachMasterList.tsv', ...
             '%s %s %s %s %s %[^\n]', ...
             'delimiter', '\t');

tfInd=union(union(find(strcmp(roachAnnot, 'DEFINITELY')), ...
		find(strcmp(roachAnnot, 'LIKELY'))), ...
		find(strcmp(roachAnnot, 'PROBABLY')));
roachGeneSymbols = unique(geneSymbols(tfInd));

tfGOIDs = [3700];
trGOIDs = [3712 9161 45449 6355];
load 'GeneOntology';

if nargin == 0
  GO=geneont('File','gene_ontology.obo');
end

N = length(masterInds);

if strict
  fid=fopen('ExpGenesTFStatusStrict.tsv','w+');
else
  fid=fopen('ExpGenesTFStatus.tsv','w+');
end  
tfInds = [];
psIsTR = zeros(1,N);
psIsTF = zeros(1,N);
for i=1:N
  isTF = 0;
  isTR = 0;
  geneName = psGeneSymbols{i};
  geneNameShort = geneName;
  k = strfind(geneName, ' /// ');
  if length(k) > 0
    geneNameShort = geneName(1:(k(1) - 1));
  end
  entrezGeneID = psEntrezIDs{i};
  if length(find(strcmp(geneNameShort, thorssonUniqueGeneSymbols))) > 0 ...
	|| length(intersect(entrezGeneID, thorssonUniqueGeneIDs)) > 0 ...
	|| length(find(strcmp(upper(geneNameShort), roachGeneSymbols))) ...
	> 0
    isTF = 1;
    isTR = 1;
  else
    goTerms = [psGOFunction{i}' psGOProcess{i}'];
    T = length(goTerms);
    for t=1:T
      goid = goTerms(t);
      ancestgoids = [];
      try
	ancestgoids = getancestors(GO, goid);
      catch
	warning 'illegal go id: '
	goid
	ancestgoids = [];
      end
      if length(intersect(ancestgoids, tfGOIDs)) > 0
	isTF = 1;
	isTR = 1;
      else
	if length(intersect(ancestgoids, trGOIDs)) > 0
	  isTR = 1;
	end
      end
    end
  end
  isTFStr = '';
  if isTF
    isTFStr = '1';
  end
  isTRStr = '';
  if isTR
    isTRStr = '1';
  end
  psIsTR(i) = isTR;
  psIsTF(i) = isTF;
  fprintf(fid, '%s\t%s\t%s\n', geneName, isTFStr, isTRStr);
end
fclose(fid);

if strict
  save 'ExpGenesTFStatusStrict' 'psIsTR' 'psIsTF' 'masterInds';
else
  save 'ExpGenesTFStatus' 'psIsTR' 'psIsTF' 'masterInds'; 
end 
