% Computes the rate at which a given set of GO terms appears in a
% list of genes.  
function [stats, goCount, goFreq]=ComputeGeneOntFreq(geneList, ...
						     psGOProcess, ...
						     psGOComponent, ...
						     psGOFunction, ...
						     ontCode, ...
						     GO)

% get the connection to the Gene Ontology database

N=length(geneList);
M = length(ontCode);

goCount = zeros(M,1);

haveProcess = 0;
haveFunction = 0;
haveComponent = 0;

haveAny = 0;
for i=1:N
  geneInd = geneList(i);

  anyFlag = 0;
  
  gp = psGOProcess{geneInd};
  if length(gp) > 0
    haveProcess = haveProcess + 1;
    anyFlag = 1;
  end
  
  gc = psGOComponent{geneInd};
  if length(gc) > 0
    haveComponent = haveComponent + 1;
    anyFlag = 1;
  end
  
  gf = psGOFunction{geneInd};
  if length(gf) > 0
    haveFunction = haveFunction + 1;
    anyFlag = 1;
  end

  if anyFlag
    haveAny = haveAny + 1
  end
end

haveProcess
haveComponent
haveFunction
haveAny

ancestCache = cell(M,1);

%tic
for i=1:N
%  if 1000*ceil(i/1000) == i
%    toc
%    i
%    tic
%  end
  
  geneInd = geneList(i);
  
  gp = psGOProcess{geneInd};
  gc = psGOComponent{geneInd};
  gf = psGOFunction{geneInd};  

  ga = cat(1, gp(:), cat(1, gf(:), gc(:)));
  
  gta = [];
  
  if size(ga, 2) > 1
    error 'invalid size for ga';
  end
  
  for j=ga'
    if length(ancestCache{j})==0
      try
	ancest = getancestors(GO, j);
      catch
	warning 'illegal GO ID: '
	j
	j=0;
	ancest = [];
      end
      if j ~= 0
	ancestCache{j}=ancest;
      end
    else
      ancest = ancestCache{j};
    end
   
    gta = cat(1, gta, ancest);
    
  end
  gta = unique(gta);
  goCount(gta)=goCount(gta)+1;
end

stats = [haveProcess, haveComponent, haveFunction];

goFreq = zeros(M,1);

% If there are any GO terms for which the count is greater
% than the number of genes "N" that have *any* annotation for
% that particular ontology, set the count to N.  Ideally this
% should not ever occur, but there must be some GO terms that
% are members of more than one ontology, or something, because
% some of the GO terms, like "biological process", end up having
% a count that is a few greater than N.
ind1 = find(ontCode > 0);
ind = find(goCount(ind1) > stats(ontCode(ind1))');
goCount(ind1(ind)) = stats(ontCode(ind1(ind)));

% compute GO term occurrence frequencies
ind = find(ontCode == 1);
goFreq(ind) = goCount(ind) / N;

ind = find(ontCode == 2);
goFreq(ind) = goCount(ind) / N;

ind = find(ontCode == 3);
goFreq(ind) = goCount(ind) / N;

% haveProcess
% haveFunction
% haveComponent

