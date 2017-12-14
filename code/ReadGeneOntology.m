function [goTerms, ontCode, ontLevel, ontNames] = ReadGeneOntology(GO)

gt = GO.Terms;
M = gt(end).id;
goTerms = cell(1,M);
ontCode = zeros(1,M);
ontLevel = zeros(1,M);
Q = length(gt);
maxLevel = 15;
ontNames = {'biological process';
            'cellular component';
            'molecular function'};
for i=1:Q
  if 100*ceil(i/100)==i
    i
  end
  gti = gt(i);
  id = gti.id;
  goTerms{id}=gti.name;
  ont = gti.ontology;
  lastnum = 0;
  for i=0:maxLevel
    nextnum = length(getrelatives(GO, id, 'Height', i));
    if nextnum == lastnum
      break;
    end
    lastnum = nextnum;
  end
  ontLevel(id)=i;
  ontId = find(strcmp(ont, ontNames));
  if length(ontId) == 0
    ontId = 0;
  end
  ontCode(id) = ontId;
end

