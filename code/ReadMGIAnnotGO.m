% Add the MGI (Jackson Lab) mouse genome GO annotations to the
% annotations obtained from the Mouse GeneChip

function [psgp, psgc, psgf]=ReadMGIAnnotGO(mgiAnnotFile, ...
					   psGeneSymbols, ...
					   psGOProcess, ...
					   psGOComponent, ...
					   psGOFunction, ...
					   strict)

annot = goannotread(mgiAnnotFile);
annot = struct2cell(annot);

N = length(annot);

psgp = psGOProcess;
psgc = psGOComponent;
psgf = psGOFunction;

for i=1:N
  geneName = annot{3,i};
  type = annot{9,i};
  goid = annot{5,i};
  evidence = annot{7,i};
  
  if ~strict || (~strcmp(evidence, 'IEA') && ~strcmp(evidence, 'RCA'))
    geneid = find(strcmp(geneName, psGeneSymbols)==1);
    if length(geneid) > 0
      for j=[geneid]
	if strcmp(type, 'P') == 1
	  gp = psgp{j};
	  psgp{j} = [];
	  psgp{j} = [ unique(cat(1,gp,goid)) ];
	else
	  if strcmp(type, 'C') == 1
	    gc = psgc{j};
	    psgc{j} = [];
	    psgc{j} = [ unique(cat(1,gc,goid)) ];
	  else
	    if strcmp(type, 'F')== 1
	      gf = psgf{j};
	      psgf{j} = [];
	      psgf{j} = [ unique(cat(1,gf,goid)) ];
	    end
	  end
	end
      end
      sprintf('adding GO id %d for gene %s\n', goid, geneName);
    else
      %%%%%    sprintf('could not find gene name: %s\n', geneName)
    end
  end
end

