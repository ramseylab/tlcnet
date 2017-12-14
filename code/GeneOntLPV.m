function lpv=GeneOntLPV(countSmallSet, ...
			allBigSet, ...
			countBigSet, ...
			allSmallSet)

pv = BetterHygeccdf(countSmallSet, ...
		    allBigSet, ...
		    countBigSet, ...
		    allSmallSet);
 
ind = find(isnan(pv));
if length(ind) > 0
  error 'invalid result from hygepdf';
end

ind = find(countBigSet == 0);
pv(ind) = 1;

ind = find(pv == 0);
pv(ind)=min(pv(find(pv)));

lpv = -log10(pv);

