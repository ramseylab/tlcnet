function AddOntToGO

load GeneOntology;

fid=fopen('GOAnalysis2.tsv','w+');

[clust, goids, goterm, lpv, level, smallcount, superfreq, smallfreq]=...
    textread('GOAnalysis.tsv','%d GO:%d %s %f %d %d %f %f', ...
	     'commentstyle','matlab', ...
	     'delimiter', '\t');

L = length(clust);

for i=1:L
  goid = goids(i);
  
  ont = ontCode(goid);
  ontStr = '';
  if ont == 1
    ontStr = 'process';
  else 
    if ont == 2
      ontStr = 'component';
    else
      if ont == 3
	ontStr = 'function';
      end
    end
  end
  
  fprintf(fid, '%d\tGOID:%0.7d\t%s\t%s\t%f\t%d\t%d\t%f\t%f\n', ...
	  clust(i), ...
	  goid, ...
	  goterm{i}, ...
	  ontStr, ...
	  lpv(i), ...
	  level(i), ...
	  smallcount(i), ...
	  superfreq(i), ...
	  smallfreq(i));
end

fclose(fid);

