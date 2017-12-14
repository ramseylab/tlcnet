function GenTFComplexSIF

load 'TimeLaggedCorr2';
load 'Network2';

T = size(networkMat,1);
C = size(networkMat,2);

[p1s, int, p2s] = textread('Interactions2.sif','%s %s %s');
P = length(p1s);
pairsInt = zeros(T,T);
for i=1:P
  p1 = p1s{i};
  p2 = p2s{i};
  p1i = find(strcmp(tfNames, p1));
  if ~p1i
    p1
    error 'could not find p1';
  end
  p2i = find(strcmp(tfNames, p2));
  if ~p2i
    p2
    error 'could not find p2';
  end
  pairsInt(p1i,p2i)=1;
end


pairsGot = zeros(T,T);

fid=fopen('TFsCoregulatingClusters.sif','w+');
for i=1:C
  tfInClust = find(networkMat(:,i));
  N = length(tfInClust);
  for j=1:N
    for k=[ 1:(j-1), (j+1):N]
      ti1 = tfInClust(j);
      ti2 = tfInClust(k);
      if ~pairsGot(ti1,ti2) && ~pairsGot(ti2,ti1)
	t1 = tfNames{ti1};
	t2 = tfNames{ti2};
	pairsGot(ti1,ti2) = 1;
	pairsGot(ti2,ti1) = 1;
	if pairsInt(ti1,ti2)
	  fprintf(fid, '%s\tcc\t%s\n', t1, t2);
	end
      end
      
    end
  end
end
fclose(fid);
