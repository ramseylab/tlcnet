function BuildNetworkNew

load 'TimeLaggedCorr2';
load 'ClusterScans2';
load 'ClusterAnalysis';
load 'DiffExpGenes' 'psGeneSymbols';

pcsLPV = -log10(pcsPV);
tlcLPVC = -log10(tlcPVC);
toptLPVC = -log10(toptPVC);

chi2 = -2*log(pcsPV .* tlcPVC);

combLPVFinal = zeros(size(pcsLPV));
combLPVFinal(:) = -log10(gammainc(chi2/2, 2, 'upper')/gamma(2));

fid=fopen('NetworkTemp.tsv','w+');
fprintf(fid, '%%Cluster\tTF\tTF Clust\tLPV(scanning)\tMatrix\tPct Bind\tHits\tLPV(expression)\tTopt PV\tmean shift\tLPV(combined)\taverage corr\n');

M = size(tlcPVC,1);

outDeg = zeros(1,M);
avgCorr = zeros(M,K);
for k=1:K
  kidx = find(k==bestKidx);
  avgCorr(:,k)=mean(tlcOpt(:,kidx),2);
end

masterTFClust = GetClusterForGenes(tfNames);

ictr = 0;

networkMat = zeros(M,K);

% get a vector of the combined P values
combPVs = sort(10.^(-combLPVFinal(:)));

pvCutoff = 10^(-1.6);

meanTLCShifts = zeros(M,K);
for k=1:K
  kind = find(bestKidx == k);
  for m=1:M
    meanTLCShifts(m,k) = mean(tlcShifts(m,kind));
  end
end

for k=1:K
  ind = find(combLPVFinal(:,k) > -log10(pvCutoff));
  ind2 = find(pcsLPV(:,k) > 1.3);
  ind3 = find(meanTLCShifts(:,k) > 10);
  ind4 = find(abs(avgCorr(:,k)) > 0.0);
  ind = intersect(intersect(intersect(ind, ind2),ind3),ind4);

  kind = find(bestKidx == k);
  nk = length(kind);
  
  T = length(ind);
  for t=1:T
    ti = ind(t);
    
    tfName = tfNames{ti};

    tfind = find(1==strcmp(psGeneSymbols, tfName));
    tfclust = masterTFClust(ti);
    
    fprintf(fid, '%d\t%s\t%d\t%f\t%s\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n', ...
	    k, ...
	    tfNames{ti}, ...
	    tfclust, ...
	    pcsLPV(ti,k), ...
	    matNames{ti,k}, ...
	    pctBind(ti,k), ...
	    hits(ti,k), ...
	    tlcLPVC(ti,k), ...
	    -log10(toptPVC(ti,k)), ...
	    meanTLCShifts(ti,k), ...
	    combLPVFinal(ti,k), ...
	    avgCorr(ti,k));

    outDeg(ti) = outDeg(ti) + nk*pctBind(ti,k);

    networkMat(ti,k)=combLPVFinal(ti,k);
    ictr = ictr + 1;
  end
end
fclose(fid);

% fid=fopen('TFOutDeg2.tsv','w+');
% for m=1:M
%   fprintf(fid, '%s\t%d\n', tfNames{m}, outDeg(m));
% end
% fclose(fid);

'statistical power of time-lagged correlation (number of positives): '
length(find(tlcLPVC(:) > 3))

'statistical power of promoter scanning (number of positive): '
length(find(pcsLPV(:) > 3))

'number of interactions found: '
ictr

%save 'Network2' 'networkMat' 'avgCorr';
