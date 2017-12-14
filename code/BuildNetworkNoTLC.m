function BuildNetworkNoTLC

load 'TimeLaggedCorr2';
load 'ClusterScans2';
load 'ClusterAnalysis';
load 'DiffExpGenes' 'psGeneSymbols';

tlcCombPVC = (ones(size(tlcPVC)) .* toptPVC);

pcsLPV = -log10(pcsPV);
tlcCombLPVC = -log10(tlcCombPVC);
toptLPVC = -log10(toptPVC);

p1 = pcsPV;
p2 = ones(size(tlcPVC));
w1 = 1;
w2 = 1;
pprod = p1.*p2;
chi2 = -2*log(pprod);

combLPVFinal = -log10(pcsPV);

%fid=fopen('Network2.tsv','w+');
%fidcsv=fopen('NetworkBiotap.csv', 'w+');
%fprintf(fid, '%%Cluster\tTF\tTF Clust\tLPV(scanning)\tMatrix\tPct Bind\tHits\tLPV(expression)\tTopt PV\tmean shift\tLPV(combined)\taverage corr\n');

M = size(tlcPVC,1);

outDeg = zeros(1,M);
avgCorr = zeros(M,K);
for k=1:K
  kidx = find(k==bestKidx);
  avgCorr(:,k)=mean(tlcOpt(:,kidx),2);
end

masterTFClust = GetClusterForGenes(tfNames);

%fprintf(fidcsv, '\"# Model Commands\",,,,,,,,,,\n');
%fprintf(fidcsv, '\"# Command Type\",\"Model Name\",\"Parent Model\",,,,,,,,\n');
%fprintf(fidcsv,'\"model\",\"root\",,,,,,,,,\n');
%fprintf(fidcsv,'\"model\",\"Clusters\",\"root\",,,,,,,,\n');
%fprintf(fidcsv, ',,,,,,,,,,\n');
%for k=1:K
%  fprintf(fidcsv, '\"region\",\"Clusters\",\"C%d\",\"C%0.2d\",,,,,,,\n', ...
%	  k, k);
%end

ictr = 0;

networkMat = zeros(M,K);

% get a vector of the combined P values
combPVs = sort(10.^(-combLPVFinal(:)));
pvCutoff = EstimateFDR(combPVs', 0.025)
%pvCutoff = 10^(-1.6);

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
  ind = intersect(intersect(ind, ind2),ind3);

  kind = find(bestKidx == k);
  nk = length(kind);
  
  T = length(ind);
  for t=1:T
    ti = ind(t);
    
    tfName = tfNames{ti};

    tfind = find(1==strcmp(psGeneSymbols, tfName));
    tfclust = masterTFClust(ti);
  
    
%    fprintf(fid, '%d\t%s\t%d\t%f\t%s\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n', ...
% 	    k, ...
% 	    tfNames{ti}, ...
% 	    tfclust, ...
% 	    w1*pcsLPV(ti,k), ...
% 	    matNames{ti,k}, ...
% 	    pctBind(ti,k), ...
% 	    hits(ti,k), ...
% 	    w2*tlcCombLPVC(ti,k), ...
% 	    -log10(toptPVC(ti,k)), ...
% 	    meanTLCShifts(ti,k), ...
% 	    combLPVFinal(ti,k), ...
% 	    avgCorr(ti,k));

    outDeg(ti) = outDeg(ti) + nk*pctBind(ti,k);

%    fprintf(fidcsv, '\"general\",\"Clusters\",\"gene\",\"%s\",\"gene\",\"C%d\",\"positive\",\"C%0.2d\",\"C%0.2d\",,\n', ...
%	    tfNames{ti}, ...
%	    k, ...
%	    tfclust, ...
%	    k);
    
    networkMat(ti,k)=combLPVFinal(ti,k);
    ictr = ictr + 1;
  end
end
%fclose(fid);
%fclose(fidcsv);

%fid=fopen('TFOutDeg2.tsv','w+');
%for m=1:M
%  fprintf(fid, '%s\t%d\n', tfNames{m}, outDeg(m));
%end
%fclose(fid);

'statistical power of time-lagged correlation (number of positives): '
length(find(w2*tlcCombLPVC(:) > 3))

'statistical power of promoter scanning (number of positive): '
length(find(w1*pcsLPV(:) > 3))

'number of interactions found: '
ictr

%save 'Network2No' 'networkMat' 'avgCorr';
