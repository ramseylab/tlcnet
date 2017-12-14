function tlc=RunTimeLaggedCorr

shifts = [0:10:80];

[ WT_LPS_Cols, WT_LPS_Times] = GetColsets('WT', 'LPS', 0);
[ WT_PIC_Cols, WT_PIC_Times] = GetColsets('WT', 'PIC', 0);
[ WT_PAM3_Cols, WT_PAM3_Times] = GetColsets('WT', 'PAM3', 0);
[ WT_R848_Cols, WT_R848_Times] = GetColsets('WT', 'R848', 0);
[ ATF3_LPS_Cols, ATF3_LPS_Times] = GetColsets('ATF3', 'LPS', 0);
[ WT_PAM2_Cols, WT_PAM2_Times]=GetColsets('WT','PAM2',0);		
[ CREM_LPS_Cols, CREM_LPS_Times] = GetColsets('CREM', 'LPS', 0);
[ CREM_PIC_Cols, CREM_PIC_Times ] = GetColsets('CREM', 'PIC', 0);
[ ATF3_CPG_Cols, ATF3_CPG_Times ] = GetColsets('ATF3', 'CPG', 0);
[ ATF3_PIC_Cols, ATF3_PIC_Times ] = GetColsets('ATF3', 'PIC', 0);
[ ATF3_PAM2_Cols, ATF3_PAM2_Times ] = GetColsets('ATF3', 'PAM2', ...
						 0);

tlcCols = { WT_LPS_Cols; 
	    WT_PIC_Cols;
	    WT_PAM3_Cols;
	    WT_PAM2_Cols;
	    WT_R848_Cols;
	    ATF3_LPS_Cols;
	    CREM_LPS_Cols;
	    CREM_PIC_Cols;
	    ATF3_CPG_Cols;
	    ATF3_PIC_Cols; 
	    ATF3_PAM2_Cols};

tlcTimes = { WT_LPS_Times; 
	     WT_PIC_Times;
	     WT_PAM3_Times;
	     WT_PAM2_Times;
	     WT_R848_Times;
	     ATF3_LPS_Times;
	     CREM_LPS_Times; 
	     CREM_PIC_Times;
	     ATF3_CPG_Times;
	     ATF3_PIC_Times;
	     ATF3_PAM2_Times};

K = length(tlcCols);

for k=1:K
  [tvec, tind]=sort(tlcTimes{k});
  tlcTimes{k}=tvec;
  tlcCols{k}=tlcCols{k}(tind);
end

% load the set of differentially expressed genes, from which we are
% pulling our list of transcription factors
load 'DiffExpTFGenes' 'psGeneSymbols' 'means' 'psEntrezIDs' 'refCol';

% read transcription factors
[f1, f2, f3, tfEntrezIDs, tfGeneNames]=textread('ThorssonTFsWithMatrices.tsv', ...
					        '%s %s %s %d %s', ...
					        'delimiter', '\t');


% get unique list of TF gene names from Thorsson's spreadsheet
[uniqueTFs, utfInd] = unique(tfGeneNames);

numTFs = length(uniqueTFs)

tfInds = zeros(1,numTFs);

N = length(psEntrezIDs);

% get the index for each TF gene name
for t=1:numTFs
  t
  tfGeneName = uniqueTFs{t};
  tfEntrezID = tfEntrezIDs(utfInd(t));
  
  tfInd = find(strcmp(psGeneSymbols, tfGeneName)==1);
  if length(tfInd)==0
    for n=1:N
      if length(intersect(psEntrezIDs{n}, tfEntrezID))>0
	tfInd = n;
      end
    end
  end
  
  if length(tfInd) > 0
    tfInds(t) = tfInd;
  end
end

tfInds = tfInds(find(tfInds));
tfNames = psGeneSymbols(tfInds);

sprintf('number of transcription factors: %d', length(tfInds))

tfMeans = means(tfInds,:);
M=size(means,2);
tfMeansRef = repmat(tfMeans(:,refCol),1,M);
tfDat = (tfMeans-tfMeansRef)./repmat(max(abs(tfMeans-tfMeansRef),[],2),1,M);
clear 'means';

% load the set of differentially expressed genes, from which we are
% pulling our list of targets
load 'DiffExpGenes' 'psGeneSymbols' 'psGOProcess' 'psGOFunction';

% build a list of non-TF genes
[tfNamesFile] = textread('ProbableTransRegs.txt', '%s');
allTfInds = [];
ntf = length(tfNamesFile);
for i=1:ntf
  ind = find(1==strcmp(psGeneSymbols, tfNamesFile{i}));
  if length(ind) > 0
    allTfInds = [allTfInds, ind];
  end
end

N = length(psGeneSymbols);
nonTFInds = setxor([1:N], allTfInds);

nonTFIndsToUse = [];
numNon = length(nonTFInds);
for i=1:numNon
  ind = nonTFInds(i);
  if length(psGOProcess{ind}) > 2 && length(psGOFunction{ind}) > 2 ...
	&&  1 ~= strncmp(psGeneSymbols{ind}, 'Zfp', 3)
    nonTFIndsToUse = [nonTFIndsToUse, ind];
  end
end

nntf = length(nonTFIndsToUse);
sprintf('number of high-confidence non-TFs: %d', nntf)

% load expression clusters
load 'ClusterAnalysis' 'normRelMeans' 'bestKidx';

% get a list of non-transcription factor genes of which 
% half are upregulated, and half are downregulated
updown = sign(Extreme(normRelMeans(nonTFIndsToUse,:)'));
numUp = length(find(updown > 0));
numDown = length(find(updown < 0));
sprintf('number upregulated non-TFs:  %d  number downregulated: %d', ...
	numUp, numDown)
minUpDown = min(numUp, numDown);
balancedNonTFInds = [];
upctr = 0;
downctr = 0;
i = 0;
while i < 2*minUpDown
  ti = ceil(nntf*rand());
  if length(intersect(ti, balancedNonTFInds))==0
    if updown(ti) > 0
      if upctr < minUpDown
	balancedNonTFInds = [balancedNonTFInds, ti];
	upctr = upctr + 1;
      end
    else
      if downctr < minUpDown
	balancedNonTFInds = [balancedNonTFInds, ti];
	downctr = downctr + 1;
      end
    end
  end
  i = upctr + downctr;
end

balancedNonTFInds = unique(balancedNonTFInds);
sprintf('number of non-TF inds after balancing: %d', ...
	length(balancedNonTFInds))

nonTFIndsToUse = balancedNonTFInds;
nntf = length(balancedNonTFInds);

% compute the time-lagged correlation between TFs and all genes
tlc = TimeLaggedCorr(tfDat, ...
		     normRelMeans, ...
		     shifts, ...
		     tlcCols, ...
		     tlcTimes);

datNonTFs = normRelMeans(nonTFIndsToUse, :);

clear 'normRelMeans';

% compute the time-lagged correlation among the non-TFs
tlcNonTFs = TimeLaggedCorr(datNonTFs, ...
			   datNonTFs, ...
			   shifts, ...
			   tlcCols, ...
			   tlcTimes);

save 'TimeLaggedCorrBackgd' 'tlcNonTFs';

allCols = [];
K=length(tlcCols);
for k=1:K
  allCols = [allCols, tlcCols{k}'];
end
allCols = unique(allCols);
K=length(allCols);

% set the self-correlation to zero, for each non-TF gene
for i=1:nntf
  tlcNonTFs(:,i,i)=0;
end

B = 200;
sl = 0.005;

[tlcPV,tlcShifts,toptPV,tlcOpt] = TimeLaggedCorrPValues(tlc, ...
						  tlcNonTFs, ...
						  B, ...
						  shifts, ...
						  sl);

% the P-value for a TF-to-self association is always 1
M=size(tlc,2);
tfGeneInds = cell(1,M);
for m=1:M
  tfGeneName = tfNames{m};
  ind = find(1==strcmp(psGeneSymbols, tfGeneName));
  tfGeneInds{m}=ind;
  if length(ind)==1
    tlcPV(m,ind)=1;
  end
end

size(toptPV)

K = max(bestKidx);
tlcPVC = zeros(M, K);
tempPV = zeros(1,M);
toptPVC = zeros(M, K);

for k=1:K
  d = 0.18*K;
  ind = find(bestKidx == k);
  for m=1:M
    subt = length(intersect(ind, tfGeneInds{m}));

    chi2=-2*d*sum(log(tlcPV(m,ind)))/(length(ind)-subt);
    tlcPVC(m,k) = gammainc(chi2/2, d, 'upper');
    toptPVC(m,k) = 10.^(sum(log10(toptPV(m,ind)))/(length(ind)-subt));
  end
end

'statistical power (at 10^-3):'
length(find(tlcPVC(:).*toptPVC(:) < 0.001))

'statistical power (at 10^-2):'
length(find(tlcPVC(:).*toptPVC(:) < 0.01))

figure;
hist(-log10(tlcPVC(:)), 100);

save 'TimeLaggedCorr2' 'tlc' 'shifts' 'tlcCols' 'tlcTimes' 'tfInds' ...
     'tlcNonTFs' 'tlcPV' 'tlcShifts' 'tlcPVC' 'tfNames' 'tfDat' ...
    'toptPVC' 'tlcOpt';
