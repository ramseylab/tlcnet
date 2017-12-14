function EstimateTLCUncertainty(gene1, gene2)

load 'DiffExpGenes';
load 'ClusterAnalysis' 'normRelMeans';

shifts = [0:10:90];

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

g1ind = find(strcmp(gene1, psGeneSymbols));
g2ind = find(strcmp(gene2, psGeneSymbols));
g1dat = normRelMeans(g1ind,:);
g2dat = normRelMeans(g2ind,:);

g1noise = 2*stdDevs(g1ind,:)./means(g1ind,:);
g2noise = 2*stdDevs(g2ind,:)./means(g2ind,:);
for i=1:K
  ci = tlcCols{i}(1);
  g1noise(ci) = 0;
  g2noise(ci) = 0;
end

iter = 100;
tlcopts = zeros(1,iter);
for i = 1:iter
  g1datnew = g1dat + g1noise.*randn(size(g1dat));
  g2datnew = g2dat + g2noise.*randn(size(g2dat));
  tlc = TimeLaggedCorr(g1datnew, g2datnew, shifts, tlcCols, tlcTimes);
  [maxcc,maxi] = max(abs(tlc));
  tlcopt = shifts(maxi);
  tlcopts(i) = tlcopt;
end

std(tlcopts)
mean(tlcopts)
