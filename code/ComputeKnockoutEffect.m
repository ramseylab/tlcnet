function ComputeKnockoutEffect(clust)

load 'DiffExpGenes';
load 'ClusterAnalysis';

[myd88Cols,myd88Times] = GetColsets('MYD88','LPS',0,120)
[trifCols,trifTimes] = GetColsets('TRIF','LPS',0,120)
[wtCols,wtTimes] = GetColsets('WT','LPS',0,120);


[myd88Times,iwt,im] = intersect(wtTimes, myd88Times);

wtMCols = wtCols(iwt)
myd88Cols = myd88Cols(im)
myd88Times

[trifTimes,iwt,it] = intersect(wtTimes, trifTimes);
wtTCols = wtCols(iwt)
trifCols =  trifCols(it)
trifTimes

kind = find(bestKidx == clust);

myd88Relexp = median(means(kind,myd88Cols) - means(kind,wtMCols),1)

trifRelExp = median(means(kind,trifCols) - means(kind,wtTCols),1)



