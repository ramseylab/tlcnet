function ComputeKnockoutTFEffect(tfName)

load 'DiffExpGenes';
load 'ClusterAnalysis';
load 'TimeLaggedCorr';


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

tind = find(strcmp(tfNames, tfName));
if ~tind
  error 'could not find tf';
end

myd88Relexp = tfDat(tind,myd88Cols) - tfDat(tind,wtMCols)

trifRelExp = tfDat(tind,trifCols) - tfDat(tind,wtTCols)



