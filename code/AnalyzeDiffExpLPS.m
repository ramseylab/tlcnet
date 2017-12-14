function AnalyzeDiffExpLPS

load 'DiffExpGenes';

[lpsCols] = GetColsets('WT','LPS',0);

M = length(lpsCols);
refExp = means(:,lpsCols)-repmat(means(:,refCol),1,M);
length(find(sign(Extreme(refExp'))>0))
