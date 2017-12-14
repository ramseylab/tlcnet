function AnalyzeTFMotifMap

[mats, tfgs] = textread('TFMotifMap','%s %s', 'delimiter', '\t');
umats = unique(mats);
ngs = [];
T = length(umats);
for i=1:T
  umat = umats{i};
  inds = find(strcmp(mats, umat));
  ngs = [ngs length(unique(tfgs(inds)))];
end
mean(ngs)

