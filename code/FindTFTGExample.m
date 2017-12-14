function FindTFTGExample

% find a pair of genes (TF, TG) for which the standard  measure of
% TLC optimality gives a different time lag than the improved
% measure of optimality, but for which the rho-sq looks reasonable

load TimeLaggedCorr;
load DiffExpGenes;


T = length(tfNames);
N = length(psGeneSymbols);

maxRhoSqLags = zeros(T,N);

for t=1:T
  for i=1:N
    specTlc = squeeze(tlc(:,t,i).^2);
    [mactlc, maxtlci] = max(specTlc);
    maxRhoSqLags(t,i) = shifts(maxtlci);
  end
end

ind = find( (tlcOpt.^2  > 0.9 ) &  ...
	    abs(maxRhoSqLags - tlcShifts) > 0);

[i,j]=ind2sub(size(tlcOpt),ind);

Q = length(i);
for k=1:Q
  disp(sprintf('TF: %s  Target: %s  Shift: %d', ...
	       tfNames{i(k)}, ...
	       psGeneSymbols{j(k)}, ...
	       tlcShifts(i(k),j(k))))
end

