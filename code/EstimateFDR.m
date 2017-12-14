function pd=EstimateFDR(pv, fdr)

N = length(pv);
pvs = sort(pv);

j = [1:N];
d = max(find(pvs(:) < j(:)*fdr/N));
if length(d) == 1
  pd = pvs(d);
else
  error 'unable to estimate the P value cutoff';
end


