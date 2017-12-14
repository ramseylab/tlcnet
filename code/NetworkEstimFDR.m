function NetworkEstimFDR

load 'ClusterScans';
load 'TimeLaggedCorr';

pcomb = gammainc(-log(pcsPV .* tlcPVC), 2);
fdr=EstimateFDR(pcomb(find(pcomb < 1))', 0.03)

fracFiltered = 106/length(find(pcomb < 1e-2));

fdr / fracFiltered

