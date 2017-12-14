% K:  the number of independent time points for computing the correlation
%
% B:  the number of bins to use, for the kernel density estimation
%      of the probability
%
% sl:  smoothing length, in units of correlation coefficient
%
% st:  vector of time lags
%
% 
function [tlcPV, tlcShifts, toptPV, tlcOpt, tlcPVAll]=TimeLaggedCorrPValues( tlc, ...
						  tlcNoInt, ...
						  B, ...
						  st, ...
						  sl )
% T:  number of time lags
% M:  number of transcription factors
% N:  number of target genes
[T, M, N] = size(tlc);

tlcOpt = zeros(M,N);

% Nb:  number of genes for background
Nb = size(tlcNoInt, 2);

if T ~= size(tlcNoInt, 1)
  error 'invalid size for tlcNoInt (first dimension)';
end

if Nb ~= size(tlcNoInt, 3)
  error 'invalid size for tlcNoInt (third dimension)');
end

% build a vector of r values for the background distribution, at
% each time value
L=Nb*(Nb-1);

% time-lagged correlation values for background
tlcNoIntVec = zeros(T, L);
ctr = 1;
for i=1:Nb
  tlcNoIntVec(:, [ctr:(ctr + Nb - 2)])=squeeze(tlcNoInt(:, i, setxor(i, [1:Nb])));
  ctr = ctr + Nb - 1;
end

tlcPVAll = zeros(size(tlc));
tmpmat = zeros([M,N]);

tlcPVecNonTFAll = zeros(T, L);

for t=1:T
  % for each time lag

  % get the real TLCs for this time lag
  tlct = squeeze(tlc(t,:,:));
  
  % get the background TLCs for this time lag
  tlcNoIntT = tlcNoIntVec(t,:);

  % compute P values for the real TLCs
  tlcPVec=GaussKernDensCDF(tlcNoIntT.^2, ...
			  (tlct(:).^2)', ...
			  sl, ...
			  B, ...
			  0);

  % compute P values for the background TLCs
  tlcPVecNonTFAll(t,:)=GaussKernDensCDF(tlcNoIntT.^2, ...
				tlcNoIntT.^2, ...
				sl, ...
				B, ...
				0);
  tmpmat(:) = tlcPVec';
  tlcPVAll(t,:,:)=tmpmat(:,:);
end

% compute the background distribution of optimal TLCs
[tlcNonTFPVMin, tlcNonTFPVMinInd]=min(tlcPVecNonTFAll, [], 1);

% get the optimal lags for the background
tlcNonTFPVMinLags = st(tlcNonTFPVMinInd);

tlcNonTFShifts = st(tlcNonTFPVMinInd);
pthist = hist(tlcNonTFShifts(:), st);
pthist = pthist/sum(pthist);
stfull = [st, st+max(st)+st(2)-st(1)]
pthistsm = zeros(size(stfull));

L1=length(pthist);
L2=length(pthistsm);
for l=1:L2
  pthistsm(l) = sum(pthist.*exp(-(stfull(l) - st).^2 / (200))) / (sqrt(2*pi)*10);
end
pthistsm = pthistsm / sum(pthistsm);

pthistsm

Ptopt = pthistsm(1:L1);
k=8
theta=45/k;
Pint = 0.08;
PtoptInt = st.^(k-1) .* exp(-st/theta)/(theta^k * gamma(k));
PtoptInt = PtoptInt / sum(PtoptInt)
PNotIntTopt = 1 - (Pint*PtoptInt./Ptopt)

% get the optimal time lag P-values and indices
[tlcPVMin, tlcPVMinInd]=min(tlcPVAll, [], 1);

bckgdSigmas = -squeeze(log10(PNotIntTopt(tlcNonTFPVMinInd) .* tlcNonTFPVMin));
sigmas = -squeeze(log10( PNotIntTopt(tlcPVMinInd) .* tlcPVMin) );


tlcPVMinInd = squeeze(tlcPVMinInd);

for m=1:M
  for n=1:N
    tlcOpt(m,n)=tlc(tlcPVMinInd(m,n),m,n);
  end
end

tlcPV = GaussKernDensCDF(bckgdSigmas, ...
			 sigmas, ...
			 sl, ...
			 B, ...
			 0);

tlcShifts = st(tlcPVMinInd);
tlcNonTFPVMinInd = squeeze(tlcNonTFPVMinInd);


toptPV = squeeze(PNotIntTopt(tlcPVMinInd));

