% BetterHygecdf
%
% Computes the CDF of the hypergeometric distribution, but with
% better accuracy than the built-in MATLAB "hygecdf".  Requires
% the statistics toolbox.
%
% Stephen Ramsey
% Institute for Systems Biology
% Aug 2006
function res=BetterHygecdf(X, M, K, N)

Q = length(X);

res=zeros(size(X));

for i=1:Q
  xi = X(i);
  for j=0:xi
    res(i)=res(i)+hygepdf(j, M(i), K(i), N(i));
  end
end

