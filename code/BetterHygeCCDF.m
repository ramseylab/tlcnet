% BetterHygecdf
%
% Computes the complementary CDF of the hypergeometric distribution,
% but with better accuracy than the built-in MATLAB "hygecdf".
% Requires the statistics toolbox.
%
% Stephen Ramsey
% Institute for Systems Biology
% Aug 2006
function res=BetterHygeccdf(X, M, K, N)

Q = length(X);

res=zeros(size(X));

for i=1:Q
  xi = X(i);
  for j=xi:min(N(i),K(i))
    res(i)=res(i)+hygepdf(j, M(i), K(i), N(i));
  end
end

