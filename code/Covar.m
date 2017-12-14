% Covar
%
% Stephen Ramsey, Institute for Systems Biology
% Jun. 2006
%
% Computes the covariance for a pair of random variables.
function r = Covar( x, y )

n = length(x);
meanx = mean(x);
meany = mean(y);
r=mean(((x-meanx).*(y-meany)));

