% StudentCDF.m
%
% Stephen Ramsey, Institute for Systems Biology
% Feb. 2006
%
% This is a vectorized function for computing the CDF of Student's 
% t distribution.  The code is based on the "student_cdf.m" program
% included in the "PROB" MATLAB library written by John Burkardt
% (Florida State University).  It uses the incomplete Beta function
% that is built-in to MATLAB ("betainc"), and thus bypasses the
% need to use the MATLAB Statistics Toolbox.  Please note that only
% vectorized x, a, and b are allowed.  
%  x:  the mean value of the reference distribution
%  a:  the mean value of the data we are comparing to the reference
%      distribution (the null hypothesis is that the observations
%      whose mean value is "a", come from the reference
%      distribution)
%  b:  the uncertainty (standard deviation) of the reference
%      distribution
%  c:  the number of degrees of freedom (n-1).  Please note that
%  "c" must be a scalar!  Also, the parameters "x", "a", and "b"
% must all have the exact same length and dimensions.

function cdf=StudentCDF(x, a, b, c)

y=(x-a)./b;
a2=0.5*c;
b2=0.5;
c2=c ./(c + y.*y);
cdf = betainc(c2, a2, b2);

% ----------- uses the Statistics toolbox --------------
% t = (x-a)./b;
% cv = c*ones(size(x));
% p = tcdf(t, cv);
% indp = find(p > 0.5);
% p(indp)=1-p(indp);
% cdf=2.0*p;
% ----------- uses the Statistics toolbox --------------
