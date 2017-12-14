% GaussKernDensCDF
%
% Stephen Ramsey, Institute for Systems Biology
% Mar. 2006
%
% Computes the complementary right-single-tailed CDF of the value x, from
% the estimated nonparametric distribution of the collection "b" of
% background measurements.  Larger "x" values will lead to a
% smaller p-value, but small "x" values will not lead to a small
% p-value (because we are computing only the right-single-tailed CDF).
%
% Function arguments are as follows:
%  b:    a row vector of background measurements from which the
%        nonparametric probability distribution is constructed
%  x:    a row vector of measurements for which the p-values are to
%        be computed, based on the nonparametric probability
%        density
%  smoothingLength:  the length scale (in the same units as "x" and
%        "b") for doing the Gaussian kernel density smoothing of 
%        the data points in "b".  Must be a positive-definite real
%        number.  Alternatively, a row vector of the same size as
%        "b", containing the smoothing length for each element of "b".
%  numBins:   the number of points at which to evaluate the
%        nonparametric PDF of the values in "b".  This should be
%        set to something reasonably large, like 100 or more.  It
%        cannot be set to a value less than 8.
%  minX:  an "optional" parameter indicating the minimum value for
%        any values in "x" and "b" (i.e., a value such that all 
%        "x" measurements will always be greater than this value).  
%        This parameter is typically used to specify that the "x" 
%        and "b" values are always positive-definite.  In this
%        case, the Gaussian smearing is performed with a
%        multiplicative correction for the amount of area under the
%        Gaussian that is less than (to the left of) the minX parameter.
%
% The p-values are returned as a vector.  Any resulting p-values
% greater than 1.0 (which might occur due to error in the numerical
% estimation of the integral, or roundoff error) are set to 1.0
% exactly.  The code is very "quick and dirty" but it is completely
% vectorized, and should scale well to large sizes for "x" and "b"
% vectors.  The vectors "x" and "b" do not have to have the same
% length.  The p-values are returned in the return vector "p".  The
% optional return cell array "d" is used to return the probability
% density.  The elements of "d" are as follows: d{1} is the set of
% density (PDF) values, and d{2} is the set of x values of the bin
% edges at which the probability distribution is estimated.  The
% numerical integration is performed using the Alternative Extended
% Simpson's Rule for most values in "x", but for values in "x" that
% are within the last 9 bins (highest values of "x"), the Extended
% Trapezoidal Rule is used.  Linear interpolation is used to estimate
% the integrand between bin edges.  Integration formulas used her
% are as described in "Numerical Recipes in C" (First Edition).

function [p,d]=GaussKernDensCDF(b, ...
				x, ...
				smoothingLength, ...
				numBins, ...
				minX)


if numBins <= 0
  error 'number of bins must be poisitive';
end

numX=length(x);
p=zeros(1,numX);

numB=length(b);

extensionFactor=5;
% get the minimum and maximum values of x and b

if length(smoothingLength)==1
  sig = repmat(smoothingLength,[1 numB]);
else
  if length(smoothingLength) ~= numB
    error ['the smoothing length, if passed as a vector, must have' ...
	   ' the same size as b'];
  else
    sig = smoothingLength;
  end
end

maxX=max([max(x) max(b)])+extensionFactor*sig(numB);

if nargin < 5
  minX=max([(min([min(x) min(b)])-extensionFactor*sig(1))]);
end

binSize =(maxX-minX)/(numBins-1);

binEdges=[minX:binSize:maxX];
probDist = zeros(1,numBins);
intProbDist = zeros(1,numBins);


% compute the nonparametric (smoothed) probability density function
% for the set of points "b", over the interval from minX to maxX
for j=1:numBins
  % need to correct for the fact that some of the area of each
  % Gaussian lies outside the [minX, maxX] interval
  areaInsideInterval = 1.0 - ...
      0.5*erfc( abs(maxX-b)./sig ) - ...
      0.5*erfc( abs(b-minX)./sig );
  
  probDist(j) = sum( (exp(-(b-binEdges(j)).^2./(2*sig.^2))) ./ ...
		              areaInsideInterval ./...
                (numB*sqrt(2*pi)*sig) );
end

intProbDist(numBins)=0.0;
intProbDist(numBins-1)=0.5*binSize*(probDist(numBins)+ ...
				    probDist(numBins-1));
intProbdist((numBins-7):(numBins-2))=0.5*binSize*...
                       (probDist((numBins-7):(numBins-2))+...
			probDist(numBins));
cs1=cumsum(probDist((numBins-1):-1:(numBins-6)));
intProbDist((numBins-7):(numBins-2))=intProbDist((numBins-7):(numBins-2))+...
                       binSize*cs1(6:-1:1);

cs2=cumsum(probDist((numBins-4):-1:5));
cs2((numBins-8):-1:1);

intProbDist(1:(numBins-8))=binSize*(...
    (17/48)*probDist(numBins) +...
    (59/48)*probDist(numBins-1) + ...
    (43/48)*probDist(numBins-2) + ...
    (49/48)*probDist(numBins-3) + ...
    (17/48)*probDist(1:(numBins-8)) + ...
    (59/48)*probDist(2:(numBins-7)) + ...
    (43/48)*probDist(3:(numBins-6)) + ...
    (49/48)*probDist(4:(numBins-5)) + ...
    cs2((numBins-8):-1:1));

% make sure that the first bin is exactly 1.0
intProbDist(1)=1.0;

% for each element in "x", compute the index of the element of
% binEdges that is just to the left of "x"
xInd=1+floor((x-minX)/binSize);

% for each element in "x", compute the fractional distance
% of "x" between the closest-left bin edge, and the closest-right
% bin edge
xSep=(x-binEdges(xInd))/binSize;

p = (xSep.*intProbDist(xInd+1) + (1-xSep).*intProbDist(xInd));

% for any p-values greater than 1.0, just set them to 1.0
p(find(p > 1.0))=1.0;

if nargout > 1
  d{1}=probDist;
  d{2}=binEdges;
end


