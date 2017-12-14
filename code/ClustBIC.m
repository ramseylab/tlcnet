function [bestKidx, bestKC, kvals, kBIC]=ClustBIC(X, krange, numRep, varK)
%ClustBIC - k-means clustering, using Shwarz Criterion to select K
%
% Perform a k-means clustering for a range of K values specified by
% the user.  For each K value, the optimality of the clustering is
% computed using the Schwarz Criterion (which is equivalent to the
% BIC, up to a constant factor).  The K value at which the
% clustering is most optimal, is selected and the clustering
% results for this K value are returned.  The formula for computing
% the Schwarz Criterion is taken from pp. 206-207 of the book "The 
% Elements of Statistical Learning", by Trevor Hastie et al. (New York,
% Springer, 2001), as well as Andrew Moore's online tutorial on
% k-means clustering at Carnegie Mellon University
% (http://www.autonlab.org/tutorials/kmeans.html, October 19,
% 2006).  This clustering code does not allow for "outliers" that
% are excluded from the clustering (when the code calls the
% underlying k-means clustering algorithm, the fraction of
% allowed outliers is intentionally set to zero).  This code uses
% the k-means implementation from Piotr Dollar's MATLAB toolbox.
%
% Syntax:  [bestKidx, bestKC] = ClustBIC(X, min
%
% Inputs:
%    X - a NxM matrix of data to be clustered.  N is the number of
%    features, and M is the "dimensionality of the space" (the
%    number of measurements made, for each feature).
%
%    krange - row vector containing the possible numbers of clusters
%
%    numRep - the number of replicates for running the k-means
%    algorithm, with random initial cluster centers.  The bigger
%    this number is, the slower the overall code runs (but
%    asymptotically, the closer to truly optimal the clustering
%    results will be).
% 
%    varK - (optional) the K value at which the variance in
%    X-clust(X) is to be estimated (otherwise, min(krange) is 
%    used)
%
% Outputs:
%    bestKidx - the N vector indicating the cluster to which each
%    of the N features is assigned, for the clustering based on the
%    optimal number of clusters, Kbest.
%
%    bestKC - the KxM matrix indicating the cluster centers, for
%    the clustering based on the optimal number of clusters,
%    "Kbest". Therefore, Kbest = size(bestKC, 1).
%
% Example: 
%    
%
% Other m-files required: kmeans2, in Piotr Dollar's Video & Image
% Analysis Toolbox (http://vision.ucsd.edu/~pdollar/toolbox/doc).
%
% Subfunctions: none
%
% MAT-files required: none
%
% See also: 

% Author: Stephen A. Ramsey
% Institute for Systems Biology
% 1441 N 34th St
% Seattle, WA 98103 USA
% email: sramsey@systemsbiology.org
% Website: http://magnet.systemsbiology.net/publicRepository/sramsey
% Oct. 2006; Last revision: 2006/10/19

[N,M]=size(X);

bestKbic = inf;
bestKidx = [];
bestKC = [];
bestK = 0;

meanSqError = -1;

if nargin < 4
  varK = min(krange);
end
  
numK = length(krange);
bicSave = zeros(1,numK);
kvals = zeros(1,numK);

if length(intersect(varK, krange)) == 0
  krange = [varK, krange];
  numk = numK + 1;
else
  krange = [varK, setxor(varK, krange)];
end

kctr = 1;
for k=krange
  kvals(kctr)=k;
  k
  [idx, C, sumd]=kmeans(X, k, ...
		  'replicates', numRep, ...
		  'emptyaction', 'singleton');

%  [idx, C, sumd]=kmeans2(X, k, ...
%			 'replicates', numRep, ...
%			 'outlierfrac', 0);
  
  % C is a KxM matrix of cluster centers
  [idxsort, idxsorti]=sort(idx);

  % sort the matrix of values in cluster order
  Xsort = X(idxsorti,:);
  
  % build a NxM matrix of cluster centers
  XsortModel = zeros(N,M);

  kused = size(C,1);
  
%   ctr = 1;
%   for j=1:kused
%     numMembers = length(find(idx == j));
%     XsortModel([ctr:(ctr + numMembers - 1)],:)=...
% 	repmat(C(j,:),numMembers,1);
%     ctr = ctr + numMembers;
%   end
  
  % compute the distortion term  
  D = sum(sumd);

  % compute the mean-squared error (if first iteration)
  if meanSqError < 0
    meanSqError = D / (N*M);
  end
  
  % normalize the distortion by the mean-squared error
  D = D / meanSqError;
  
  % compute the penalty term
  bic = D + M*k*log(N)
  
  % if this BIC is the best we have observed yet, save the k and
  % the clustering results
  if bic < bestKbic
    bestK = kused;
    bestKbic = bic;
    bestKidx = idx;
    bestKC = C;
  end
  
  bicSave(kctr)=bic;
  kctr = kctr + 1;
end

if nargout > 2
  [ksrt, ksrtind] = sort(kvals);
  kvals = ksrt;
  kBIC = bicSave(ksrtind);
end

