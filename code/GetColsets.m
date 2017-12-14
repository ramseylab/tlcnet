% defines colsets for 

function [cols,times]=GetColsets(strain, stim, mintime, maxtime)

[gotStrain, gotStim, gotTimept]=textread('UniqueExptNames.txt', ...
				       '%s %s %d', ...
				       'delimiter', '_');

M = size(gotStrain, 1);

if strcmp(strain, '*')==1
  strainInds = ones(M,1);
else 
  if strncmp(strain, '!', 1)==1
    filter = strain(2:length(strain));
    strainInds = (strcmp(gotStrain, filter)~=1);
  else
    strainInds = (strcmp(gotStrain, strain)==1);
  end
end

if strcmp(stim, '*')==1
  stimInds = ones(M,1);
else
  % stimulus filter is not a wildcard
  if strncmp(stim, '!', 1)==1
    filter = stim(2:length(stim));
    stimInds = (strcmp(gotStim, filter)~=1);
  else
    if length(strfind(stim, '|'))==0
      % stimulus filter is not a wildcard, and not a negated filter
      if nargin < 3 || (length(mintime)==1 && mintime > 0) || ...
	    (length(mintime)>1 && mintime(1)==0)
	stimInds = (strcmp(gotStim, stim)==1);
      else
	stimInds = (strcmp(gotStim, stim)==1) | (strcmp(gotStim, 'UNSTIM')==1);
      end	
    else
      k = strfind(stim, '|');
      filter1 = stim(1:(k-1));
      filter2 = stim((k+1):length(stim));
      stimInds = (strcmp(gotStim, filter1)==1) | (strcmp(gotStim, filter2)==1);
    end
  end
end

if nargin < 3
  timeptInds = ones(M,1);
else
  timeptInds = zeros(M,1);
  if length(mintime)==1
    inds = find(gotTimept >= mintime);
    if nargin > 3
      inds = intersect(inds, find(gotTimept <= maxtime));
    end
  else
    if length(mintime)>1
      inds = [];
      for j = 1:length(mintime)
	inds = cat(1, inds, find(mintime(j)==gotTimept));
      end
    else
      error 'invalid mintime, it is empty';
    end
  end
  timeptInds(inds) = 1;
end

allInds = strainInds & stimInds & timeptInds;

cols = find(allInds);

				       
if nargout > 1
  times = gotTimept(cols);
end