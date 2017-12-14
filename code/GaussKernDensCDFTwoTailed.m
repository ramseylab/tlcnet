% GaussKernDensCDFTwoTailed
%
% Stephen Ramsey, Institute for Systems Biology
% Jun 2006
%
function p=GaussKernDensCDFTwoTailed(b, ...
				     x, ...
				     sl, ...
				     numBins, ...
				     center)

if nargin > 4
  medb = center;
else
  medb = median(b);
end

ind = find(b <= medb);
indx = find(x <= medb);

p = zeros(size(x));

if length(sl==1)
  sluse=sl*ones(size(b));
else
  sluse=sl;
end

p(indx)=GaussKernDensCDF(medb-b(ind), medb-x(indx), sluse(ind), numBins, ...
			 0);
indx = find(x >= medb);
ind = find(b >= medb);
p(indx)=GaussKernDensCDF(b(ind)-medb, x(indx)-medb, sluse(ind), numBins, ...
			 0);




