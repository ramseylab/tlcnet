% TimeLaggedCorrList
%
% Compute optimal time-shifted correlation for a collection of
% vectors representing data for pairs of elements (x and y).
% 
% Stephen Ramsey
% Institute for Systems Biology
% May 2006
% 
% For a specified vector of time-shifts, compute the
% pearson correlation of x and a time-shifted y. The variables 
% x, y, and shiftTime should have the same number of rows.
%
% "shiftTime" is the time lag (must be positive) at which the
% time-lagged correlation is to be computed
%
% x is the putative regulator, and y is the putative target;
% this code depends on all rows of y being the *same*
%
% cols is a cell array of vectors.  Each vector within "cols"
% represents a time-course.  The elements of each vector are
% indices into the vector "x", indicating which column of "x"
% corresponds to the time value at the corresponding entry of the
% corresponding vector in the "times" cell array
%
% times is a cell array of vectors (same number of vectors, and
% same sizes for all vectors, as the cell array "cols").  Each
% vector within "times" represents a time-course.  the elements of
% each vector give the time values for the data points in the time
% course. 
function sc=TimeLaggedCorrList(shiftTime, ...
			       x, ...
			       y, ...
			       cols, ...
			       times)

[N,M]=size(x);

Q=length(times);
sc=zeros(N,1);

tkl=zeros(1,Q);
tk=[];
ctl=0;
numPts = 0;
for k=1:Q
  tk=times{k};
  tkl(k)=length(tk);
  ctl = ctl + tkl(k) - 1;
  numPts = numPts + 1 + length(find(tk - tk(1) > shiftTime));
end

% yboost is set of time-boosted values of y
yboost=zeros(1,ctl);
xi=zeros(1,ctl);
ck=[];

ctr = 1;
for k=1:Q
  ck=cols{k};
  tk = times{k};
  ti = tk(1:(tkl(k)-1)) + shiftTime;
  yboostk = interp1(tk, y(ck), ti, 'linear');
  yboost(ctr:(ctr+tkl(k)-2))=yboostk;
  ctr = ctr + tkl(k) - 1;
end
if length(find(isnan(yboost))) > 0
  error 'invalid yboost';
end

for i=1:N
  % for each row of x
  ctr = 1;
  
  for k=1:Q
    % get the list of columns
    ck=cols{k};

    % pull the xi time-series data out of the full array of x data
    xi(ctr:(ctr + tkl(k) - 2)) = x(i, ck(1:(tkl(k)-1)));

    % increment the counter
    ctr = ctr + tkl(k) - 1;
  end

  n=length(xi);
  cc = Covar(xi, yboost)/sqrt(Covar(xi, xi)*Covar(yboost, yboost));
  
  sc(i)=cc;

  if isnan(sc(i))
    error('invalid correlation');
  end
  
end


  
