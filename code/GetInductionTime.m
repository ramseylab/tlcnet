function indTime=GetInductionTime(times, vals, thresh, maxtime)

minval = min(abs(vals));
maxval = max(abs(vals));
if thresh > maxval
  error 'value never gets big enough';
end
if thresh < minval
  error 'value never gets small enough';
end

mintime = min(times);
if nargin < 4
  maxtime = max(times);
end
tbest = -1;
fbest = inf;
trange = mintime:1:maxtime;
q = thresh - abs(interp1(times, vals, trange));
sq = sign(q);
[sqf, sqfi]=find( sq ~= sq(2) );

msqfi = min(sqfi);
indTime = 0.5*(trange(msqfi-1)+trange(msqfi));
