function FindMedianGamma(theta, k)

fminsearch(@(x) cdfMinusZeroPointFive(x, theta, k), 1)

function f = cdfMinusZeroPointFive(x, theta, k)
f = abs(0.5 - gamcdf(x, k, theta));
return;
