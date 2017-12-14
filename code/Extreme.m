% returns the most extremal value of each column of x, which may be
% either the max() or the min() depending on the whether
% the most extremal value is positive or negative. 
% 

function ex=Extreme(x)

if isvector(x)
  [mx,mi]=max(abs(x));
  ex = sign(x(mi))*mx;
else
  [N,M]=size(x);
  [mx, mi]=max(abs(x));
  j = [1:M];
  ind = sub2ind([N,M],mi,j);
  ex = mx .* sign(x(ind));
end


