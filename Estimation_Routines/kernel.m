function res = kernel(dist,bandwidth)
% dist =j <= q
% bandwidh =q+1 if Bartlett kernel.
xx = 1 - abs(dist/(bandwidth));

if xx>=0
    res = xx;
else
    res = 0;
end
%
%res = exp(-(dist/bandwidth)^2);