function xnorm = vectorTwoNorm(x,block_size)
%
% compute 2-norm of a vector with precision given by chop
%
if nargin < 2
    block_size = 256; % default block size
end
a = abs(max(x));
xnorm = vv_blocked(x/a,x/a,block_size);
xnorm = chop(a * chop(sqrt(xnorm)));