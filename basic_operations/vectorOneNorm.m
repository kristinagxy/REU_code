function [xnorm] = vectorOneNorm(x,block_size)
%
% compute one-norm of a vector with precision given by chop
%

if nargin < 2
    block_size = 256; % default block size
end

x = abs(chop(x));
l = length(x);
k = floor(l/block_size);

xnorm = 0;
for i = 1:k
    a=0;
    for j = (i-1)*block_size+1 : i*block_size
        a = chop(a+x(j));
    end
    xnorm = chop (xnorm + a);
end

if l-k*block_size ~= 0
    b=0;
    for i = k*block_size+1:l
        b = chop(b+ x(i));
    end
    xnorm = chop(xnorm+b);
end




