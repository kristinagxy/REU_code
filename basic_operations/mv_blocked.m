function C = mv_blocked(X,y,block_size)
%
% compute the product of a matrix and a vector with chop
%

A = chop(chop(X).*chop(y'));

if nargin < 3
    block_size = 256; % default block size
end

[m, n] = size(X);

k = floor(n/block_size);
C = zeros(m,1);

for i = 1:k
    a=zeros(m,1);
    for j = (i-1)*block_size+1 : i*block_size
        a = chop(a+A(:,j));
    end
    C = chop(C + a);
end

if n-k*block_size ~= 0
    b=zeros(m,1);
    for i = k*block_size+1:n
        b = chop(b+ A(:,i));
    end
    C = chop(C + b);
end