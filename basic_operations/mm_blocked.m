function C = mm_blocked(X,Y,block_size)
%
% copmute the product of two matrices with chop
%

if nargin < 3
    block_size = 256; % default block size
end

[m, n] = size(X);
[~, p] = size(Y);

k = floor(n/block_size);
C = zeros(m,p);

for i = 1:k
    a=zeros(m,p);
    for j = (i-1)*block_size+1 : i*block_size
        a = chop(a+chop(chop(X(:,j))*chop(Y(j,:))));
    end
    C = chop(C + a);
end

if n-k*block_size ~= 0
    b=zeros(m,p);
    for i = k*block_size+1:n
        b = chop(b+ chop(chop(X(:,i))*chop(Y(i,:))));
    end
    C = chop(C + b);
end