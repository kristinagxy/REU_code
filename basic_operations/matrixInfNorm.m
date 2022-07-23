function [anorm] = matrixInfNorm(A,block_size)
%
% compute inf-norm of a matrix with precision given by chop
%

if nargin < 2
    block_size = 256; % default block size
end

%first chop the matrix
A = abs(chop(A));
[m,~] = size(A);

s = zeros(m,1);

for i=1:m
    s(i)=vectorOneNorm(A(i,:)',block_size);
end

anorm = max(s);