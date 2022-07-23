function [anorm] = matrixOneNorm(A,block_size)
%
% compute one-norm of a matrix with precision given by chop
%

if nargin < 2
    block_size = 256; % default block size
end

%first chop the matrix
A = abs(chop(A));
[~,n] = size(A);

s = zeros(n,1);

for i=1:n
    s(i)=vectorOneNorm(A(:,i),block_size);
end

anorm = max(s);