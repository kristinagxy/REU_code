function [anorm] = matrixFroNorm(A,block_size)
%
% compute Frobenius-norm of a matrix with precision given by chop
%

if nargin < 2
    block_size = 256; % default block size
end

A = reshape(chop(A),1,[]); % flatten the matrix into a vector

anorm = vectorTwoNorm(A,block_size);