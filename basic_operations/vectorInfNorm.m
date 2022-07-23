function [xnorm] = vectorInfNorm(x)
%
% compute inf-norm of a vector with precision given by chop
%

%first chop the vector
x = abs(chop(x));
xnorm = max(x);