function [X,info,ProbInfo,lambda] = run_cs(tomo,size,noise,blur,precision,tol)
%
% tomo: true->tomo, false->blur
% size: 64,32
% noise: 0, 0.001, 0.01, 0.1
% blur: true->mild, false->default
% precision: double, single, fp16
% tol:tolerance

options.format = precision;
chop([],options);

prob_options = PRset('BlurLevel',[]);
if blur
    prob_options = PRset('BlurLevel','mild');
end

[A, b, xtrue, ProbInfo] = PRblur(size,prob_options);

if tomo
    [A, b, xtrue, ProbInfo] = PRtomo(size,prob_options);
end

if noise~=0
    [b, ~] = PRnoise(b, noise);   
end

cgls_options = IRset('x_true', xtrue);

if strcmp(precision,'double')
    [X,info,lambda] = cs_o(full(A),b,tol,1:500,cgls_options);
else
    [X,info,lambda] = cs_chop(full(A),b,tol,1:500,cgls_options);
end