function [X,info,ProbInfo] = run_cgls(tomo,size,noise,blur,precision)
%
% tomo: true->tomo, false->blur
% size: 64,32
% noise: 0,0.01, 0.1, 0.001
% blur: true->mild, false->default
% precision: double, single, fp16
%

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
    [X,info] = IRcgls(full(A),b,1:100, cgls_options);
else
    [X,info] = IRcgls_chop(full(A),b,1:100, cgls_options);
end