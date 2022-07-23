function p = vv_blocked(x,y,block_size)

if nargin < 3
    block_size = 256;
end

x = chop(x);
y = chop(y);
L = length(x);
k = floor(L/block_size); 

z = chop(x.*y);
p = 0;
for i = 1:k
    a=0;
     for j = (i-1)*block_size+1 : i*block_size
        a = chop(a+z(j));
    end
    p = chop (p + a);
end

if L-k*block_size ~= 0
    b=0;
    for i = k*block_size+1:L
        b = chop(b+z(i));
    end
    p = chop(p+b);
end