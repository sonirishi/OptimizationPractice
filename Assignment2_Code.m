x = rand(1000,1);

px = x.^3 + x.^2 + 3*x + 1;

y = px + randn(1000,1);

const = ones(1000,1);

A = cat(2,const,x);

cvx_begin
    variable b(2)
    minimize (transpose(A*b - y)*(A*b-y))
cvx_end
    