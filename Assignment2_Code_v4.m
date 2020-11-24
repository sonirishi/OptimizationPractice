x = rand(1000,1);

px = x.^3 + x.^2 + 3*x + 1;

y = px ;

const = ones(1000,1);

A = cat(2,const,x,x.^2,x.^3);

cvx_begin
    variables b(4) u(1000)
    minimize (sum(u))
    subject to
        A*b - y >= -u;
        A*b - y <= u;
cvx_end
    