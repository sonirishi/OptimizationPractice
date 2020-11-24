x = rand(1000,1);

px = x.^3 + x.^2 + 3*x + 1;

y = px + randn(1000,1);

const = ones(1000,1);

A = cat(2,const,x,x.^2,x.^3);

cvx_begin
    variables b(4) u(1000,2)
    minimize (sum(u(:)))
    subject to
        A*b - y >= sum(u,2);
        A*b - y <= sum(u,2);
        sum(u,2) >= 0;
cvx_end
    