tic;
rho=0.1;
A = [110 205 160 160 420 260 -1 0 0; 4 32 13 8 4 14 0 -1 0; 2 12 54 285 22 80 0 0 -1];
b = [2000; 55;800];
c = [3 ;24 ;13 ;9 ;24 ;13 ;0 ;0 ;0];
lambda = ones(9,1);
x = ones(9,1);
y = zeros(9,1);
while norm(x(1:6)-y(1:6)) >= pow2(10,-14)
    xbar = (y - c/rho - lambda/rho);
    for i = 1:length(x)
        if xbar(i) < 0
            x(i) = 0;
        else
            x(i) = xbar(i);
        end
    end
    ybar = x + lambda/rho;
    mu = (A*transpose(A))\(b-A*ybar);
    y = ybar + transpose(A)*mu;
    lambda = lambda + rho*(x-y);
end      
toc;