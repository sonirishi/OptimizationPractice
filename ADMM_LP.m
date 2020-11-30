tic;
pricedata = readtable('E:\Documents\IISc\Optimization\ProjectMaterial\Implementation\final_price_data.csv');

rng(42);

scenarios=100;

return_target = 0.00045;

beta=0.05;

sizeof = size(pricedata);

covariance_matrix = cov(table2array(pricedata(:,2:end)));

mean_vector = mean(table2array(pricedata(:,2:end)));

return_matrix = zeros([scenarios,sizeof(2)-1]);
    
for i = 1:sizeof(2)-1
    return_matrix(:,i) = sqrt(covariance_matrix(i,i))*randn([scenarios,1]) + mean_vector(1,i);
end
    
mainc = 1/(scenarios*beta);
rho = 1;
A = zeros(scenarios+2,2*(scenarios)+sizeof(2)+1);
A(1,1:end) = [mean_vector repelem(0,scenarios+1) -1 repelem(0,scenarios)];
A(2,1:end) = [repelem(1,sizeof(2)-1) repelem(0,2*scenarios+2)];
A(3:end,1:end) = [return_matrix eye(scenarios,scenarios) ones(scenarios,1) zeros(scenarios,1) -1*eye(scenarios,scenarios)];
b = [return_target; 1; zeros(scenarios,1)];
c = transpose([repelem(0,sizeof(2)-1) repelem(mainc,scenarios) 1 repelem(0,1+scenarios)]);
lambda = ones(2*(scenarios)+sizeof(2)+1,1);
x = zeros(2*(scenarios)+sizeof(2)+1,1);
y = ones(2*(scenarios)+sizeof(2)+1,1);

for k = 1:1000
    xbar = (y - c/rho - lambda/rho);
    x(1:30) = xbar(1:30);  % no constraint on weights which are for 30 stocks
    for i = 31:30+scenarios  %% This the clipping for a defined in LP
        if xbar(i) < 0
            x(i) = 0;
        else
            x(i) = xbar(i);
        end
    end
    x(31+scenarios) = xbar(31+scenarios); % no constraint on tau which is the VaR value defined in LP
    for i = 32+scenarios:length(x) %% This the clipping for slacks
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

portret = sort(return_matrix*x(1:30));

cvar = mean(portret(1:beta*scenarios));

sprintf("Optimal CVaR %d",cvar)
toc;
