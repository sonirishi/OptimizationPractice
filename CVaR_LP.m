tic;
pricedata = readtable('E:\Documents\IISc\Optimization\ProjectMaterial\Implementation\final_price_data.csv');

rng(42);

scenarios=10000;

return_target = 0.009;

beta=0.05;

sizeof = size(pricedata);

covariance_matrix = cov(table2array(pricedata(:,2:end)));

mean_vector = mean(table2array(pricedata(:,2:end)));

return_matrix = zeros([scenarios,sizeof(2)-1]);
    
for i = 1:sizeof(2)-1
    return_matrix(:,i) = sqrt(covariance_matrix(i,i))*randn([scenarios,1]) + mean_vector(1,i);
end
    
A = [mean_vector;repelem(1,(sizeof(2)-1))];

b = [return_target;1];

prob = repelem(1/scenarios,scenarios);  % probability of each scenario is exactly same as its simulated

cvx_begin
    variables tau a(scenarios,1) w(30,1)
    minimize tau + (1/(scenarios*beta))*sum(a)
    subject to
        A*w==b;
        a >= 0;
        return_matrix*w + tau + a >= 0;
cvx_end

portret = sort(-return_matrix*w);

cvar = mean(portret(1:beta*scenarios));

sprintf("Optimal CVaR %d",cvar)
toc;
