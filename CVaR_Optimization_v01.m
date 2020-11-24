tic;
loc = 'E:\Documents\IISc\Optimization\ProjectMaterial\Implementation\';

file = 'final_price_data.csv';

return_target = 0.009;

scenarios = 10000;

accuracy = 0.00001;

eta = 0.01;

beta = 0.05;

margin = 1;

[weight_matrix,qhat,CVaR,k] = portfolio_construction(loc,file,return_target,scenarios,accuracy,eta,beta,margin);

final_asset_allocation = weight_matrix(:,k+1);

sprintf("Optimal CVaR %d",CVaR)

toc;
function [weight_matrix,qhat,CVaR,counter] = portfolio_construction(loc,file,return_target,scenarios,accuracy,eta,beta,margin)
    
    rng(42);
    
    pricedata = readtable(strcat(loc,file));

    sizeof = size(pricedata);
    
    covariance_matrix = cov(table2array(pricedata(:,2:end)));
    
    mean_vector = mean(table2array(pricedata(:,2:end)));
    
    sim_return_matrix = zeros([scenarios,sizeof(2)-1]);
    
    for i = 1:sizeof(2)-1
        sim_return_matrix(:,i) = sqrt(covariance_matrix(i,i))*randn([scenarios,1]) + mean_vector(1,i);
    end
    
    D2 = log(1/beta)*(1/beta) - log(1/beta-1)*(1/beta-1);
    
    sigma2 = (1/beta)/((1/beta)-1);
    
    D1 = 1 + 2/margin;
    
    iterations = ceil((4/accuracy)*max(vecnorm(sim_return_matrix,2,2))*sqrt(D1*D2/sigma2));
    
    iterations = ceil(iterations/100); %To avoid memory issues
    
    sprintf("Maximum Iterations %d",iterations)
    
    weight_matrix = zeros(sizeof(2)-1,iterations+1);

    weight_matrix(:,1) = 1/(sizeof(2)-1);  %%%% Initializing by equal weighted portfolio
    
    mu = accuracy/(2*D2);
    
    prob = repelem(1/scenarios,scenarios);  % probability of each scenario is exactly same as its simulated
    
    q = zeros([scenarios,iterations]);

    qhat = zeros([scenarios,iterations]);
    
    A = [mean_vector;repelem(1,(sizeof(2)-1))];

    b = [return_target;1];

    L = eta + max(vecnorm(sim_return_matrix,2,2))/(sigma2*mu);

    y = zeros([sizeof(2)-1,iterations]);

    z = zeros([sizeof(2)-1,iterations]);

    for counter = 1:iterations
        
        sprintf("Iteration No %d Started",counter)
        
        alpha = wrap_alphasolve(prob,beta,weight_matrix,sim_return_matrix,mu,scenarios,counter);
        
        for i = 1:scenarios
            q(i,counter) = (prob(i)/beta)/(1+exp(((transpose(weight_matrix(:,counter))*transpose(sim_return_matrix(i,:)))/mu)+(alpha/mu)));
        end
        
        h = wrap_hsolve(A,weight_matrix,sim_return_matrix,b,L,eta,q,counter);
        
        y(:,counter) = (transpose(A)*h - (-transpose(sim_return_matrix)*q(:,counter) + eta*weight_matrix(:,counter)))*(1/L) + weight_matrix(:,counter);
        
        dum_vec = zeros([1,counter]);

        for i = 1:counter
            dum_vec(1,i) = i/2;
        end
        
        hhat = wrap_hhatsolve(A,weight_matrix,sim_return_matrix,b,L,eta,q,counter,dum_vec);
        
        z(:,counter) = (transpose(A)*hhat-(transpose(dum_vec*transpose(-transpose(sim_return_matrix)*q(:,1:counter)+eta*weight_matrix(:,1:counter)))))*(1/L);
        
        new_vec = zeros([1,counter]);

        for i = 1:counter
            new_vec(1,i) = (2*i)/((counter)*(counter+1));
        end

        qhat(:,counter) = sum((transpose(new_vec*transpose(q(:,1:counter)))),2);
        
        fmax = -transpose(q(:,counter))*sim_return_matrix*y(:,counter) + (eta/2)*sum_square(y(:,counter));
        
        gmin = -transpose(qhat(:,counter))*sim_return_matrix*weight_matrix(:,counter) + (eta/2)*sum_square(weight_matrix(:,counter));

        sprintf("Primal Dual Difference %d", fmax - gmin)
        
        if fmax - gmin >= accuracy
            weight_matrix(:,counter+1) = (2/(counter+3))*z(:,counter) + ((counter+1)/(counter+3))*y(:,counter);
            portret = sort(-sim_return_matrix*weight_matrix(:,counter+1));
            CVaR = mean(portret(1:beta*scenarios));
        else
            weight_matrix(:,counter+1) = y(:,counter);
            portret = sort(-sim_return_matrix*weight_matrix(:,counter+1));
            CVaR = mean(portret(1:beta*scenarios));
            sprintf("Stopping Here, Found the Optimal Solution")
            break;
        end
        
        sprintf("Iteration No %d Done",counter)
        
    end
end

function x1 = wrap_alphasolve(prob,beta,weight_matrix,sim_return_matrix,mu,scenarios,counter)
    fun = @alpha_solve;
    alpha0 = 0;
    x1 = fsolve(fun,alpha0);
    function F = alpha_solve(alpha)
        F = sum((prob(1:scenarios)./beta)./(1+exp(((transpose(weight_matrix(:,counter))*transpose(sim_return_matrix(1:scenarios,:)))./mu)+(alpha/mu)))) - 1;
    end
end

function h = wrap_hsolve(A,weight_matrix,sim_return_matrix,b,L,eta,q,counter)
    fun = @h_solve;
    h0 = randn([2,1]);
    h = fsolve(fun,h0);
    function F = h_solve(h)
        F = A*transpose(A)*h - A*(-transpose(sim_return_matrix)*q(:,counter) + eta*weight_matrix(:,counter)) - b*L + A*weight_matrix(:,counter)*L;
    end
end

function hhat = wrap_hhatsolve(A,weight_matrix,sim_return_matrix,b,L,eta,q,counter,dum_vec)
    fun = @hhat_solve;
    hhat0 = randn([2,1]);
    hhat = fsolve(fun,hhat0);
    function F = hhat_solve(hhat)
        F = A*transpose(A)*hhat-A*(transpose(dum_vec*transpose(-transpose(sim_return_matrix)*q(:,1:counter) + eta*weight_matrix(:,1:counter)))) - b*L;
    end
end

