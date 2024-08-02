function [R_lower, R_upper] = HPD(parameter_MCMC, alpha)
%The function "HPD" returns the upper and lower bound of the Highest 
%Posterior Density (HPD) for a parameter based on the Chen-Shao Highest 
%Posterior Density (HPD) Estimation Algorithm
%(ex. the smallest 95% credible interval will be given by the HPD using 
%alpha=0.05)

parameter_MCMC_ordered = sort(parameter_MCMC);

n = length(parameter_MCMC);

R_j = zeros(2, n - fix((1 - alpha)*n));

for j = 1:(n - fix((1 - alpha)*n))
    
    R_j(1,j) = parameter_MCMC_ordered(1, j);
    R_j(2,j) = parameter_MCMC_ordered(1, j + fix((1 - alpha)*n));
    
end

R_j_diff = zeros(1, n - fix((1 - alpha)*n));

for j = 1:(n - fix((1 - alpha)*n))
    
    R_j_diff(1, j) = R_j(2,j) - R_j(1,j);
    
end

[~, ind] = min(R_j_diff);

R_lower = R_j(1,ind);

R_upper = R_j(2,ind);

end

