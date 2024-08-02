function [ out ] = multnbinrnd(mean, p, t)
%The function "multnbinrnd" returns random samples from a negative binomial
%distribution with mean "mean" and dispersion parameter "p" across the time 
%t.

out = zeros(t, 1);

for i = 1:t
    
    out(i,1) = nbinrnd((mean(i,1)*p)/(1 - p), p);
    
end


end

