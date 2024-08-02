function [ out ] = multnormrnd_sigma_t(mean, sigma, rep, t)

out = zeros(t, rep);

for i = 1:t

    for j = 1:rep
    
        out(i,j) = normrnd(mean(1,i), sigma(1,i));
    
    end;

end;

end

