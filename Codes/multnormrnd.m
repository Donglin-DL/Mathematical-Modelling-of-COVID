function [ out ] = multnormrnd(mean, sigma, rep, t)

out = zeros(t, rep);

for i = 1:t

    for j = 1:rep
    
        out(i,j) = normrnd(mean(1,i), sigma);
    
    end;

end;

end