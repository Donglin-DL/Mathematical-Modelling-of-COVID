function [ out ] = beta_like(par, time)

b = par(1);
q = par(6);
s = par(7);

time1beta = par(15);
time2beta = par(16);

beta_end = beta_t(time, b, q, s, time1beta, time2beta);

out = 0;

if(beta_end > b)
    
    out = -1000000;
    
end;

if(q > s)
    
    out = -1000000;
    
end;

end
