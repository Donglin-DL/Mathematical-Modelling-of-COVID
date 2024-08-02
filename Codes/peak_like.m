function [ out ] = peak_like(par, max_time)

peakTime = par(8);

t2 = floor(peakTime);

out = 0;

if(~(t2 == max_time))
    
    out = -1000000;
    
end;


end