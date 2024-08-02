function [out] = penalty(par)

out = 0;

q = par(6);
s = par(7);

time1 = par(9);
time2 = par(10);
time3 = par(11);
time4 = par(12);

height_val = par(13);
const_val = par(14);

time1beta = par(15);
time2beta = par(16);

time1b =par(17);
m0_val =par(18);

if((q >= s) || ~(time1 < time1b && time1b < time2 && time2 < time3 && time3 < time4) || (height_val <= const_val) ||(height_val <= m0_val) || (time1beta >= time2beta))
    
    out = -Inf;
    
end

end

