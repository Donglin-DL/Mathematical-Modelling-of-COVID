function [ rho_val ] = rhoFunction(t, time1, time1b, time2, time3, time4, height_val, const_val, min_rho, m0_val)

t1 = floor(time1);

t1b = floor(time1b);

t2 = floor(time2);

t3 = floor(time3);

t4 = ceil(time4);

if(t <= t1)
    
    rho_val = min_rho;
    
end

if(t > t1 && t <= t1b)
    
    m0 = (m0_val - min_rho)/(t1b - t1);
    
    rho_val = m0*(t - t1) + min_rho;
    
end

if(t > t1b && t <= t2)
    
    m1 = (height_val - m0_val)/(t2 - t1b);
    
    rho_val = m1*(t - t1b) + m0_val;
    
end

if(t > t2 && t <= t3)
    
    rho_val = height_val;
    
end

if(t > t3 && t < t4)
    
    m2 = (const_val - height_val)/(t4 - t3);
    
    rho_val  = m2*(t - t3) + height_val;
    
end

if(t >= t4)
    
    rho_val = const_val;
    
end

end

