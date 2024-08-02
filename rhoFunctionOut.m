function [ v ] = rhoFunctionOut(time1, time1b, time2, time3, time4, height_val, const_val, min_rho, max_time, m0_val)

v = repmat(min_rho,1,max_time);

t1 = floor(time1);

t1b = floor(time1b);

t2 = floor(time2);

t3 = floor(time3);

t4 = ceil(time4);

for i = (t1+1):t1b
    
    m0 = (m0_val - min_rho)/(t1b - t1);
    
    v(1, i) = m0*(i - t1) + min_rho;
    
end;

for i = (t1b+1):t2
    
    m1 = (height_val - m0_val)/(t2 - t1b);
    
    v(1, i) = m1*(i - t1b) + m0_val;
    
end;

for i = (t2+1):t3
    
    v(1, i) = height_val;
    
end;

for i = (t3+1):(t4-1)
    
    m2 = (const_val - height_val)/(t4 - t3);
    
    v(1, i) = m2*(i - t3) + height_val;
    
end;

for i = t4:max_time
    
    v(1, i) = const_val;
    
end;

end

