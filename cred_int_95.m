function [ out ] = cred_int_95(k1, lelik)

%The function cred_int works provided that the posterior distribution of the parameter is symmetric

lower = 0;
m = 0;
upper = 0;

b = k1;

e = lelik;

[tmp, ind] = sort(b);
e_sort = e(ind);

CDF = e_sort;

for i = 1:length(tmp)
    
    if(i == 1)
        
        CDF(1,i) = CDF(1,i);
        
    end
    
    if(i ~= 1)
        
        CDF(1,i) = CDF(1,i) + CDF(1,i-1);
        
    end
    
end

CDF2 = CDF*(1/CDF(1,length(CDF)));

v = [];
j = 1;

for i = 1:length(CDF2)
    
    if(CDF2(1,i) >= 0.0230 && CDF2(1,i) <= 0.027)
        
        v(1, j) = i;
        
        j = j+1;
        
    end
    
end

lower = tmp(1,min(v));

v = [];
j = 1;

for i = 1:length(CDF2)
    
    if(CDF2(1,i) >= 0.480 && CDF2(1,i) <= 0.520)
        
        v(1, j) = i;
        
        j = j+1;
        
    end
    
end

m = tmp(1,min(v));

v = [];
j = 1;

for i = 1:length(CDF2)
    
    if(CDF2(1,i) >= 0.9730 && CDF2(1,i) <= 0.9770)
        
        v(1, j) = i;
        
        j = j+1;
        
    end
    
end

upper = tmp(1,max(v));

out = [lower m upper];

end

