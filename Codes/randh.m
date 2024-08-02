function [out] = randh()
%The function "randh" chooses a proposed value to perturb a parameter's 
%current position scaled by the difference in the parameter's lowest and 
%highest value set by its specified uniform distribution. (Please see line 
%48 in the "move" function.)

choose = rand;

a = randn;
n = randn;
b = rand;

if(choose < 1/2)
    
    t = a/sqrt(-log(b));
    
     out = (10^(1.5 - 3*abs(t)))*n;
    
else
    
    c = tan(pi*(b - 0.5));
    
    out = (10^(1 - abs(c)))*n;
    
end

end

