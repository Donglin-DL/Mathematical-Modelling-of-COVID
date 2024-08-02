function [sol] = model_all_out_tirebreaker2(par,n,Pop,min_rho)

beta1 = par(1);
f = par(2);
I0_1 = par(4);
gamma = par(5);
q = par(6);
s=par(7);

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

%sol holds values of S, I, cases, rho
sol = zeros(4,n);

%sol_incidence holds values of new referrals for testing

H = @(t,x) [-beta_t(t, beta1, q, s, time1beta, time2beta)*x(1)*x(2);
   beta_t(t, beta1, q, s, time1beta, time2beta)*x(1)*x(2)-(gamma)*x(2)-f*rhoFunction(t, time1, time1b, time2, time3, time4, height_val, const_val, min_rho, m0_val)*x(2);
   f*rhoFunction(t, time1, time1b, time2, time3, time4, height_val, const_val, min_rho, m0_val)*x(2)];

options1 = odeset('NonNegative', [1 2 3]);
[T, Y] = ode45(H, [0:1:n], [Pop-I0_1,I0_1,0], options1);

S = transpose(Y(2:(n+1),1));
I = transpose(Y(2:(n+1),2));
% f*rho*I cul case

sol(1,:) = S;
sol(2,:) = I;

X3 = Y(:,3);

for i = 2:(n+1)
    
    sol(3,i-1) = X3(i)-X3(i-1);
    
    
    %Ideally, these differences will always be positive because X3 are 
    %cumulative cases. However, since we are solving the ODEs
    %numerically, when these differences approach zero, it is possible to
    %get negative numbers that are very small in magnitude. So, we make 
    %these negative numbers very small positive numbers.
    if(sol(3,i-1) <= 0)
        
        sol(3, i-1) = 1e-15;
        
    end;
  
end;

sol(4,:) = rhoFunctionOut(time1, time1b, time2, time3, time4, height_val, const_val, min_rho, n, m0_val);

end

