function [logLike] = LogLikeFunc_tirebreaker2(par,n, Pop, data_case, data_rho, min_rho)

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

sol = zeros(1, n);

H = @(t,x) [-beta_t(t, beta1, q, s, time1beta, time2beta)*x(1)*x(2);
   beta_t(t, beta1, q, s, time1beta, time2beta)*x(1)*x(2)-(gamma)*x(2)-f*rhoFunction(t, time1, time1b, time2, time3, time4, height_val, const_val, min_rho, m0_val)*x(2);
   f*rhoFunction(t, time1, time1b, time2, time3, time4, height_val, const_val, min_rho, m0_val)*x(2)];

options1 = odeset('NonNegative', [1 2 3]);
[T, Y] = ode45(H, [0:1:n], [Pop-I0_1,I0_1,0], options1);


S = transpose(Y(2:(n+1),1));
I = transpose(Y(2:(n+1),2));

X3 = Y(:,3);

for i = 2:(n+1)
    
    sol(1,i-1) = X3(i)-X3(i-1);
    
    
    %Ideally, these differences will always be positive because X3 are 
    %cumulative cases. However, since we are solving the ODEs
    %numerically, when these differences approach zero, it is possible to
    %get negative numbers that are very small in magnitude. So, we make 
    %these negative numbers very small positive numbers.
    if(sol(1,i-1) <= 0)
        
        sol(1, i-1) = 1e-15;
        
    end;
  
end;

sol_rho = rhoFunctionOut(time1, time1b, time2, time3, time4, height_val, const_val, min_rho, length(data_rho), m0_val);

data_all = [data_case; data_rho];
sol_all = [sol; sol_rho];

sol_rho_T = sol_rho';
sol_T = sol';

% logLike = log(prod((par(3)*(sol_rho_T.^(-1))).^(1/2))) + log(prod((par(8)*(sol_rho_T.^(-1))).^(1/2))) - (1/2)*trace(((data_all - sol_all).^2)*[par(3)*(sol_rho_T.^(-1)) par(8)*(sol_rho_T.^(-1))]) + R_like(par, S, sol_rho) + beta_like(par, length(data_case));

logLike = log(prod(((par(3)*sol_T).^(-1)).^(1/2))) + log(prod(((par(8)*sol_rho_T).^(-1)).^(1/2))) - (1/2)*trace(((data_all - sol_all).^2)*[(par(3)*sol_T).^(-1) (par(8)*sol_rho_T).^(-1)]) + R_like(par, S, sol_rho) + beta_like(par, length(data_case));


end

