function [ out ] = R_like(par, S, sol_rho)

beta1 = par(1);
p1 = par(2);
%I0_1 = par(4);
gamma = par(5);
q=par(6);
s=par(7);

%R = (beta1*(Pop-I0_1))/(gamma + p1*min(rhoFunctionOut(par(7), par(8), par(9), par(10), 400)));



%R = (beta1*(Pop-I0_1))./(gamma + p1*pTime(data_rho,0));

R = (beta1*(S(1:16)))./(gamma + p1*sol_rho(1:16)); %beta1 is constant for t<16 instead of t=48 (April 25)

out = 0;

if(min(R) <= 1 || max(R) > 2) %consider changing to 3
    
    out = -1000000;
    
end;


end

