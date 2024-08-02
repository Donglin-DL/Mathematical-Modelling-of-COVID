function [ out ] = beta_t(t, b, q, s, time1beta, time2beta)
%The time of intervention is March 15th to May 14th and open up after
%So, time 7 then time 30 remain restriction then time 67 open up

%time1beta and time2beta occur sometime between (29,32) and (65,68)

% 10% increase in slope: q2 = 0.1 +0.9*q1 where q1 is the original q;
%                        s2 = 1.1s1 where s1 is the original s
% 20%


if(t < 7)
    out = b;
end;


if(t>=7 && t <= time1beta)
   out = (-(1-q)*b/(time1beta-7))*(t-7)+b;
end;

if(t>=time1beta && t <= time2beta)
   out = q*b;
end;

if (t>= time2beta)
    out = s*b;

end;

