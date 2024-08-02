function [out] = wrapper_function(lower, upper, value)
%The "wrapper_function" ensures that poposed guesses of moving a
%parameter's value keeps the parameter within the bounds specified by the 
%parameter's uniform prior distribution.

out = lower + mod(value - lower, upper - lower);

end

