function [out] = LogDiffExp_func(a, b)
%log(exp(a) - exp(b)), where a > b

out = a + log(1 - exp(b - a));

end

