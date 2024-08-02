function [out] = LogSumExp_func(vec)

        max_vec = max(vec);
        
        a = vec - max_vec;
        
        out = log(sum(exp(a))) + max_vec;

end

