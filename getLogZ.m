function [log_Z] = getLogZ(log_like_sample, log_like_best, visits_sum, C)
%The function "getZ" returns the log evidence estimate Z given the cut 
%off values of the log likelihood, currently highest log likelihood value, 
%visits matrix, and the number C, which indicates the amount of confidence 
%in the theoretical expectation.

        X_sample = zeros(1, length(log_like_sample));
        X_sample(1,1) = 1;

        for i = 2:length(log_like_sample)

            X_sample(1,i) = ((visits_sum(i-1,2) + C*exp(-1))/(visits_sum(i-1,1) + C))*X_sample(1,i-1);
     
        end
        
        log_X = log(X_sample);
        
        log_Z_vec = zeros(1, length(log_like_sample) - 1);
        
        for i = 2:(length(log_like_sample) - 1)
            
            if(X_sample(1,i) - X_sample(1,i+1) > 0)
            
                m1 = (log_like_sample(1,i)-log_like_sample(1,i+1))/(log_X(1,i) - log_X(1,i+1));
            
                n1 = log_like_sample(1,i+1) - m1*log_X(1,i+1);
                
                if(m1 == -1)
                    
                    log_Z_vec(1, i-1) = n1 + log(log(abs(X_sample(1,i))) - log(abs(X_sample(1,i+1))));
                    
                end
                
                if(m1 > -1 || m1 < -1)
                                
                    log_Z_vec(1, i-1) = n1 + log((1/(m1 + 1))*((X_sample(1,i)^(m1 + 1)) - (X_sample(1,i+1)^(m1 + 1))));
                  
                end
            
            end
            
        end
        
        log_Z_vec(1, length(log_like_sample) - 1) = log_X(1,length(log_like_sample)) + log_like_best;
        
        %Use the LogSumExp (LSE) function, which is the logarithm of the 
        %sum of the exponentials of the arguments, to find the log evidence 
        %estimate Z
        
        max_log_Z_vec = max(log_Z_vec);
        
        a = log_Z_vec - max_log_Z_vec;
        
        log_Z = log(sum(exp(a))) + max_log_Z_vec;

end

