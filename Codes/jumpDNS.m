function [theta_level_Alg, X_sample, W_sample, log_Z_final] = jumpDNS(which_action, log_like_sample, log_like_Alg, theta_level_Alg, t, m, k, C, visits_sum, lambda, second_phase, log_like_best, visits_max, log_Z_final, count, save_interval, beta_power)
%The function "jumpDNS" decides whether or not to jump a particle from its
%current level to the proposed level.

%%
%Compute X and W values
        
        X_sample = zeros(1, m);
        X_sample(1,1) = 1;
        
        if(m >= 2)

        for i = 2:m
     
            if(second_phase == 0)
            
            X_sample(1,i) = ((visits_sum(i-1,2) + C*exp(-1))/(visits_sum(i-1,1) + C))*X_sample(1,i-1);
            
            else
                
                X_sample(1,i) = ((visits_max(i-1,2) + C*exp(-1))/(visits_max(i-1,1) + C))*X_sample(1,i-1);
                
            end
     
        end
        
%Estimate log Z using Trapezoidal Rule

        X_diff = zeros(1,m - 1);
        
        for i = 2:m
            
            X_diff(1,i-1) = X_sample(1,i-1) - X_sample(1,i);
            
        end
        
        Z = 0;

        for i = 1:(m - 1)

            Z = Z + (exp(log_like_sample(1,i))+exp(log_like_sample(1,i+1)))*(X_diff(1,i)/2);
    
        end
        
        log_Z = log(Z);

%If Z estimate from Trapezoidal is zero, use the log scale estimate of
%Z.

if(Z == 0)

        log_X = log(X_sample);
        
        log_Z_vec = zeros(1, m - 1);
        
        for i = 2:(m - 1)
            
            if(X_sample(1,i) - X_sample(1,i+1) > 0)
            
                m_i = (log_like_sample(1,i)-log_like_sample(1,i+1))/(log_X(1,i) - log_X(1,i+1));
            
                b_i = log_like_sample(1,i+1) - m_i*log_X(1,i+1);
                
                if(m_i == -1)
                    
                    log_Z_vec(1, i-1) = b_i + log(log(abs(X_sample(1,i))) - log(abs(X_sample(1,i+1))));
                    
                end
                
                if(m_i > -1 || m_i < -1)
                                
                    log_Z_vec(1, i-1) = b_i + log((1/(m_i + 1))*((X_sample(1,i)^(m_i + 1)) - (X_sample(1,i+1)^(m_i + 1))));
                  
                end
            
            end
            
        end
        
        log_Z_vec(1, m - 1) = log_X(1,m) + log_like_best;
        
        %Use the LogSumExp (LSE) function, which is the logarithm of the 
        %sum of the exponentials of the arguments, to find the log evidence 
        %estimate Z
        log_Z = LogSumExp_func(log_Z_vec);
        
end
        
        if(second_phase == 0)
            
            log_Z_final = log_Z;
            
        else

            if(mod(count+2, save_interval) == 0)
            
            log_Z_final(1,(count+2)/save_interval) = log_Z;
            
            end
            
        end
        
        end
        
        if(second_phase == 0)
        
        W_sample = zeros(1,m);
        
        W_sample(1,1) = 1;
        
        if(m >= 2)
        
        for i = 1:m
 
            W_sample(1,i) = exp(((i-1) - (m-1))/(lambda));
                
        end
        
        end
        
        else
            
            W_sample = zeros(1,m);
            
            for i = 1:m
                
                W_sample(1,i) = 1;

            end
            
        end


%%

level_tminus1 = theta_level_Alg{k}(:,t(k)-1);

if(which_action < 0.5)
    
    log_like = log_like_Alg{k}(:,t(k)-1);
    
else
    
    log_like = log_like_Alg{k}(:,t(k));
    
end

if(m == 1)
    
    level_proposal = 0;
    
else

u = rand(1);
b = rand;

c = tan(pi*(b - 0.5));

mag = 1 + floor(c);

if(u < 0.5)
    
    level_proposal = level_tminus1 - 1*mag;
    
else
    
    level_proposal = level_tminus1 + 1*mag;
    
end

end

if(level_proposal >= m || level_proposal < 0 || level_proposal == level_tminus1)
    
    theta_level_Alg{k}(:,t(k)) = level_tminus1;
    
else
    
    if(log_like < log_like_sample(1,level_proposal+1))
        
        theta_level_Alg{k}(:,t(k)) = level_tminus1;
        
    else
        
        X_level_proposal = X_sample(1, level_proposal+1);
        W_level_proposal = W_sample(1, level_proposal+1);
        
        X_level_tminus1 = X_sample(1, level_tminus1+1);
        W_level_tminus1 = W_sample(1, level_tminus1+1);
        
        if(second_phase == 0)
        
        log_alpha = log(W_level_proposal) - log(W_level_tminus1);
        
        if(level_proposal < level_tminus1)
            
            log_alpha = log_alpha + log(X_level_tminus1) - log(X_level_proposal);
            
        end
        
        else
            
            log_visits_ratio = log(C + visits_max(level_tminus1+1,5)) - log(C + visits_max(level_proposal+1,5));
            
            if(log_visits_ratio < log(1.1) || log_visits_ratio > log(0.9))
                
            log_alpha = log(X_level_tminus1) - log(X_level_proposal);
            
            else
                
                log_alpha = (beta_power*log_visits_ratio) + (log(X_level_tminus1) - log(X_level_proposal));
                
            end
            
        end
        
        u2 = rand(1);
        
        if(log_alpha > 0)
            
            log_alpha = 0;
            
        end
        
        if(u2 < exp(log_alpha))
            
            theta_level_Alg{k}(:,t(k)) = level_proposal;
            
        else
            
            theta_level_Alg{k}(:,t(k)) = level_tminus1;
            
        end
        
    end
    
end

end

