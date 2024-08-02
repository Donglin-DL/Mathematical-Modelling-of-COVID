function [theta_Alg, log_prior_Alg, log_like_Alg, my_changed_theta, log_like_best_changed, par_best_changed] = move(which_action, log_like_sample, theta_level_Alg, theta_Alg, log_prior_Alg, log_like_Alg, logfuns, t, k, log_like_best, widths, n1, par_best)
%The function "move" decides whether or not to move a particle from its
%current position to the proposed position within its current level.

        log_like_best_changed = log_like_best;
        par_best_changed = par_best;

        my_changed_theta = 0;
            
        if(which_action < 0.5)
            
            log_L_s = log_like_sample(1,1+theta_level_Alg{k}(:,t(k)));
            
        else
            
            log_L_s = log_like_sample(1,1+theta_level_Alg{k}(:,t(k)-1));
            
        end
        
        theta_k_tminus1 = theta_Alg{k}(:, t(k)-1);
        
        log_prior_k_tminus1 = log_prior_Alg{k}(1,t(k)-1);
        log_like_k_tminus1 = log_like_Alg{k}(1,t(k)-1);
        
        amount_to_move = 1;
        
        amount_to_move_u = rand(1);
        
        if(amount_to_move_u <= 0.5)
            
            amount_to_move = round(n1^(amount_to_move_u));
            
        end
        
        Y = zeros(1,n1);
        
        for i = 1:n1

                Y(1,i) = theta_k_tminus1(i);
                
        end
        
        %choose parameter(s) to move
        par_move = datasample(1:n1,amount_to_move,'Replace',false);
        
        for i = 1:amount_to_move
            
            temp_value = theta_k_tminus1(par_move(1,i)) + (widths(par_move(1,i),2) - widths(par_move(1,i),1))*randh();
            
            Y(1,par_move(1,i)) = wrapper_function(widths(par_move(1,i),1), widths(par_move(1,i),2), temp_value);
           
        end
        
        prior_theta_g = logfuns{1}(Y);
        
        if(isinf(prior_theta_g))
            
            theta_Alg{k}(:,t(k)) = theta_k_tminus1;
            log_prior_Alg{k}(1,t(k)) = log_prior_k_tminus1;
            log_like_Alg{k}(1,t(k)) = log_like_k_tminus1;
            
        else
            
            alpha = min([1 exp(prior_theta_g - log_prior_k_tminus1)]);
             
            u = rand(1);
            
            if(u < alpha)
            
            like_theta_g = logfuns{2}(Y);
            
            if(like_theta_g <= log_L_s)
            
                theta_Alg{k}(:,t(k)) = theta_k_tminus1;
                log_prior_Alg{k}(1,t(k)) = log_prior_k_tminus1;
                log_like_Alg{k}(1,t(k)) = log_like_k_tminus1;
                
            else
                
                my_changed_theta = 1;
                
                theta_Alg{k}(:,t(k)) = Y;
                log_prior_Alg{k}(1,t(k)) = prior_theta_g;
                log_like_Alg{k}(1,t(k)) = like_theta_g;
                
                if(like_theta_g > log_like_best)
                    
                    log_like_best_changed = like_theta_g;
                    par_best_changed = Y;
                    
                end
            
            end
            
            else
                
                theta_Alg{k}(:,t(k)) = theta_k_tminus1;
                log_prior_Alg{k}(1,t(k)) = log_prior_k_tminus1;
                log_like_Alg{k}(1,t(k)) = log_like_k_tminus1;
                
            end
        
        end

end

