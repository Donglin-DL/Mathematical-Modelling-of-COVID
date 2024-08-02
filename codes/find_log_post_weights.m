function [l_post_weight_sort, l_X_sort, index_sort, l_p_w_sort_c] = find_log_post_weights(l_like_save, log_like_star, X_sample_save_median, level_assign)
%The function "find_log_post_weights" estimates the log prior mass and the 
%log posterior weight for each sample from the DNS output (output is sorted 
%by increasing likelihood value), and this information is used to obtain 
%the representative log posterior sample

n = size(l_like_save);

%number of samples
n2 = n(2);

%length of log_like_star
n_star = length(log_like_star);

l_post_weight_sort = zeros(1, n2);

l_p_w_sort_c = zeros(1, n2);

l_X_sort = zeros(1, n2);

[l_like_save_sort, index_sort] = sort(l_like_save);

C = unique(l_like_save_sort);

Ncount = histc(l_like_save, C);

c_sum_Ncount = cumsum(Ncount);

for i = 1:length(Ncount)
    
    if(Ncount(i) > 1)
        
        if(i == 1)
            
        index_sort_temp = index_sort(1:Ncount(i));
        index_sort(1:Ncount(i)) = index_sort_temp(datasample(1:Ncount(i), Ncount(i), 'Replace', false));
        
        end
        
        if(i > 1)
            
        low_index = c_sum_Ncount(i-1)+1;
        upper_index = (c_sum_Ncount(i-1)+Ncount(i));
            
        index_sort_temp = index_sort(low_index:upper_index);
        index_sort(low_index:upper_index) = index_sort_temp(datasample(1:Ncount(i), Ncount(i), 'Replace', false));
        
        end
    end
    
end

l_like_save_sort = l_like_save(index_sort);

level_assign_sort = level_assign(index_sort);

j = 1;

h = waitbar(0,'Initialize...');
for i = 2:(n_star+1)
    
    n_in_level = length(level_assign_sort(level_assign_sort == i));
    
    if(i ~= n_star+1)
        
        u_sample = zeros(1, n_in_level);
    
        for l = 1:n_in_level
            
            upper_index = i;
            lower_index = i - 1;
        
            u_sample(1, l) = X_sample_save_median(1, upper_index) + (X_sample_save_median(1, lower_index) - X_sample_save_median(1, upper_index))*rand(1);
        
        end
        
        u_sample_sort = sort(u_sample, 'descend');
        
        for l = 1:n_in_level
            
            l_X_sort(1, j) = log(u_sample_sort(1, l));
            
            if(j == 1)
                
                l_post_weight_sort(1, j) = l_like_save_sort(1,j) + LogDiffExp_func(0, l_X_sort(1, j));
                
                if(imag(l_post_weight_sort(1, j)))
                    
                    l_post_weight_sort(1, j) = -inf;
                    
                end
                
                l_p_w_sort_c(1, j) = l_post_weight_sort(1, j);
                
            else
            
                l_post_weight_sort(1, j) = l_like_save_sort(1,j) + LogDiffExp_func(l_X_sort(1, j-1), l_X_sort(1, j));
            
                
                if(imag(l_post_weight_sort(1, j)))
                    
                    l_post_weight_sort(1, j) = -inf;
                    
                end
                
                l_p_w_sort_c(1, j) = LogSumExp_func(l_post_weight_sort(1, 1:j));
                
            end
            
            j = j + 1;
            
            waitbar(j/n2,h,sprintf('Finding log posterior weights %d%%',(j/n2)*100))
            
        end
        
    end
    
    if(i == n_star+1)
        
        u_sample = zeros(1, n_in_level);
        
        for l = 1:n_in_level
            
            u_sample(1, l) = X_sample_save_median(1, n_star)*rand(1);
        
        end
        
        u_sample_sort = sort(u_sample, 'descend');
        
        for l = 1:n_in_level
            
            l_X_sort(1, j) = log(u_sample_sort(1, l));

                l_post_weight_sort(1, j) = l_like_save_sort(1,j) + LogDiffExp_func(l_X_sort(1, j-1), l_X_sort(1, j));
            
                if(imag(l_post_weight_sort(1, j)))
                    
                    l_post_weight_sort(1, j) = -inf;
                    
                end
                
                l_p_w_sort_c(1, j) = LogSumExp_func(l_post_weight_sort(1, 1:j));
                
            j = j + 1;
            
            waitbar(j/n2,h,sprintf('Finding log posterior weights %d%%',(j/n2)*100))
            
        end
        
    end

end
close(h)

l_p_w_sort_c = l_p_w_sort_c - LogSumExp_func(l_post_weight_sort);

l_post_weight_sort = l_post_weight_sort - LogSumExp_func(l_post_weight_sort);

end