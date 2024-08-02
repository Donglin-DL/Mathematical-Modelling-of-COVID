function [theta_samples, l_like_samples, l_prior_samples, l_post_samples, l_post_weight_samples, l_X_samples, index_used] = find_log_post(theta_save, l_like_save, l_prior_save, log_Z_median, l_post_weight_sort, l_X_sort, index_sort, num_iterations, l_p_w_sort_c)
%The "find_log_post" function determines the representative log posterior 
%sample, and outputs the following information for each sample in the 
%representative log posterior sample: theta value, log likelihood value,
%log prior value, log posterior value, log posterior weight, and log prior 
%mass value

index_used = [];

theta_save_sort = theta_save(:, index_sort);
l_like_save_sort = l_like_save(index_sort);
l_prior_save_sort = l_prior_save(index_sort);

theta_samples = [];
l_like_samples = [];
l_prior_samples = []; 
l_post_samples = [];
l_post_weight_samples = [];
l_X_samples = [];

h = waitbar(0,'Initialize...');
for i = 1:num_iterations
    
    %Randomly choose a sample based off of the cumulative posterior weight 
    %of the samples.
    
    w1 = rand(1);
    
    my_diff = w1 - exp(l_p_w_sort_c);
    
    [~, min_index] = min(abs(my_diff));
    
    index_used(1, i) = min_index;
    
    theta_samples(:, i) = theta_save_sort(:, min_index);
    l_like_samples(1, i) = l_like_save_sort(1, min_index);
    l_prior_samples(1, i) = l_prior_save_sort(1, min_index); 
    l_post_samples(1, i) = l_like_save_sort(1, min_index) - log_Z_median + l_prior_save_sort(1, min_index);
    l_post_weight_samples(1, i) = l_post_weight_sort(1, min_index);
    l_X_samples(1, i) = l_X_sort(1, min_index);
    
    waitbar(i/num_iterations,h,sprintf('Finding representative posterior sample %d%%',(i/num_iterations)*100))
end
close(h)

end

