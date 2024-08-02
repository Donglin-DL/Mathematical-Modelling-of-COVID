function [log_Z_median, visits_median, X_sample_save_median] = median_function(log_Z_final,visits_save, X_sample_save, burn_in)
%The function "median_function" returns the median of Z, visits matrix, and
%the X values.

log_Z_final_ordered = log_Z_final(1, burn_in:length(log_Z_final));
[log_Z_final_ordered, ind] = sort(log_Z_final_ordered);
log_Z_length = length(log_Z_final_ordered);

visits_save_ordered = visits_save(:,:,burn_in:length(log_Z_final));
visits_save_ordered = visits_save_ordered(:,:,ind);

X_sample_save_ordered = X_sample_save(burn_in:length(log_Z_final),:);
X_sample_save_ordered = X_sample_save_ordered(ind,:);

is_even = 0;

if(mod(log_Z_length, 2) == 0)
    
    is_even = 1;
    
end

if(is_even == 0)
    
    median_index = ceil(log_Z_length/2);
    
    log_Z_median = log_Z_final_ordered(1,median_index);

    visits_median = visits_save_ordered(:,:,median_index);

    X_sample_save_median = X_sample_save_ordered(median_index,:);
    
else
    
    median_index_1 = log_Z_length/2;
    
    median_index_2 = (log_Z_length/2)+1;
    
    log_Z_median_1 = log_Z_final_ordered(1,median_index_1);

    visits_median_1 = visits_save_ordered(:,:,median_index_1);

    X_sample_save_median_1 = X_sample_save_ordered(median_index_1,:);
    
    log_Z_median_2 = log_Z_final_ordered(1,median_index_2);

    visits_median_2 = visits_save_ordered(:,:,median_index_2);

    X_sample_save_median_2 = X_sample_save_ordered(median_index_2,:);
    
    log_Z_median = (log_Z_median_1+log_Z_median_2)/2;

    visits_median = (1/2)*(visits_median_1+visits_median_2);

    X_sample_save_median = (1/2)*(X_sample_save_median_1+X_sample_save_median_2);
    
end

end

