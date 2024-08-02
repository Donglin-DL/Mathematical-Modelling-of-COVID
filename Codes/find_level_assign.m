function [level_assign] = find_level_assign(l_like_save, log_like_star)
%The function "find_level_assign" assigns a level to each sample based off
%of its likelihood value

n = size(l_like_save);

%number of samples
n2 = n(2);

%length of log_like_star
n_star = length(log_like_star);

level_assign = zeros(1, n2);

h = waitbar(0,'Initialize...');
for i = 1:n2
    
    value = l_like_save(1, i);
    
    [~, close_index] = min(abs(log_like_star - value));
    
    if((log_like_star(1, close_index) - value) == 0)
        
        u1 = rand;
        
        if(u1 < 0.5)
            
            level_assign(1, i) = close_index;
            
        else
            
            level_assign(1, i) = close_index + 1;
            
        end
        
    end
    
    if((log_like_star(1, close_index) - value) > 0)
        
        level_assign(1, i) = close_index;
        
    end
        
    if(close_index == n_star && (log_like_star(1, close_index) - value) < 0)

        level_assign(1, i) = close_index + 1;
        
    end
        
    if(close_index ~= n_star && (log_like_star(1, close_index) - value) < 0)
        
        level_assign(1, i) = close_index + 1;
        
    end
    
    waitbar(i/n2,h,sprintf('Assigning levels to samples %d%%',(i/n2)*100))
end
close(h)

end