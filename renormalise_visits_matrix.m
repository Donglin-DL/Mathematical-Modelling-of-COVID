function [visits] = renormalise_visits_matrix(visits, C)
%The "renormalise_visits_matrix" resets the visit matrix in the following 
%two ways:
%1) to have the same number of visits to all levels but keep the knowledge 
%about the number of visits to level i larger than L_{i+1}
%2) to have the same number of number of proposals in total to all levels
%but keep the knowledge about the number of proposals accepted to level i
%
%This adjustment to the visits matrix allows the second phase of the 
%algorithm to find better estimates of the evidence Z faster.

C_int = round(C);

n = size(visits);

n1 = n(1);

for i = 1:n1
    
    if(visits(i, 5) >= C_int)
        
        visits(i, 4) = ((visits(i, 4) + 1)/(visits(i, 5) + 1))*C_int;
        visits(i, 5) = C_int;
        
    end
    
    if(visits(i, 1) >= C_int)
        
        visits(i, 2) = ((visits(i, 2) + 1)/(visits(i, 1) + 1))*C_int;
        visits(i, 1) = C_int;
        
    end
    
end

end

