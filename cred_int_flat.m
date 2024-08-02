function [ out ] = cred_int_flat( k1, lelik, maxlelik)

lower = 0;
upper = 0;

v = [];
j = 1;

for i = 1:length(lelik)
    
    if(lelik(1, i) >= maxlelik)
        
        v(1, j) = k1(1, i);
        
        j = j+1;
        
    end;
    
end;

hold on
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 1;    % AxesLineWidth
fsz = 12;      % Fontsize
lw = 1.5;      % LineWidth
msz = 9;  
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot(k1, lelik)
plot([min(v) min(v)], [0 max(lelik)], 'r')
plot([max(v) max(v)], [0 max(lelik)], 'r')
plot([min(v) max(v)], [maxlelik maxlelik], 'g')
xlabel('par', 'FontSize', 12), ylabel('Unnormalized Posterior Probability', 'FontSize', 12)

hold off

lower = min(v);
upper = max(v);

out = [lower upper];

end

