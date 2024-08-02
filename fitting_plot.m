function [] = fitting_plot(par_best, data,n,Pop,min_rho)
%The function "fitting_plot" displays the current model fitting as the DNS 
%sampler runs

%x_func was written separately from the DNS sampler
%(The user can write their own model function here.)
my_model_sol = @(par) model_all_out(par,n,Pop,min_rho);

model_sol = my_model_sol(par_best);



%The data used for this model fitting
%(The data is passed in as an input into the "Diff_Nest_Alg" function.)
my_data = data;

figure(2); clf

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot([1:1:n], model_sol(3,:), 'r')
plot([1:1:n], my_data(1,:), 'ro')

xlabel('Days'), ylabel('Cases')

hold off

figure(4); clf

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot([1:1:n], model_sol(1,:), 'r')

xlabel('Days'), ylabel('S')

hold off

figure(5); clf

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot([1:1:n], model_sol(2,:), 'r')

xlabel('Days'), ylabel('I')

hold off

figure(6); clf

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot([1:1:n], model_sol(4,:), 'r')
plot([1:1:n], my_data(2,:), 'ro')


xlabel('Days'), ylabel('rho')


hold off

end

