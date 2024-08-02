function [] = post_silhouette_plot(theta_samples,l_post_samples, par_names, par_num)

b1 = theta_samples(par_num,:);

e = exp(l_post_samples);

[tmp, ind] = sort(b1);
e_sort = e(ind);

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot(tmp, e_sort)

xlabel(par_names(1,par_num), 'FontSize', 12), ylabel('Posterior Density', 'FontSize', 12)

hold off

end

