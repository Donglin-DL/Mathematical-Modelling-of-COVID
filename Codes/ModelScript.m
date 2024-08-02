
%Alberta First Wave (March 9th to May 19th)

%case data
data_case = [1 0 0 2 0 12 7 16 6 17 17 22 22 19 24 38 36 37 52 36 18 22 79 61 72 47 47 42 29 42 29 40 52 53 68 65 128 130 145 215 186 196 189 272 294 349 238 229 207 180 261 236 224 139 118 69 76 63 64 82 75 85 73 45 70 63 61 71 54 46 38 44];

%testing data
data_rho = [0.001085187 0 0 0.000573394 0 0.002282888 0.001616628 0.002987304 0.000982238 0.002951645 0.002162575 0.004269358 0.005236225 0.005142084 0.006166495 0.010412385 0.011131725 0.009521359 0.014016173 0.013675214 0.006497022 0.008865605 0.027724162 0.030108588 0.030271179 0.020968102 0.030470016 0.029587883 0.015508021 0.017287508 0.010736764 0.010372099 0.023471 0.019857625 0.030720578 0.01753676 0.029547553 0.024602574 0.03616411 0.061393489 0.074984882 0.085626911 0.052209945 0.08220006 0.084519189 0.102316036 0.072983747 0.116509794 0.099209202 0.049382716 0.091354568 0.090944123 0.085774459 0.063937443 0.07375 0.044145873 0.025270158 0.017080114 0.020259576 0.028691393 0.040149893 0.054140127 0.057277364 0.017989206 0.028305702 0.024896266 0.028365496 0.036902287 0.031671554 0.031039136 0.030170703 0.084130019];

%take nonzero rho data
data_rho2 = data_rho(data_rho > 0); 

%find minimum of the rho data and use this as the minimum value for the rho
%function
min(data_rho2)
min_rho = min(data_rho2);

%fit guass3 model for death data and prepare for later. 
x=[1:1:72];
y=Aver_death_new_report;
f=fit(x.',y.','gauss3');
plot(f,x,y)
close

death_mean = [];
for i=1:72
    death_mean(i)=f(i);
end

save death_mean.mat death_mean

widths = [1e-9 1e-7; 
    0 1.5; 
    1 79; 
    1 300; 
    1/11 1/5; 
    0.01 1; 
    0.01 1;
    1e-4 0.03; 
    1 10; 
    39 55; 
    52 60; 
    57 length(data_case); 
    1e-3 0.1; 
    1e-3 0.1; 
    29 32; 
    65 68; 
    10 37; 
    0.01 0.04];

prior1 = @(beta) unifpdf(beta,widths(1,1),widths(1,2));    
prior2 = @(f) unifpdf(f,widths(2,1),widths(2,2));            
prior3 = @(p) unifpdf(p, widths(3,1),widths(3,2));      
prior4 = @(I0) unifpdf(I0,widths(4,1),widths(4,2));         
prior5 = @(gamma) unifpdf(gamma,widths(5,1),widths(5,2));         
prior6 = @(q) unifpdf(q,widths(6,1),widths(6,2));            
prior7 = @(s) unifpdf(s,widths(7,1),widths(7,2));            
prior8 = @(p2) unifpdf(p2,widths(8,1),widths(8,2));      
prior9 = @(time1) unifpdf(time1,widths(9,1),widths(9,2));      
prior10 = @(time2) unifpdf(time2,widths(10,1),widths(10,2));      
prior11 = @(time3) unifpdf(time3,widths(11,1),widths(11,2));      
prior12 = @(time4) unifpdf(time4,widths(12,1),widths(12,2));   
prior13 = @(height_val) unifpdf(height_val,widths(13,1),widths(13,2));      
prior14 = @(const_val) unifpdf(const_val,widths(14,1),widths(14,2));      
prior15 = @(time1beta) unifpdf(time1beta,widths(15,1),widths(15,2));  
prior16 = @(time2beta) unifpdf(time2beta,widths(16,1),widths(16,2)); 
prior17 = @(t1b) unifpdf(t1b,widths(17,1),widths(17,2)); 
prior18 = @(m0_val) unifpdf(m0_val,widths(18,1),widths(18,2)); 
     

%Population
Pop = 4371323.74;

data = [data_case; data_rho];

par_names = ["beta","f","p","I0","gamma","q","s","p2","time1","time2","time3","time4","height_val","const_val","time1beta","time2beta","t1b","m0_val"];

%'loglike' is the natural logorithm of the likelihood function
loglike = @(par) LogLikeFunc(par,length(data_case), Pop, data_case, data_rho, min_rho);

logprior = @(par) log(prior1(par(1))) + log(prior2(par(2))) + log(prior3(par(3))) + log(prior4(par(4))) + log(prior5(par(5))) + log(prior6(par(6))) + log(prior7(par(7)))  + log(prior8(par(8))) + log(prior9(par(9))) + log(prior10(par(10))) + log(prior11(par(11))) + log(prior12(par(12))) + log(prior13(par(13))) + log(prior14(par(14))) + log(prior15(par(15))) + log(prior16(par(16))) + log(prior17(par(17))) + log(prior18(par(18))) + penalty(par);

%Generate the vector of initial conditions
theta0 = [6.1278e-08 0.2189 30 277.3939 0.1340 0.3551 0.3888 0.01 9.2119 40 54.0604 59 0.0725 0.0314 30 67 35 0.03; 
    6.1278e-08 0.2189 30 277.3939 0.1340 0.3551 0.3888 0.01 9.2119 40 54.0604 59 0.0725 0.0314 30 67 35 0.03; 
    6.1278e-08 0.2189 30 277.3939 0.1340 0.3551 0.39 0.01 9.2119 40 54.0604 59 0.0725 0.0314 30 67 35 0.03; 
    6.1278e-08 0.2189 30 277.3939 0.1340 0.3551 0.3888 0.01 9.2119 40 54.0604 59 0.0725 0.0314 30 67 35 0.03; 
    6.1278e-08 0.2189 30 277.3939 0.1340 0.3551 0.39 0.01 9.2119 40 54.0604 59 0.0725 0.0314 30 67 35 0.03]';

save theta0.mat theta0

logfuns = {@(m)logprior(m) @(m)loglike(m)};

logprior_test = {@(m)log(prior1(m)) @(m)log(prior2(m)) @(m)log(prior3(m)) @(m)log(prior4(m)) @(m)log(prior5(m)) @(m)log(prior6(m)) @(m)log(prior7(m)) @(m)log(prior8(m)) @(m)log(prior9(m)) @(m)log(prior10(m)) @(m)log(prior11(m)) @(m)log(prior12(m)) @(m)log(prior13(m)) @(m)log(prior14(m)) @(m)log(prior15(m)) @(m)log(prior16(m)) @(m)log(prior17(m)) @(m)log(prior18(m))};

for i = 1:18
    
    disp(i)
    logprior_test{i}(theta0(i,1))
    
end


%Diffusive Nested Sampling (DNS)
Num_to_Create_Level = 100;
T = 101;
J = 500;
lambda = 90;
howMany = 1e6;
Tol = 0.0015;
C = 100;
save_interval = 10;
visualize = true;

%Fitting plot reference
fitting_plot_fun = @(par_best, data) fitting_plot(par_best, data,72,Pop,min_rho);

%In the second phase of the Diffusive Nested sampling, the parameter 
%beta_power controls the strength of the effect to correct the mass X
%values

beta_power = 100;
[log_like_sample, visits_save, log_like_save, log_prior_save, theta_save, theta_level_save, X_sample_save, log_Z_final, part_max] = Diff_Nest_Alg(theta0, logfuns, T, J, C, Num_to_Create_Level, lambda, howMany, widths, Tol, save_interval, visualize, par_names, data, fitting_plot_fun, beta_power);

save log_like_sample.mat log_like_sample
save visits_save.mat visits_save
save log_like_save.mat log_like_save
save log_prior_save.mat log_prior_save
save theta_save.mat theta_save
save theta_level_save.mat theta_level_save
save X_sample_save.mat X_sample_save
save log_Z_final.mat log_Z_final
save part_max.mat part_max


%%

%Find the median of log Z, visits matrix, and the X values
[log_Z_median, visits_median, X_sample_save_median] = median_function(log_Z_final,visits_save, X_sample_save, 3.8e3);

save log_Z_median.mat log_Z_median
save X_sample_save_median.mat X_sample_save_median

[level_assign] = find_level_assign(log_like_save{part_max}, log_like_sample);

save level_assign.mat level_assign

%The function 'find_log_post_weights' estimates the log prior mass and the 
%log posterior weight for each sample from the DNS output (output is sorted 
%by increasing likelihood value), and this information is used to obtain 
%the representative log posterior sample
[l_post_weight_sort, l_X_sort, index_sort, l_p_w_sort_c] = find_log_post_weights(log_like_save{part_max}, log_like_sample, X_sample_save_median, level_assign);

save l_post_weight_sort.mat l_post_weight_sort
save l_X_sort.mat l_X_sort
save index_sort.mat index_sort
save l_p_w_sort_c.mat l_p_w_sort_c

%View CDF
plot([1:1:length(l_p_w_sort_c)], exp(l_p_w_sort_c), 'o')

%The 'find_log_post' function determines the representative log posterior 
%sample, and outputs the following information for each sample in the 
%representative log posterior sample: theta value, log likelihood value,
%log prior value, log posterior value, log posterior weight, and log prior 
%mass value
[theta_samples, l_like_samples, l_prior_samples, l_post_samples, l_post_weight_samples, l_X_samples, index_used] = find_log_post(theta_save{part_max}, log_like_save{part_max}, log_prior_save{part_max}, log_Z_median, l_post_weight_sort, l_X_sort, index_sort, length(index_sort), l_p_w_sort_c);

save theta_samples.mat theta_samples
save l_like_samples.mat l_like_samples
save l_prior_samples.mat l_prior_samples
save l_post_samples.mat l_post_samples
save l_post_weight_samples.mat l_post_weight_samples
save l_X_samples.mat l_X_samples
save index_used.mat index_used

%%

%ki holds the parameter sample values for parameter i
k1 =  theta_samples(1,:);
k2 =  theta_samples(2,:);
k3 =  theta_samples(3,:);
k4 =  theta_samples(4,:);
k5 =  theta_samples(5,:);
k6 =  theta_samples(6,:);
k7 =  theta_samples(7,:);
k8 =  theta_samples(8,:);
k9 =  theta_samples(9,:);
k10 =  theta_samples(10,:);
k11 =  theta_samples(11,:);
k12 =  theta_samples(12,:);
k13 =  theta_samples(13,:);
k14 =  theta_samples(14,:);
k15 =  theta_samples(15,:);
k16 =  theta_samples(16,:);
k17 =  theta_samples(17,:);
k18 =  theta_samples(18,:);
%%
%This plot displays the posterior weights over the log(X) values. There
%should be a clear peak of the posterior weights and the samples to the 
%left of the peak should have small posterior weights in comparison to the
%peak weight values.
hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot(l_X_samples, exp(l_post_weight_samples), 'o')

ylabel('Posterior Weights')
xlabel('log X')

hold off

savefig("disease_model_log_X_post_weights.fig");
close
%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k1)

savefig("freq_beta.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k1, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k1_median = median(k1);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k2)

savefig("freq_f.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k2, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k2_median = median(k2);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k3)

savefig("freq_p.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k3, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k3_median = median(k3);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k4)

savefig("freq_I0.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k4, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k4_median = median(k4);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k5)

savefig("freq_gamma.fig");
close 

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k5, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k5_median = median(k5);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k6)

savefig("freq_q.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k6, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k6_median = median(k6);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k7)

savefig("freq_s.fig");
close 

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k7, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k7_median = median(k7);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k8)

savefig("freq_p2.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k8, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k8_median = median(k8);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k9)

savefig("freq_time1.fig");
close 

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k9, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k9_median = median(k9);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k10)

savefig("freq_time2.fig");
close 

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k10, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k10_median = median(k10);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k11)

savefig("freq_time3.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k11, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k11_median = median(k11);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k12)

savefig("freq_time4.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k12, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k12_median = median(k12);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k13)

savefig("freq_height_val.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k13, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k13_median = median(k13);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k14)

savefig("freq_const_val.fig");
close 

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k14, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k14_median = median(k14);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k15)

savefig("freq_time1beta.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k15, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k15_median = median(k15);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k16)

savefig("freq_time2beta.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k16, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k16_median = median(k16);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k17)

savefig("freq_time1b.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k17, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k17_median = median(k17);

%%
%View the frequency of the parameter sample values during the second stage
%of sampling (proportionally this frequency provides the estimated marginal 
%posterior distribution for this parameter)
histogram(k18)

savefig("freq_start_val.fig");
close

%Find the Highest Posterior Density (HPD) with alpha = 0.05 (This gives the 
%lower and upper values of the parameter that contains 95% of the density. 
%It is a 95% credible interval.)
[l_HPD, u_HPD] = HPD(k18, 0.05);

%Find the median of the parameter sample values during the second stage of 
%sampling
k18_median = median(k18);
%%
%View the silhouette of the posterior surface from the x1 perspective
post_silhouette_plot(theta_samples,l_post_samples, par_names, 1)
savefig("post_surf_beta.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 2)
savefig("post_surf_f.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 3)
savefig("post_surf_p.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 4)
savefig("post_surf_I0.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 5)
savefig("post_surf_gamma.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 6)
savefig("post_surf_q.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 7)
savefig("post_surf_s.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 8)
savefig("post_surf_p2.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 9)
savefig("post_surf_time1.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 10)
savefig("post_surf_time2.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 11)
savefig("post_surf_time3.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 12)
savefig("post_surf_time4.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 13)
savefig("post_surf_height_val.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 14)
savefig("post_surf_const_val.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 15)
savefig("post_surf_time1beta.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 16)
savefig("post_surf_time2beta.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 17)
savefig("post_surf_time1b.fig");
close

%%
post_silhouette_plot(theta_samples,l_post_samples, par_names, 18)
savefig("post_surf_start_val.fig");
close

%%
%parameter values at the maximum posterior value
[M, I] = max(exp(l_post_samples));

my_model_sol = @(par) model_all_out(par,length(data_case),Pop,min_rho);

model_sol = my_model_sol(theta_samples(:,I));

my_data = data;

n1=length(data_case);

figure(2); clf

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot([1:1:n1], model_sol(1,:), 'r')

xlabel('Days'), ylabel('S')
savefig('best_fit_S')

hold off
close

figure(3); clf

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot([1:1:n1], model_sol(2,:), 'r')

xlabel('Days'), ylabel('I')
savefig('best_fit_I')

hold off
close

figure(4); clf

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot([1:1:n1], model_sol(3,:), 'r')
plot([1:1:n1], data_case, 'ro')

xlabel('Days'), ylabel('cases')
savefig('best_fit_cases')

hold off
close
%%

%The following 'for loop' estimates the posterior predictive distribution
%for the model solution x

post_samples=l_post_samples;

my_model_sol = @(par) model_all_out(par,n1,Pop,min_rho);

data_predict = zeros(length(data_case),length(post_samples));

S_predict = zeros(length(data_case),length(post_samples));

I_predict = zeros(length(data_case),length(post_samples));

betat_predict = zeros(length(data_case),length(post_samples));

case_out = zeros(length(data_case),length(post_samples));

rho_out = zeros(length(data_case),length(post_samples));

D_predict = zeros(length(data_case),length(post_samples));

rng(99)

h = waitbar(0,'Initialize...');
for i = 1:length(post_samples)
    
    theta = theta_samples(:,i);
    
    model_sol = my_model_sol(theta);
    
    out_case = model_sol(3,:);
    out_rho = model_sol(4,:);
    
    case_out(:,i) = out_case.';
    rho_out(:,i) = out_rho.';
    
    sd_case = sqrt(theta(3)*out_case);
    sd_rho = sqrt(theta(8)*out_rho);
    D_predict(:,i) = poissrnd(death_mean);
    
    data_predict(:,i) = multnormrnd_sigma_t(out_case, sd_case, 1, length(data_case));
    rho_predict(:,i) = multnormrnd_sigma_t(out_rho, sd_rho, 1, length(data_rho));
    f_rho_predict(:,i) = theta(2).*rho_predict(:,i);

    S_predict(:, i) = model_sol(1,:)';
    I_predict(:, i) = model_sol(2,:)';
    
    waitbar(i/length(post_samples),h,sprintf('%d%%',(i/length(post_samples))*100))
    
end;
close(h)


save data_predict.mat data_predict
save S_predict.mat S_predict
save I_predict.mat I_predict
save rho_predict.mat rho_predict
save f_rho_predict.mat f_rho_predict
save D_predict.mat D_predict
save case_out.mat case_out
save rho_out.mat rho_out
%%

%Predictive mean for infected population
mean_data_predict = (1/(length(post_samples)))*sum(data_predict');

%Generate 95% prediction interval
lower_predict = zeros(1, length(data_case));

upper_predict = zeros(1, length(data_case));

h = waitbar(0,'Initialize...');
for i = 1:length(data_case)
    
    lower_predict(1, i) = prctile(data_predict(i,:),2.5);
    upper_predict(1, i) = prctile(data_predict(i,:),97.5);
    
    
    waitbar(i/length(data_case),h,sprintf('%d%%',(i/length(data_case))*100))
end;
close(h)

lower_predict(lower_predict < 0) = 0;
%Plot the model fitting and prediction
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


plot([1:1:n1], data_case, 'ro')
plot([1:1:n1], mean_data_predict, 'k')
plot([1:1:n1], lower_predict, '--k')
plot([1:1:n1], upper_predict, '--k')

legend('Data','Mean', '95% Prediction Interval')

xlabel('Days'), ylabel('Case')

hold off

savefig("post_predictive_solution_Case.fig");
close
%%
%Predictive mean for infected population
mean_rho_predict = (1/(length(post_samples)))*sum(rho_predict');

%Generate 95% prediction interval
lower_predict = zeros(1, length(data_rho));

upper_predict = zeros(1, length(data_rho));

h = waitbar(0,'Initialize...');
for i = 1:length(data_rho)
    
    lower_predict(1, i) = prctile(rho_predict(i,:),2.5);
    upper_predict(1, i) = prctile(rho_predict(i,:),97.5);
    
    
    waitbar(i/length(data_rho),h,sprintf('%d%%',(i/length(data_rho))*100))
end;
close(h)



mean_rho_predict(mean_rho_predict < 0) = 0;
lower_predict(lower_predict < 0) = 0;

%Plot the model fitting and prediction
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


plot([1:1:n1], data_rho, 'ro')
plot([1:1:n1], mean_rho_predict, 'k')
plot([1:1:n1], lower_predict, '--k')
plot([1:1:n1], upper_predict, '--k')

legend('Data','Mean', '95% Prediction Interval')

xlabel('Days'), ylabel('rho')

hold off

savefig("post_predictive_solution_rho.fig");
close
%%
mean_f_rho_predict = (1/(length(post_samples)))*sum(f_rho_predict');

%Generate 95% prediction interval
lower_predict = zeros(1, length(data_rho));

upper_predict = zeros(1, length(data_rho));

h = waitbar(0,'Initialize...');
for i = 1:length(data_rho)
    
    lower_predict(1, i) = prctile(f_rho_predict(i,:),2.5);
    upper_predict(1, i) = prctile(f_rho_predict(i,:),97.5);
    
    
    waitbar(i/length(data_rho),h,sprintf('%d%%',(i/length(data_rho))*100))
end;
close(h)

mean_f_rho_predict(mean_f_rho_predict < 0) = 0;
lower_predict(lower_predict < 0) = 0;

%Plot the model fitting and prediction
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

plot([1:1:n1], mean_f_rho_predict, 'k')
plot([1:1:n1], lower_predict, '--k')
plot([1:1:n1], upper_predict, '--k')

legend('Data','Mean', '95% Prediction Interval')

xlabel('Days'), ylabel('f_rho')

hold off

savefig("post_predictive_solution_f_rho.fig");
close

%%

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

plot([1:1:n1], mean(S_predict'), 'b', 'DisplayName','Est. Mean S')
plot([1:1:n1], prctile(S_predict', 2.5), 'b--')
plot([1:1:n1], prctile(S_predict', 97.5), 'b--')

title('S')
ylabel('Number of S')
xlabel('days')

legend('Est. Mean S','Prediction interval for S')

hold off

savefig("post_predictive_solution_S.fig");
close

%%
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

plot([1:1:n1], mean(I_predict'), 'b', 'DisplayName','Est. Mean I')
plot([1:1:n1], prctile(I_predict', 2.5), 'b--')
plot([1:1:n1], prctile(I_predict', 97.5), 'b--')

title('I')
ylabel('Number of I')
xlabel('days')

legend('Est. Mean I','Prediction interval for I')

hold off

savefig("post_predictive_solution_I.fig");
close

%%
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

for i = 1:length(data_case)
    beta(i) = beta_t(i,theta_samples(1,I),theta_samples(6,I),theta_samples(7,I), theta_samples(15,I), theta_samples(16,I));
end;

plot([1:1:length(data_case)], beta, 'b', 'DisplayName','Beta(t))')
x=[30 66];
y=[beta(30) beta(66)];
stem(x, y, 'Marker', 'none')
title('Graph for Beta(t)')
str1='$$beta$$';
str2='$$q\cdot beta$$';
str3='$$s\cdot beta$$';
text(2,6.3e-8,str1,'Interpreter','latex','FontSize',10)
text(45,4.2e-8,str2,'Interpreter','latex','FontSize',10)
text(67,4.35e-8,str3,'Interpreter','latex','FontSize',10)
xlabel('days')
ylabel('Beta(t)')
legend('show')
hold off

savefig("Beta_t.fig");
close
%%
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

%plot([1:1:length(data_rho)], data_rho, 'ok', 'DisplayName','Cases/(Average)')
plot([1:1:length(data_rho)], model_sol(4,:), 'b', 'DisplayName','Est.rho(t)')

x=[3 36 46 52 58];
yy=model_sol(4,:);
y=[yy(3) yy(36) yy(46) yy(52) yy(58)];
stem(x, y, 'Marker', 'none')
title('Graph for Rho(t)')
str1='$$min\_rho$$';
str2='$$start\_val$$';
str3='$$height\_val$$';
str4='$$const\_val$$';
text(1.5,0.005,str1,'Interpreter','latex','FontSize',10)
text(30,0.025,str2,'Interpreter','latex','FontSize',10)
text(46,0.092,str3,'Interpreter','latex','FontSize',10)
text(60,0.040,str4,'Interpreter','latex','FontSize',10)
xlabel('days')
ylabel('Rho(t)')
ylim([0 0.1])

legend('show')

hold off

savefig("most_like_soln_rho.fig");
close