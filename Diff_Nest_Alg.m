function [log_like_sample, visits_save, log_like_save, log_prior_save, theta_save, theta_level_save, X_sample_save, log_Z_final, part_max] = Diff_Nest_Alg(theta0, logfuns, T, J, C, Num_to_Create_Level, lambda, howMany, widths, Tol, save_interval, visualize, par_names, data, fitting_plot_fun, beta_power)
%Copyright (C) 2022 Weston C. Roda and Donglin Han
%Diffusive Nested Sampling algorithm (Matlab implementation)
%
%License
%GNU General Public License v3.0
%
%Derived Work
%This software is a matlab implementation of the Brewer, Partay, Csanyi,
%Foreman-Mackey Diffusive Nested Sampling algorithm based off of elements 
%of DNest5 ((c) 2020 Brendon J. Brewer) and the following papers:
%
%Brewer BJ, Partay LB, Csanyi G (2011) Diffusive Nested Sampling. Stat
%Comput 21(4):649-656. https://doi.org/10.1007/s11222-010-9198-8
%
%Brewer BJ, Foreman-Mackey D (2016) DNest4: Diffusive Nested Sampling in 
%C++ and Python. J Stat Softw 86(7):1-33. 
%http://dx.doi.org/10.18637/jss.v086.i07
%
%Skilling J (2006) Nested Sampling for General Bayesian Computation.
%Bayesian Anal. 1(4): 833-859. https://doi.org/10.1214/06-BA127
%
%Inputs:
%
%Num_to_Create_Level: the number of samples that are needed above the 
%current likelihood cut off to create another level
%
%T: the number of samples that are at least needed after the creation of a 
%level to create another level. In this code it is assumed that 
%Num_to_Create_Level < T.
%
%J: the maximum number of levels that are allowed for the DNS sampler.
%If the sampler has not reached the log likelihood tolerance between levels
%before the Jth level creation, then the DNS sampler will stop after the 
%creation of the Jth level.
%
%lambda: the amount that the particles are allowed to backtrack during
%the sampling. This is the diffusive part of the DNS sampler.
%
%howMany: the number of samples to use during the second phase of the 
%sampler
%[The first phase of the DNS sampler builds the likelihood cut offs 
%(levels). The second phase of the DNS sampler randomly samples between the
%built levels. howMany are the number of samples to use during the second 
%phase of the sampler. The second phase of the algorithm continues to 
%sample and improve the estimate of Z.]
%
%Tol: if the absolute value of the difference between log likelihoods of 
%two subsequent levels are less than the tolerance, Tol, then the DNS 
%sampler moves to the second stage of the sampling algorithm.
%
%C: this number provides the amount of confidence in the theoretical 
%expectation X_{j+1} = exp(-1) X_{j}.
%
%save_interval: this parameter thins the MCMC chain to provide random 
%samples
%
%visualize: if the boolean parameter "visualize" is set to true, then the
%particle tracking plot, fitting plot, and Z estimates plot are displayed.
%Otherwise, these plots are not displayed.
%
%par_names: this parameter holds the string names of the parameters to be
%estimated
%
%data: if the boolean parameter "visualize" is set to true, the vector or
%matrix "data" will be passed to the "fitting_plot_fun" function to display the
%data along with the currently best fitted model.
%
%beta_power: in the second phase of the Diffusive Nested sampling, this 
%parameter controls the strength of the effect to correct the mass X values
%
%Outputs:
%
%log_like_sample: the log likelihood cut off values for each level are 
%stored in this 1x(number of levels) vector
%
%visits_save: This 
%(number of levels)x(number of particles)x(saved iterations) 
%matrix stores the visits information for each level during the second 
%phase of the sampler.
%
%log_like_save: this cell array holds the log likelihood value for the 
%thinned MCMC particles. Each particle is a separate cell in the cell 
%array.
%
%log_prior_save: this cell array holds the log prior value for the thinned 
%MCMC particles. Each particle is a separate cell in the cell array.
%
%theta_save: this cell array holds the theta vector for the thinned MCMC 
%particles. Each particle is a separate cell in the cell array.
%
%theta_level_save: this cell array holds the level value for the thinned 
%MCMC particles. Each particle is a separate cell in the cell array.
%
%X_sample_save: this (saved iterations)x(number of levels) matrix stores 
%the prior mass X for each level during the second phase of the sampler.
%
%log_Z_final: this 1x(saved iterations) vector stores the estimates of the 
%log evidence Z during the second phase of the sampler.
%
%part_max: This number is the particle used during the second stage of the
%algorithm. This particle was chosen because it was most successful in
%reaching the maximum likelihood.

%Save elements
log_like_save = [];
log_prior_save = [];
theta_save = [];
theta_level_save = [];

%If the specified tolerance, Tol, has been reached exit the level creation 
%step early
tol_reached = 0;

%log_like_sample holds the cut off values of log L star from level 0 to
%level J
log_like_sample = zeros(1, J+1);

log_like_sample(1,1) = -Inf;

n = size(theta0);

%n1 is the number of parameters
n1 = n(1);

%n2 is the number of initial points (also called particles)
n2 = n(2);

for i = 1:n2
    
    log_like_save{i} = [];
    log_prior_save{i} = [];
    theta_save{i} = [];
    theta_level_save{i} = [];
    
end

theta_current_all = theta0;
theta_current_all_level = zeros(1,n(2));

%Currently the best log likelihood value found
log_like_best = -Inf;

%Currently the best log likelihood value found
par_best = theta_current_all(:,1);

%The second phase of the algorithm
second_phase = 0;

%Iteration
t = repmat(2,1,n2);

%Number of particles initially
num_particles = n2;

for i = 1:n2
        
    theta_Alg{i} = theta_current_all(:,i);
    theta_level_Alg{i}(1,1) = theta_current_all_level(1,i);
    log_like_Alg{i}(1,1) = logfuns{2}(theta_current_all(:,i)); 
    log_prior_Alg{i}(1,1) = logfuns{1}(theta_current_all(:,i));
        
end

%Col 1 of "visits" matrix: number of visits to level i
%Col 2 of "visits" matrix: number of visits to level i larger than L_{i+1}
%Col 3 of "visits" matrix: number of visits to level i ever
%Col 4 of "visits" matrix: number of proposals accepted to level i
%Col 5 of "visits" matrix: number of proposals in total to level i
%There are J+1 rows in the visits matrix from level 0 to level J. Once the
%tolerance has been reached, the top level is reduced to the last created
%level.

for i = 1:n2
        
    visits{i} = zeros(J+1, 5);
        
end

visits_sum = zeros(J+1, 5);

%stash
stash = [];
a = 1;

%tracks the cpu time for each level creation
time_cpu = [];

%Generate the levels 
%for m = 1:J

m = 1;

while(m <= J && tol_reached == 0)
    
tStart = cputime;

disp("Level "+m);

have_enough = 0;

while(have_enough == 0)
    
    continue_iteration = 1;
    
    counter_T = 1;
        
        stash_display = waitbar(0,'Initialize...');
        
        while(continue_iteration == 1)
            
            k = randi(num_particles);
                
                which_action = rand(1);
                
                if(which_action < 0.5)
                %jump first
                %then move
            
                [theta_level_Alg,X_sample,W_sample,~] = jumpDNS(which_action, log_like_sample, log_like_Alg, theta_level_Alg, t, m, k, C, visits_sum, lambda, second_phase, log_like_best, NaN, NaN, NaN, NaN, beta_power);
            
            
                [theta_Alg,log_prior_Alg,log_like_Alg,my_changed_theta,log_like_best_changed, par_best_changed] = moveDNS(which_action, log_like_sample, theta_level_Alg, theta_Alg, log_prior_Alg, log_like_Alg, logfuns, t, k, log_like_best, widths, n1, par_best);
            
                log_like_best = log_like_best_changed;
                par_best = par_best_changed;
                
                visits{k}(theta_level_Alg{k}(:,t(k))+1, 5) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 5) + 1;
                visits_sum(theta_level_Alg{k}(:,t(k))+1, 5) = visits_sum(theta_level_Alg{k}(:,t(k))+1, 5) + 1;
            
                if(my_changed_theta == 1)
                
                    visits{k}(theta_level_Alg{k}(:,t(k))+1, 4) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 4) + 1;
                    visits_sum(theta_level_Alg{k}(:,t(k))+1, 4) = visits_sum(theta_level_Alg{k}(:,t(k))+1, 4) + 1;
                
                end
            
                else
            
                [theta_Alg,log_prior_Alg,log_like_Alg,my_changed_theta,log_like_best_changed, par_best_changed] = moveDNS(which_action, log_like_sample, theta_level_Alg, theta_Alg, log_prior_Alg, log_like_Alg, logfuns, t, k, log_like_best, widths, n1, par_best);
            
                log_like_best = log_like_best_changed;
                par_best = par_best_changed;
                
                visits{k}(theta_level_Alg{k}(:,t(k)-1)+1, 5) = visits{k}(theta_level_Alg{k}(:,t(k)-1)+1, 5) + 1;
                visits_sum(theta_level_Alg{k}(:,t(k)-1)+1, 5) = visits_sum(theta_level_Alg{k}(:,t(k)-1)+1, 5) + 1;
                
                if(my_changed_theta == 1)
                
                    visits{k}(theta_level_Alg{k}(:,t(k)-1)+1, 4) = visits{k}(theta_level_Alg{k}(:,t(k)-1)+1, 4) + 1;
                    visits_sum(theta_level_Alg{k}(:,t(k)-1)+1, 4) = visits_sum(theta_level_Alg{k}(:,t(k)-1)+1, 4) + 1;
                
                end
                
                [theta_level_Alg,X_sample,W_sample,~] = jumpDNS(which_action, log_like_sample, log_like_Alg, theta_level_Alg, t, m, k, C, visits_sum, lambda, second_phase, log_like_best, NaN, NaN, NaN, NaN, beta_power);
            
                end
        
                visits{k}(theta_level_Alg{k}(:,t(k))+1, 3) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 3) + 1;
                visits_sum(theta_level_Alg{k}(:,t(k))+1, 3) = visits_sum(theta_level_Alg{k}(:,t(k))+1, 3) + 1;
        
            if(theta_level_Alg{k}(:,t(k)) < m-1)
                
                visits{k}(theta_level_Alg{k}(:,t(k))+1, 1) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 1) + 1;
                visits_sum(theta_level_Alg{k}(:,t(k))+1, 1) = visits_sum(theta_level_Alg{k}(:,t(k))+1, 1) + 1;
                
                if(log_like_Alg{k}(:,t(k)) > log_like_sample(1,theta_level_Alg{k}(:,t(k))+2))
                    
                    visits{k}(theta_level_Alg{k}(:,t(k))+1, 2) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 2) + 1;
                    visits_sum(theta_level_Alg{k}(:,t(k))+1, 2) = visits_sum(theta_level_Alg{k}(:,t(k))+1, 2) + 1;
                    
                end
            
            end
            
            if(log_like_Alg{k}(:,t(k)) > log_like_sample(1,m))
                
                stash(1, a) = log_like_Alg{k}(:,t(k));
                
                a = a + 1;
                
            end
            
            if(a <= Num_to_Create_Level)
                
                waitbar(a/Num_to_Create_Level,stash_display,sprintf('Stash Display'))
                
            end
            
            if(a > Num_to_Create_Level)
                
                waitbar(counter_T/T,stash_display,sprintf('Counter T Finishing'))
                
            end
            
            if(mod(sum(t-1), save_interval) == 0)
            
            log_like_save{k} = [log_like_save{k} log_like_Alg{k}(:,t(k))];
            log_prior_save{k} = [log_prior_save{k} log_prior_Alg{k}(:,t(k))];
            theta_save{k} = [theta_save{k} theta_Alg{k}(:,t(k))];
            theta_level_save{k} = [theta_level_save{k} theta_level_Alg{k}(:,t(k))];
            
            end
            
            t(k) = t(k) + 1;
            counter_T = counter_T + 1;
            
            if(length(stash) >= Num_to_Create_Level && counter_T > T)
                
                continue_iteration = 0;
                
            end
             
        end
    
    log_L_s_temp = quantile(stash,1-exp(-1));
    
    stash = stash(stash > log_L_s_temp);
    
    a = length(stash) + 1;
    
    if(isempty(stash))
        
        stash = [];
        a = 1;
        
    end
    
    close(stash_display)
        
        disp("Currently best Log Likelihood: "+log_like_best)
        
        tEnd = cputime - tStart;
        
        disp("Time to create current level: "+datestr(tEnd/(60*60*24), 13))
        
        time_cpu(1, m) = tEnd;
        
        disp("Total time to create up to "+m+" levels: "+datestr(sum(time_cpu)/(60*60*24), 13))
        
        disp("Estimated time left to reach the maximum level of "+J+": "+datestr((mean(time_cpu)*(J-m))/(60*60*24), 13))
        
        have_enough = 1;
    
end

log_L_s = log_L_s_temp;

log_like_sample(1,m+1) = log_L_s;

disp("Log Likelihood level cutoff: "+log_L_s)

good = zeros(1, n2);
num_bad = 0;

max_log_W = -inf;

for l = 1:n2
    
    if(theta_level_Alg{l}(1, length(theta_level_Alg{l}(1,:))) ~= 0)
    
    log_W = log(W_sample(theta_level_Alg{l}(1, length(theta_level_Alg{l}(1,:)))));
    
    if(log_W > max_log_W)
        
        max_log_W = log_W;
        
    end
    
    end
    
end

for l = 1:n2
    
    if(theta_level_Alg{l}(1, length(theta_level_Alg{l}(1,:))) ~= 0)
    
    log_W = log(W_sample(theta_level_Alg{l}(1, length(theta_level_Alg{l}(1,:)))));
    
    prune_prob = (1 - 1/(1 + exp(-1*log_W - 4)))^3;
    
    rand_num = rand;
    
    if(rand_num <= prune_prob)
        
        good(1,l) = 1;
        num_bad = num_bad + 1;
        
    end
    
    end
    
end

if(num_bad < n2)

for l = 1:n2
    
    if(good(1,l) == 1)
        
        l_rand = randi([1 n2]);
        
        while(good(1,l_rand) == 1 || rand >= exp(log(W_sample(theta_level_Alg{l_rand}(1, length(theta_level_Alg{l_rand}(1,:))))) - max_log_W))
            
            l_rand = randi([1 n2]);
            
        end
        
        visits_sum = visits_sum - visits{l};
        
        visits{l} = visits{l_rand};
        log_like_Alg{l} = log_like_Alg{l_rand};
        log_prior_Alg{l} = log_prior_Alg{l_rand};
        theta_Alg{l} = theta_Alg{l_rand};
        theta_level_Alg{l} = theta_level_Alg{l_rand};
        t(l) = t(l_rand);
        
        visits_sum = visits_sum + visits{l};
        
        disp('prunned');
        
    end
    
end

end

s = " ";

disp("Currently best parameters")
        for i = 1:n1
            
            s = s+"("+i+") "+par_names(i)+" "+par_best(i)+"; ";
    
        end

disp(s)
        
if(visualize == true)

figure(1); clf

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 12, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

for i = 1:n2

if(good(1,i) == 0)
    
    h = plot([1:1:length(theta_level_save{i})], theta_level_save{i}, 'b');
    
end

if(good(1,i) == 1)
    
    h = plot([1:1:length(theta_level_save{i})], theta_level_save{i}, 'r');
    
end

end

xlabel('Iterations'), ylabel('Levels')

hold off

figure(2); clf
fitting_plot_fun(par_best, data)

end

if(abs(log_like_sample(1,m+1) - log_like_sample(1,m)) < Tol)
    
    tol_reached = 1;
    
else
    
    m = m+1;
    
end
    
end

%Find the particle with highest last level
    part_max_level = theta_level_Alg{1}(1, length(theta_level_Alg{1}(1,:)));
    part_max = 1;
    
    if(num_particles > 1)
    
    for i = 2:num_particles
        
        if(theta_level_Alg{i}(1, length(theta_level_Alg{i}(1,:))) > part_max_level)
            
            part_max_level = theta_level_Alg{i}(1, length(theta_level_Alg{i}(1,:)));
            part_max = i;
            
        end
        
    end
    
    end
    
    visits_save = visits{part_max};

if(tol_reached == 1)

log_like_sample = log_like_sample(1,1:m+1);

X_sample = zeros(1, m+1);
X_sample(1,1) = 1;
        
        if(m+1 >= 2)

        for i = 2:m+1
            
            X_sample(1,i) = ((visits_save(i-1,2) + C*exp(-1))/(visits_save(i-1,1) + C))*X_sample(1,i-1);
            
        end
        
        end

for i = 1:n2
    
    visits{i} = visits{i}(1:m+1,:);
        
end

visits_sum = visits_sum(1:m+1,:);

else
    
    log_like_sample = log_like_sample(1,1:m);

X_sample = zeros(1, m);
X_sample(1,1) = 1;
        
        if(m >= 2)

        for i = 2:m
            
            X_sample(1,i) = ((visits_save(i-1,2) + C*exp(-1))/(visits_save(i-1,1) + C))*X_sample(1,i-1);
            
        end
        
        end

for i = 1:n2
    
    visits{i} = visits{i}(1:m,:);
        
end

visits_sum = visits_sum(1:m,:);

end

log_Z_final = zeros(1, floor((howMany+1)/save_interval));

log_Z_final(1,1) = getLogZ(log_like_sample, log_like_best, visits{part_max}, C);

[visits{part_max}] = renormalise_visits_matrix(visits{part_max}, C);

if(tol_reached == 1)

visits_save = zeros(m+1, 5, floor((howMany+1)/save_interval));

visits_save(:, :, 1) = visits{part_max};

X_sample_save = zeros(floor((howMany+1)/save_interval), m+1);

X_sample_save(1,:) = X_sample;

m = m+1;

else
    
    visits_save = zeros(m, 5, floor((howMany+1)/save_interval));

    visits_save(:, :, 1) = visits{part_max};

    X_sample_save = zeros(floor((howMany+1)/save_interval), m);

    X_sample_save(1,:) = X_sample;
    
end

second_phase = 1;

%continue sampling between levels

%tracks the cpu time for each iteration in second phase
time_cpu = [];

disp("Sampling between levels");

count = 0;

while(count <= (howMany-1))
    
tStart = cputime;

if(mod(count, 100000) == 0 && count ~= 0)
    
    disp((count+1)+" out of "+howMany+" in second stage")
        
    disp("Total time to reach "+(count+1)+": "+datestr(sum(time_cpu)/(60*60*24), 13))
        
    disp("Estimated time left to reach "+howMany+": "+datestr((mean(time_cpu)*(howMany-(count+1)))/(60*60*24), 13))
        
        disp("Currently best Log Likelihood: "+log_like_best)
        
        if(floor((count+1)/save_interval) >= 1)
        
        disp("Last log Z Saved: "+log_Z_final(1, floor((count+1)/save_interval)))
        
        end
        
        s = " ";

        disp("Currently best parameters")
        for i = 1:n1
            
            s = s+"("+i+") "+par_names(i)+" "+par_best(i)+"; ";
    
        end

        disp(s)
        
if(visualize == true)
    
    figure(1); clf
    
    hold on
    set(0,'defaultLineLineWidth',1.5);   
    set(0,'defaultLineMarkerSize',9);
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultTextFontName', 'Arial');

    for i = 1:n2

        if(good(1,i) == 0)
            
            h = plot([1:1:length(theta_level_save{i})], theta_level_save{i}, 'b');
    
        end
        
        if(good(1,i) == 1)
            
            h = plot([1:1:length(theta_level_save{i})], theta_level_save{i}, 'r');
    
        end

    end
    
    xlabel('Iterations'), ylabel('Levels')

    hold off
    
    figure(2); clf
    fitting_plot_fun(par_best, data)
    
    figure(3); clf
    
    hold on
    set(0,'defaultLineLineWidth',1.5);   
    set(0,'defaultLineMarkerSize',9);
    set(gca, 'FontSize', 12, 'LineWidth', 1);
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultTextFontName', 'Arial');

        %plot([1:1:(count+1)], log_Z_final(1, 1:(count+1)))
        plot([1:1:length(log_Z_final)], log_Z_final)
        
        xlabel('Iterations'), ylabel('log Z estimate')

    hold off
    
end

end

                k = part_max;
                
                which_action = rand(1);
                
                if(which_action < 0.5)
                %jump first
                %then move
            
                [theta_level_Alg,X_sample,~,log_Z_final] = jumpDNS(which_action, log_like_sample, log_like_Alg, theta_level_Alg, t, m, k, C, visits_sum, lambda, second_phase, log_like_best, visits{part_max}, log_Z_final, count, save_interval, beta_power);
            
                if(mod(count+2, save_interval) == 0)
                    
                    visits_save(:, :, (count+2)/save_interval) = visits{part_max};
                    X_sample_save((count+2)/save_interval,:) = X_sample;
                    
                end
            
                [theta_Alg,log_prior_Alg,log_like_Alg,my_changed_theta,log_like_best_changed, par_best_changed] = moveDNS(which_action, log_like_sample, theta_level_Alg, theta_Alg, log_prior_Alg, log_like_Alg, logfuns, t, k, log_like_best, widths, n1, par_best);
            
                log_like_best = log_like_best_changed;
                par_best = par_best_changed;
                
                visits{k}(theta_level_Alg{k}(:,t(k))+1, 5) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 5) + 1;
            
                if(my_changed_theta == 1)
                
                    visits{k}(theta_level_Alg{k}(:,t(k))+1, 4) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 4) + 1;
                
                end
            
                else
            
                [theta_Alg,log_prior_Alg,log_like_Alg,my_changed_theta,log_like_best_changed, par_best_changed] = moveDNS(which_action, log_like_sample, theta_level_Alg, theta_Alg, log_prior_Alg, log_like_Alg, logfuns, t, k, log_like_best, widths, n1, par_best);
            
                log_like_best = log_like_best_changed;
                par_best = par_best_changed;
                
                visits{k}(theta_level_Alg{k}(:,t(k)-1)+1, 5) = visits{k}(theta_level_Alg{k}(:,t(k)-1)+1, 5) + 1;
                
                if(my_changed_theta == 1)
                
                    visits{k}(theta_level_Alg{k}(:,t(k)-1)+1, 4) = visits{k}(theta_level_Alg{k}(:,t(k)-1)+1, 4) + 1;
                    
                end
                
                [theta_level_Alg, X_sample,~,log_Z_final] = jumpDNS(which_action, log_like_sample, log_like_Alg, theta_level_Alg, t, m, k, C, visits_sum, lambda, second_phase, log_like_best, visits{part_max}, log_Z_final, count, save_interval, beta_power);
            
                if(mod(count+2, save_interval) == 0)
                
                    visits_save(:, :, (count+2)/save_interval) = visits{part_max};
                    X_sample_save((count+2)/save_interval,:) = X_sample;
                
                end
            
                end
        
                visits{k}(theta_level_Alg{k}(:,t(k))+1, 3) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 3) + 1;
                
                
            if(theta_level_Alg{k}(:,t(k)) < m-1)
                
                visits{k}(theta_level_Alg{k}(:,t(k))+1, 1) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 1) + 1;
                
                if(log_like_Alg{k}(:,t(k)) > log_like_sample(1,theta_level_Alg{k}(:,t(k))+2))
                    
                    visits{k}(theta_level_Alg{k}(:,t(k))+1, 2) = visits{k}(theta_level_Alg{k}(:,t(k))+1, 2) + 1;
                    
                end
            
            end
            
            if(mod(sum(t-1), save_interval) == 0)
                
            log_like_save{k} = [log_like_save{k} log_like_Alg{k}(:,t(k))];
            log_prior_save{k} = [log_prior_save{k} log_prior_Alg{k}(:,t(k))];
            theta_save{k} = [theta_save{k} theta_Alg{k}(:,t(k))];
            theta_level_save{k} = [theta_level_save{k} theta_level_Alg{k}(:,t(k))];
            
            end
            
            tEnd = cputime - tStart;
        
            time_cpu(1, count+1) = tEnd;
            
            t(k) = t(k) + 1;
            count = count +1;
            
end

end

