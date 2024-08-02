%%After run the ModelScript, run this script for more results from model. 

%The following 'for loop' estimates the posterior predictive distribution
%for the model solution x

post_samples=l_post_samples;

my_model_sol = @(par) model_SI_out(par,n1,Pop,min_rho);

S_out_predict = [];

I_out_predict = [];

betat_out_predict = [];

rng(99)

h = waitbar(0,'Initialize...');
for i = 1:length(post_samples)
    
    theta = theta_samples(:,i);
    
    model_sol = my_model_sol(theta);
    
    S_out_predict(:, i) = model_sol(:,1)';
    
    I_out_predict(:, i) = model_sol(:,2)';

    
    waitbar(i/length(post_samples),h,sprintf('%d%%',(i/length(post_samples))*100))
    
end;
close(h)

%save('S_out_predict.mat', 'S_out_predict', '-v7.3')
%save('I_out_predict.mat', 'I_out_predict', '-v7.3')

h = waitbar(0,'Initialize...');
for j = 1:length(post_samples)
    
    
    k=1;
    for i = 0:0.01:length(data_case)
        
        betat_out_predict(k,j) = beta_t(i,theta_samples(1,j),theta_samples(6,j),theta_samples(7,j), theta_samples(15,j), theta_samples(16,j));
        k = k+1;
    end;

waitbar(j/length(post_samples),h,sprintf('%d%%',(j/length(post_samples))*100))
end;
close(h)

%save('betat_out_predict.mat', 'betat_out_predict', '-v7.3')

%%
%Predictive mean for new infected population betaIS

BIS_out_predict = betat_out_predict.*S_out_predict.*I_out_predict;

%save('BIS_out_predict.mat', 'BIS_out_predict', '-v7.3')

%Integral_BIS

Infected_Proportion_Integral = sum(BIS_out_predict)*0.01;

%save Infected_Proportion_Integral.mat Infected_Proportion_Integral

%St-S0 3.0573e+04

mean(Infected_Proportion_Integral) % 3.0657e+04
prctile(Infected_Proportion_Integral,2.5) % 2.1632e+04
prctile(Infected_Proportion_Integral,97.5) % 4.6087e+04

h = waitbar(0,'Initialize...');
trapz_BIS =[];
for i = 1:length(post_samples)
    
    
    for j = 1:72
%     trapz_BIS_temp = trapz([(j-1+0.01):0.01:j],BIS_out_predict((j-1)*100+1:100*j,i));
    
    
    trapz_BIS(j,i) = trapz([0:0.01:j],BIS_out_predict(1:j*100+1,i));
    end
    
waitbar(i/length(post_samples),h,sprintf('%d%%',(i/length(post_samples))*100))
end
close(h)

save trapz_BIS.mat trapz_BIS

%St-S0 3.0573e+04

mean(trapz_BIS(72,:)) % 3.0655e+04
prctile(trapz_BIS(72,:),2.5) % 2.1631e+04
prctile(trapz_BIS(72,:),97.5) % 4.6084

%% Shadow Plot
x = 0:71; 
figure();

hold on
XDates = [datetime(2020,3,9:31) datetime(2020,4,1:30) datetime(2020,5,1:19)];
plot(XDates,mean(I_predict'),'k');
plot(XDates,prctile(I_predict', 2.5),"--k");
plot(XDates,prctile(I_predict', 97.5),"--k");

patch([x, fliplr(x)], [prctile(I_predict', 2.5),fliplr(prctile(I_predict', 97.5))],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Infections')
ylabel('Number of Infections')
xlabel('days')
xlim([datetime('9-Mar-2020') datetime('20-May-2020')]);

legend('Est. Mean Infections','95% Prediction interval for Infections')

savefig("post_predictive_solution_Infections_.fig");
close

%%
x = 0:71; 
figure();

hold on
XDates = [datetime(2020,3,9:31) datetime(2020,4,1:30) datetime(2020,5,1:19)];
plot(XDates,mean(S_predict'),'k');
plot(XDates,prctile(S_predict', 2.5),"--k");
plot(XDates,prctile(S_predict', 97.5),"--k");

patch([x, fliplr(x)], [prctile(S_predict', 2.5),fliplr(prctile(S_predict', 97.5))],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Susceptibles')
ylabel('Number of Susceptibles')
xlabel('days')
xlim([datetime('9-Mar-2020') datetime('20-May-2020')]);

legend('Est. Mean Susceptibles','95% Prediction interval for Susceptibles')

savefig("post_predictive_solution_Susceptibles_.fig");
close

%%
x = 0:71;
figure();
hold on
XDates = [datetime(2020,3,9:31) datetime(2020,4,1:30) datetime(2020,5,1:19)];
plot(XDates,mean(trapz_BIS'),'k');
plot(XDates,prctile(trapz_BIS', 2.5),"--k");
plot(XDates,prctile(trapz_BIS', 97.5),"--k");

patch([x, fliplr(x)], [prctile(trapz_BIS', 2.5),fliplr(prctile(trapz_BIS', 97.5))],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Proportion of Infected')
ylabel('Proportion of Infected')
xlabel('days')
xlim([datetime('9-Mar-2020') datetime('20-May-2020')]);

legend('Est. Mean Proportion of infected','95% Prediction interval for Proportion of Infected')

savefig("post_predictive_solution_Proportion of Infected_.fig");

%%
%Predictive mean for CFR
mean_CFR_predict = mean(CFR_predict');

mean_CFR_predict = mean_CFR_predict(15:72);

%Generate 95% prediction interval
lower_predict = zeros(1, length(data_case));

upper_predict = zeros(1, length(data_case));

h = waitbar(0,'Initialize...');
for i = 1:length(data_case)
    
    lower_predict(1, i) = prctile(CFR_predict(i,:),2.5);
    upper_predict(1, i) = prctile(CFR_predict(i,:),97.5);
    
    
    waitbar(i/length(data_case),h,sprintf('%d%%',(i/length(data_case))*100))
end;
close(h)

lower_predict(lower_predict < 0) = 0;

lower_predict = lower_predict(15:72);
upper_predict = upper_predict(15:72);

Actual_cfr = Aver_death_all_report./cumsum(data_case);
Actual_cfr = Actual_cfr(15:72);

x = 0:57;
figure();
hold on
XDates = [datetime(2020,3,23:31) datetime(2020,4,1:30) datetime(2020,5,1:19)];
plot(XDates,mean_CFR_predict,'k');
plot(XDates,lower_predict,"--k");
plot(XDates,Actual_cfr,'ob','MarkerSize',6);
plot(XDates,upper_predict,"--k");

patch([x, fliplr(x)], [lower_predict,fliplr(upper_predict)],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Case Fatality Ratio')
ylabel('Case Fatality Ratio')
xlabel('days')
xlim([datetime('23-Mar-2020') datetime('20-May-2020')]);
ylim([0 0.8])
a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
% Create a vector of '%' signs
pct = char(ones(size(a,1),1)*'%'); 
% Append the '%' signs after the percentage values
new_yticks = [char(a),pct];
% 'Reflect the changes on the plot
set(gca,'yticklabel',new_yticks)

legend('Est. Mean Case Fatality Ratio','95% Prediction interval for Case Fatality Ratio','Actual Case infection ratio')

savefig("post_predictive_solution_Case Fatality Ratio_.fig");
close

%%
%Predictive mean for IFR
mean_IFR_predict = mean(IFR_predict');

mean_IFR_predict = mean_IFR_predict(15:72);

%Generate 95% prediction interval
lower_predict = zeros(1, length(data_case));

upper_predict = zeros(1, length(data_case));

h = waitbar(0,'Initialize...');
for i = 1:length(data_case)
    
    lower_predict(1, i) = prctile(IFR_predict(i,:),2.5);
    upper_predict(1, i) = prctile(IFR_predict(i,:),97.5);
    
    
    waitbar(i/length(data_case),h,sprintf('%d%%',(i/length(data_case))*100))
end;
close(h)

lower_predict(lower_predict < 0) = 0;

lower_predict = lower_predict(15:72);
upper_predict = upper_predict(15:72);

x = 0:57;
figure();
hold on
XDates = [datetime(2020,3,23:31) datetime(2020,4,1:30) datetime(2020,5,1:19)];
plot(XDates,mean_IFR_predict,'k');
plot(XDates,lower_predict,"--k");
plot(XDates,upper_predict,"--k");

patch([x, fliplr(x)], [lower_predict,fliplr(upper_predict)],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Infection Fatality Ratio')
ylabel('Infection Fatality Ratio')
xlabel('days')
xlim([datetime('23-Mar-2020') datetime('20-May-2020')]);
ylim([0 0.02]);
%ylim([0 0.1]);
a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
% Create a vector of '%' signs
pct = char(ones(size(a,1),1)*'%'); 
% Append the '%' signs after the percentage values
new_yticks = [char(a),pct];
% 'Reflect the changes on the plot
set(gca,'yticklabel',new_yticks)

legend('Est. Mean Infection Fatality Ratio','95% Prediction interval for Infection Fatality Ratio')

savefig("post_predictive_solution_Infection Fatality Ratio_.fig");
close

%%
%Predictive mean for cases
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

x = 0:71;
figure();
hold on
XDates = [datetime(2020,3,9:31) datetime(2020,4,1:30) datetime(2020,5,1:19)];

plot(XDates, data_case, 'ob','MarkerSize',6)
plot(XDates, mean_data_predict, 'k')
plot(XDates, lower_predict, '--k')
plot(XDates, upper_predict, '--k')

patch([x, fliplr(x)], [lower_predict,fliplr(upper_predict)],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Cases')
ylabel('Number of Cases')
xlabel('days')
xlim([datetime('9-Mar-2020') datetime('20-May-2020')]);

legend('Data','Est. Mean Cases','95% Prediction interval for Cases')

savefig("post_predictive_solution_Cases_.fig");
close

%%
%Predictive mean for rhot
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

x = 0:71;
figure();
hold on
XDates = [datetime(2020,3,9:31) datetime(2020,4,1:30) datetime(2020,5,1:19)];
plot(XDates, data_rho, 'ob','MarkerSize',6)
plot(XDates, mean_rho_predict, 'k')
plot(XDates, lower_predict, '--k')
plot(XDates, upper_predict, '--k')

patch([x, fliplr(x)], [lower_predict,fliplr(upper_predict)],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Rho(t)')
ylabel('Rho(t)')
xlabel('days')
xlim([datetime('9-Mar-2020') datetime('20-May-2020')]);

legend('Data rho','Est. Mean Rho(t)','95% Prediction interval for Rho(t)')

savefig("post_predictive_solution_Rho_.fig");
close

%%
%Predictive mean for betat
mean_betat_predict = (1/(length(post_samples)))*sum(betat_predict');

%Generate 95% prediction interval
lower_predict = zeros(1, length(data_rho));

upper_predict = zeros(1, length(data_rho));

h = waitbar(0,'Initialize...');
for i = 1:length(data_rho)
    
    lower_predict(1, i) = prctile(betat_predict(i,:),2.5);
    upper_predict(1, i) = prctile(betat_predict(i,:),97.5);
    
    
    waitbar(i/length(data_rho),h,sprintf('%d%%',(i/length(data_rho))*100))
end;
close(h)

lower_predict(lower_predict < 0) = 0;

x = 0:71;
figure();
hold on
XDates = [datetime(2020,3,9:31) datetime(2020,4,1:30) datetime(2020,5,1:19)];
plot(XDates, mean_betat_predict, 'k')
plot(XDates, lower_predict, '--k')
plot(XDates, upper_predict, '--k')

patch([x, fliplr(x)], [lower_predict,fliplr(upper_predict)],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Beta(t)')
ylabel('Beta(t)')
xlabel('days')
xlim([datetime('9-Mar-2020') datetime('20-May-2020')]);

legend('Est. Mean Beta(t)','95% Prediction interval for Beta(t)')

savefig("post_predictive_solution_Betat_.fig");
close

%%
%Final version of CHI
%new daily identified infection: cases / new hidden infection I(t) 
mean_data_predict = (1/(length(post_samples)))*sum(data_predict');
mean_I_predict = mean(I_predict');

d = [datetime(2020,3,23) datetime(2020,3,26) datetime(2020,3,29) datetime(2020,4,1) datetime(2020,4,4) datetime(2020,4,7) datetime(2020,4,10) datetime(2020,4,13) datetime(2020,4,16) datetime(2020,4,19) datetime(2020,4,22) datetime(2020,4,25) datetime(2020,4,28) datetime(2020,5,1) datetime(2020,5,4) datetime(2020,5,7) datetime(2020,5,10) datetime(2020,5,13) datetime(2020,5,16) datetime(2020,5,19)];
y1 = [mean_data_predict(15) mean_I_predict(15); mean_data_predict(18) mean_I_predict(18); mean_data_predict(21) mean_I_predict(21); mean_data_predict(24) mean_I_predict(24); mean_data_predict(27) mean_I_predict(27); mean_data_predict(30) mean_I_predict(30) ; mean_data_predict(33) mean_I_predict(33); mean_data_predict(36) mean_I_predict(36); mean_data_predict(39) mean_I_predict(39); mean_data_predict(42) mean_I_predict(42); mean_data_predict(45) mean_I_predict(45); mean_data_predict(48) mean_I_predict(48); mean_data_predict(51) mean_I_predict(51); mean_data_predict(54) mean_I_predict(54); mean_data_predict(57) mean_I_predict(57); mean_data_predict(60) mean_I_predict(60); mean_data_predict(63) mean_I_predict(63); mean_data_predict(66) mean_I_predict(66); mean_data_predict(69) mean_I_predict(69); mean_data_predict(72) mean_I_predict(72)];

%legend('hidden infections','new identified infections')
 
days = 0:1:20;
y2 = [2797 1155 1926 3784 1851 1488 2761 2097 3092 4326 3943 4245 4335 4075 2622 3538 3226 4519 3896 2547];

figure
yyaxis left
bar(d,y1,'stacked');
colororder('default')
yyaxis right
p = plot(d,y2,'*-black');

legend('new identified infections','hidden infections')
title('New identified infections vs. hidden infections')
xlabel('Day')
yyaxis left
ylabel('Number of infections')
yyaxis right
ylabel('Number of testings','Color','black')
ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'k';

p.LineWidth = 3;
b.FaceColor = [ 0 0.447 0.741];

%%
%Final version of CCHI
f_rho_I_predict = f_rho_predict.* I_predict; 

mean_cum_case = cumtrapz(mean(f_rho_I_predict')); %total identified I
mean_cum_I = cumtrapz(mean(BIS_predict')); %total I

mean_cum_hI = mean_cum_I - mean_cum_case ;

d = [datetime(2020,3,23) datetime(2020,3,26) datetime(2020,3,29) datetime(2020,4,1) datetime(2020,4,4) datetime(2020,4,7) datetime(2020,4,10) datetime(2020,4,13) datetime(2020,4,16) datetime(2020,4,19) datetime(2020,4,22) datetime(2020,4,25) datetime(2020,4,28) datetime(2020,5,1) datetime(2020,5,4) datetime(2020,5,7) datetime(2020,5,10) datetime(2020,5,13) datetime(2020,5,16) datetime(2020,5,19)];
y = [mean_cum_case(15) mean_cum_hI(15); mean_cum_case(18) mean_cum_hI(18); mean_cum_case(21) mean_cum_hI(21); mean_cum_case(24) mean_cum_hI(24); mean_cum_case(27) mean_cum_hI(27); mean_cum_case(30) mean_cum_hI(30) ; mean_cum_case(33) mean_cum_hI(33); mean_cum_case(36) mean_cum_hI(36); mean_cum_case(39) mean_cum_hI(39); mean_cum_case(42) mean_cum_hI(42); mean_cum_case(45) mean_cum_hI(45); mean_cum_case(48) mean_cum_hI(48); mean_cum_case(51) mean_cum_hI(51); mean_cum_case(54) mean_cum_hI(54); mean_cum_case(57) mean_cum_hI(57); mean_cum_case(60) mean_cum_hI(60); mean_cum_case(63) mean_cum_hI(63); mean_cum_case(66) mean_cum_hI(66); mean_cum_case(69) mean_cum_hI(69); mean_cum_case(72) mean_cum_hI(72)];
h = bar(d,y,'stacked');
legend('Cumulative hidden infections','Cumulative identified infections')
title('Cumulative identified infections vs. hidden infections')
xlabel('Day')
ylabel('Number of infections')
%%
% CFR from April 6
mean_CFR_predict = mean(CFR_predict');

mean_CFR_predict = mean_CFR_predict(29:72);

%Generate 95% prediction interval
lower_predict = zeros(1, length(data_case));

upper_predict = zeros(1, length(data_case));

h = waitbar(0,'Initialize...');
for i = 1:length(data_case)
    
    lower_predict(1, i) = prctile(CFR_predict(i,:),2.5);
    upper_predict(1, i) = prctile(CFR_predict(i,:),97.5);
    
    
    waitbar(i/length(data_case),h,sprintf('%d%%',(i/length(data_case))*100))
end;
close(h)

lower_predict(lower_predict < 0) = 0;

lower_predict = lower_predict(29:72);
upper_predict = upper_predict(29:72);

Actual_cfr = Aver_death_all_report./cumsum(data_case);
Actual_cfr = Actual_cfr(29:72);

x = 0:43;
figure();
hold on
XDates = [datetime(2020,4,6:30) datetime(2020,5,1:19)];
plot(XDates,mean_CFR_predict,'k');
plot(XDates,lower_predict,"--k");
plot(XDates,Actual_cfr,'ob','MarkerSize',6);
plot(XDates,upper_predict,"--k");

patch([x, fliplr(x)], [lower_predict,fliplr(upper_predict)],  'k', 'FaceAlpha',0.2, 'EdgeColor','none')
hold off

title('Case Fatality Ratio')
ylabel('Case Fatality Ratio')
xlabel('days')
xlim([datetime('6-Apr-2020') datetime('20-May-2020')]);
ylim([0 0.3])
a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
% Create a vector of '%' signs
pct = char(ones(size(a,1),1)*'%'); 
% Append the '%' signs after the percentage values
new_yticks = [char(a),pct];
% 'Reflect the changes on the plot
set(gca,'yticklabel',new_yticks)

legend('Est. Mean Case Fatality Ratio','95% Prediction interval for Case Fatality Ratio','Actual Case infection ratio')

savefig("post_predictive_solution_Case Fatality Ratio_.fig");

%%
% R0 with 95%CI
% R0 = beta_0 * S_0 / (f_rho_0 + gamma_0)

beta_0 = theta_samples(1,:);
I_0 = theta_samples(4,:);
Pop_arr = zeros(1,length(l_post_samples)) + Pop;
S_0 = Pop_arr - I_0;
f_rho_0 = f_rho_predict(1,:);
gamma_0 = theta_samples(5,:);

R_0 = (beta_0.*S_0)./(f_rho_0 + gamma_0);

[M, I] = max(exp(l_post_samples));

beta0_best = theta_samples(1,I);
%6.1946e-08
S_0_best = Pop - theta_samples(4,I);
%4.3710e+06
f_best = theta_samples(2,I);
%1.1013
tim1_best = theta_samples(9,I);
%3.2212
gamma_best=theta_samples(5,I);
%0.1429
R_0_best = (beta0_best*S_0_best)./(f_best*min_rho + gamma_best);%1.8861
prctile(R_0, 2.5)%1.5663
prctile(R_0, 97.5)%2.0029%%