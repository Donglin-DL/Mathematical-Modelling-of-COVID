%PRCC for constant CIFR(take last value), infected propotion (S(T)-S(0))/S(0) 
%and case-infection ratio(take last value)

%PRCC for CIFR at last day
death = Aver_death_all_report;
load('trapz_BIS.mat')
load('theta_samples.mat')
n = length(trapz_BIS(1,:));
m = length(trapz_BIS(:,1));

cifr_sample = death(length(death))./trapz_BIS(length(death),:);
LHS_matrix = [theta_samples(1,:);theta_samples(6,:);theta_samples(7,:);theta_samples(15,:);theta_samples(16,:);theta_samples(2,:);theta_samples(9,:);theta_samples(10,:);theta_samples(11,:);theta_samples(12,:);theta_samples(13,:);theta_samples(14,:);theta_samples(17,:);theta_samples(18,:);theta_samples(4,:);theta_samples(5,:)].';
order = [1,6, 7 ,15,16,2,9,10,11,12,13,14,17,18,4,5];
PRCC_var_all = ["beta", "f" ,'p', 'I0', 'gamma', 'q', 's', 'p2', 'time1', 'time2', 'time3', 'time4', 'highval', 'constval', 'time1beta', 'time2beta', 't1b', 'm0val'];
PRCC_var = PRCC_var_all(order);
[prcc,sign,sign_label]=PRCC(LHS_matrix,cifr_sample,1,PRCC_var,0.05);
X=categorical(PRCC_var);
X=reordercats(X,PRCC_var);
h = figure;
p1=bar(X(1:5),prcc(1:5),'r');
hold on
p2=bar(X(6:14),prcc(6:14),'b');
hold on
p3=bar(X(15:16),prcc(15:16));
xlabel('Parameters')
ylabel('PRCC for CIFR')
%saveas(h,sprintf('PRCC all parameters for CIFR'));

%scatter plot
y_var1 = {"CIFR"};
%PRCC_PLOT(LHS_matrix,cifr_sample,1,PRCC_var,y_var1) 

%PRCC for infected propotion
load('S_predict.mat')
ipro_sample = [];
for i=1:n
    ipro_sample(i) = -(S_predict(m,i)-S_predict(1,i))/S_predict(1,i);
end
[prcc1,sign1,sign_label1]=PRCC(LHS_matrix,ipro_sample,1,PRCC_var,0.05);
X=categorical(PRCC_var);
X=reordercats(X,PRCC_var);
h = figure;
p1=bar(X(1:5),prcc1(1:5),'r');
hold on
p2=bar(X(6:14),prcc1(6:14),'b');
hold on
p3=bar(X(15:16),prcc1(15:16));
xlabel('Parameters')
ylabel('PRCC for Propotion of Infected')
%saveas(h,sprintf('PRCC all parameters for Propotion of Infected'));
y_var2 = {"Propotion of Infected"};
%PRCC_PLOT(LHS_matrix,ipro_sample,1,PRCC_var,y_var2) 

%PRCC for case_infection ratio constant: sum(cases)/sum(bis)
load('case_out.mat')
cir_sample = [];
for i=1:n
    cir_sample(i) = sum(case_out(:,i))/trapz_BIS(length(death),i);
end

[prcc2,sign2,sign_label2]=PRCC(LHS_matrix,cir_sample,1,PRCC_var,0.05);
X=categorical(PRCC_var);
X=reordercats(X,PRCC_var);
h = figure;
p1=bar(X(1:5),prcc2(1:5),'r');
hold on
p2=bar(X(6:14),prcc2(6:14),'b');
hold on
p3=bar(X(15:16),prcc2(15:16));
xlabel('Parameters')
ylabel('PRCC for Case-Infection Ratio')
%saveas(h,sprintf('PRCC all parameters for Case-Infection Ratio'));

y_var3 = {"Case_Infection Ratio"};
%PRCC_PLOT(LHS_matrix,cir_sample,1,PRCC_var,y_var3) 

%PRCC for peak time
peak_time = [];
for i=1:n
    [val, peak_time(i)] = max(I_predict(:,i));
end

[prcc3,sign3,sign_label3]=PRCC(LHS_matrix,peak_time,1,PRCC_var,0.05);
X=categorical(PRCC_var);
X=reordercats(X,PRCC_var);
h = figure;
p1=bar(X(1:5),prcc3(1:5),'r');
hold on
p2=bar(X(6:14),prcc3(6:14),'b');
hold on
p3=bar(X(15:16),prcc3(15:16));
xlabel('Parameters')
ylabel('PRCC for Peak Time')
saveas(h,sprintf('PRCC all parameters for Peak time case'));
ylim([-0.8 0.9])
y_var4 = {"Peak Time"};