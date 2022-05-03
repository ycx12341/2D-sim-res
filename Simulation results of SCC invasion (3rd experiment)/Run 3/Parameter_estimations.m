% Parameter_estimations.m
% Author: Yunchen Xiao

% This MATLAB file generates the plots of final parameter estimates at
% different time points during the 14 days invasion period of SCC in 
% organotypic assay (time-dependent parameters introduced to the model).

%% Environment settings
clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 18)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

%% Read in all the data.
paras_ests = readtable('Estimated final parameter values regression.txt');
paras_ests_final = readtable('Round 9 final full parameter estimates run 3.txt');

%% Estimates
% dn estimates
dn_ests = table2array(paras_ests_final(1,2:6));
dn_fitted_quad = table2array(paras_ests(:,2));

% gamma estimates
gamma_ests = table2array(paras_ests_final(2,2:6));
gamma_fitted_quad = table2array(paras_ests(:,3));

% rn estimates
rn_ests = table2array(paras_ests_final(3,2:6));
rn_fitted_quad = table2array(paras_ests(:,4));

% eta estimates
eta_ests = table2array(paras_ests_final(4,2:6));
eta_fitted_quad = table2array(paras_ests(:,5));

% dm estimates
dm_ests = table2array(paras_ests_final(5,2:6));
dm_fitted_quad = table2array(paras_ests(:,6));

% alpha estimates
alpha_ests = table2array(paras_ests_final(6,2:6));
alpha_fitted_quad = table2array(paras_ests(:,7));

% P.ext estimates
p_ext_ests = table2array(paras_ests_final(8,2:6));
p_ext_fitted_quad = table2array(paras_ests(:,9));

% P.mit estimates
p_mit_ests = table2array(paras_ests_final(9,2:6));
p_mit_fitted_quad = table2array(paras_ests(:,10));

%
x = 0:0.01:6;                        

%% Plot of dn estimates at different periods of the invasion.
figure
plot(x, dn_fitted_quad, 'k-','Linewidth',3)
title({'Fitted regression models', 'on $\hat{d_n}$ (time-dependent', 'parameters)'})
hold on;
plot(1, dn_ests(1), 'xk',  'markersize', 20)
plot(2, dn_ests(2), 'xk', 'markersize', 20)
plot(3, dn_ests(3), 'xk', 'markersize', 20)
plot(4, dn_ests(4), 'xk', 'markersize', 20)
plot(5, dn_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([-0.002 0.015])
ylabel('Mean ($\hat{d_n}$)')

%% Plot of gamma estimates at different periods of the invasion.
figure
plot(x, gamma_fitted_quad, 'k-','Linewidth',3)
title({'Fitted regression models','on $\hat{\gamma}$ (time-dependent parameters)'})
hold on;
plot(1, gamma_ests(1), 'xk', 'markersize', 20)
plot(2, gamma_ests(2), 'xk', 'markersize', 20)
plot(3, gamma_ests(3), 'xk', 'markersize', 20)
plot(4, gamma_ests(4), 'xk', 'markersize', 20)
plot(5, gamma_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([-0.05 0.16])
ylabel('Mean ($\hat{\gamma}$)')

%% Plot of rn estimates at different periods of the invasion. 
figure
plot(x, rn_fitted_quad, 'k-','Linewidth',3)
title({'Fitted regression models','on $\hat{r_{n}}$ (time-dependent parameters)'})
hold on;
plot(1, rn_ests(1), 'xk', 'markersize', 20)
plot(2, rn_ests(2), 'xk', 'markersize', 20)
plot(3, rn_ests(3), 'xk', 'markersize', 20)
plot(4, rn_ests(4), 'xk', 'markersize', 20)
plot(5, rn_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0 0.06])
ylabel('Mean ($\hat{r_{n}}$)')

%% Plot of eta estimations at different periods of invasion.
figure
plot(x, eta_fitted_quad, 'k-','Linewidth',3)
title({'Fitted regression models','on $\hat{\eta}$ (time-dependent parameters)'})
hold on;
plot(1, eta_ests(1), 'xk', 'markersize', 20)
plot(2, eta_ests(2), 'xk', 'markersize', 20)
plot(3, eta_ests(3), 'xk', 'markersize', 20)
plot(4, eta_ests(4), 'xk', 'markersize', 20)
plot(5, eta_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([10 20])
ylabel('Mean ($\hat{\eta}$)')

%% Plot of dm estimations at different periods of the invasion.
figure
plot(1, dm_ests(1), 'xk', 'markersize', 20)
title({'Fitted regression models','on $\hat{d_{m}}$ (time-dependent parameters)'})
hold on;
plot(2, dm_ests(2), 'xk', 'markersize', 20)
plot(3, dm_ests(3), 'xk', 'markersize', 20)
plot(4, dm_ests(4), 'xk', 'markersize', 20)
plot(5, dm_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0.015 0.017])
ylabel('Mean ($\hat{d_{m}}$)')

%% Plot of alpha estimations at different periods of the invasion.
figure
plot(x, alpha_fitted_quad, 'k-','Linewidth',3)
title({'Fitted regression models','on $\hat{\alpha}$ (time-dependent parameters)'})
hold on;
plot(1, alpha_ests(1), 'xk', 'markersize', 20)
plot(2, alpha_ests(2), 'xk', 'markersize', 20)
plot(3, alpha_ests(3), 'xk', 'markersize', 20)
plot(4, alpha_ests(4), 'xk', 'markersize', 20)
plot(5, alpha_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0.07 0.18])
ylabel('Mean ($\hat{\alpha}$)')

%% Plot of P_ext estimations at different periods of the invasion.
figure
plot(1, p_ext_ests(1), 'xk', 'markersize', 20)
title({'Fitted regression models','on $\hat{P_{ext,}}$ (time-dependent parameters)'})
hold on;
plot(2, p_ext_ests(2), 'xk', 'markersize', 20)
plot(3, p_ext_ests(3), 'xk', 'markersize', 20)
plot(4, p_ext_ests(4), 'xk', 'markersize', 20)
plot(5, p_ext_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0.08 0.10])
ylabel('Mean ($\hat{P_{ext.}}$)')

%% Plot of P_mit estimations at different periods of the invasion. 
figure
plot(x, p_mit_fitted_quad, 'k-','Linewidth',3)
title({'Fitted regression models','on $\hat{P_{mit.}}$ (time-dependent parameters)'})
hold on;
plot(1, p_mit_ests(1), 'xk', 'markersize', 20)
plot(2, p_mit_ests(2), 'xk', 'markersize', 20)
plot(3, p_mit_ests(3), 'xk', 'markersize', 20)
plot(4, p_mit_ests(4), 'xk', 'markersize', 20)
plot(5, p_mit_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0.2 0.9])
ylabel('Mean ($\hat{P_{mit.}}$)')