% Parameter_estimations.m
% Author: Yunchen Xiao

% This MATLAB file generates the plots of parameter estimates at
% different periods during the 14 days invasion period of SCC in 
% organotypic assay and the regression lines fitted to them.

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
paras_ests = readtable('Parameter estimates.txt');
dn_fitted_vals = readtable('Estimated values dn.txt');
gamma_fitted_vals = readtable('Estimated values gamma.txt');
rn_fitted_vals = readtable('Estimated values rn.txt');
eta_fitted_vals = readtable('Estimated values eta.txt');
dm_fitted_vals = readtable('Estimated values dm.txt');
alpha_fitted_vals = readtable('Estimated values alpha.txt');
p_ext_fitted_vals = readtable('Estimated values Pext.txt');
p_mit_fitted_vals = readtable('Estimated values Pmit.txt');

%% Estimates
% dn estimates
dn_ests = table2array(paras_ests(:,2));
dn_fitted_inter = table2array(dn_fitted_vals(:,2));
dn_fitted_linear = table2array(dn_fitted_vals(:,3));
dn_fitted_quad = table2array(dn_fitted_vals(:,4));

% gamma estimates
gamma_ests = table2array(paras_ests(:,3));
gamma_fitted_inter = table2array(gamma_fitted_vals(:,2));
gamma_fitted_linear = table2array(gamma_fitted_vals(:,3));
gamma_fitted_quad = table2array(gamma_fitted_vals(:,4));

% rn estimates
rn_ests = table2array(paras_ests(:,4));
rn_fitted_inter = table2array(rn_fitted_vals(:,2));
rn_fitted_linear = table2array(rn_fitted_vals(:,3));
rn_fitted_quad = table2array(rn_fitted_vals(:,4));

% eta estimates
eta_ests = table2array(paras_ests(:,5));
eta_fitted_inter = table2array(eta_fitted_vals(:,2));
eta_fitted_linear = table2array(eta_fitted_vals(:,3));
eta_fitted_quad = table2array(eta_fitted_vals(:,4));

% dm estimates
dm_ests = table2array(paras_ests(:,6));
dm_fitted_inter = table2array(dm_fitted_vals(:,2));
dm_fitted_linear = table2array(dm_fitted_vals(:,3));
dm_fitted_quad = table2array(dm_fitted_vals(:,4));

% alpha estimates
alpha_ests = table2array(paras_ests(:,7));
alpha_fitted_inter = table2array(alpha_fitted_vals(:,2));
alpha_fitted_linear = table2array(alpha_fitted_vals(:,3));
alpha_fitted_quad = table2array(alpha_fitted_vals(:,4));

% P.ext estimates
p_ext_ests = table2array(paras_ests(:,8));
p_ext_fitted_inter = table2array(p_ext_fitted_vals(:,2));
p_ext_fitted_linear = table2array(p_ext_fitted_vals(:,3));
p_ext_fitted_quad = table2array(p_ext_fitted_vals(:,4));

% P.mit estimates
p_mit_ests = table2array(paras_ests(:,9));
p_mit_fitted_inter = table2array(p_mit_fitted_vals(:,2));
p_mit_fitted_linear = table2array(p_mit_fitted_vals(:,3));
p_mit_fitted_quad = table2array(p_mit_fitted_vals(:,4));

%
x = 0:0.01:6;                        

%% Plot of dn estimates at different periods.
figure
plot(x, dn_fitted_inter, 'b-','Linewidth',3)
title({'Fitted regression models', 'on $\hat{d_n}$'})
hold on;
plot(x, dn_fitted_linear, 'r-','Linewidth',3)
plot(x, dn_fitted_quad, 'k-','Linewidth',7)
plot(1, dn_ests(1), 'xk',  'markersize', 20)
plot(2, dn_ests(2), 'xk', 'markersize', 20)
plot(3, dn_ests(3), 'xk', 'markersize', 20)
plot(4, dn_ests(4), 'xk', 'markersize', 20)
plot(5, dn_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([-0.002 0.007])
ylabel('Mean ($\hat{d_n}$)')

%% Plot of gamma estimates at different periods
figure
plot(x, gamma_fitted_inter, 'b-','Linewidth',3)
title({'Fitted regression models','on $\hat{\gamma}$'})
hold on;
plot(x, gamma_fitted_linear, 'r-','Linewidth',3)
plot(x, gamma_fitted_quad, 'k-','Linewidth',7)
plot(1, gamma_ests(1), 'xk', 'markersize', 20)
plot(2, gamma_ests(2), 'xk', 'markersize', 20)
plot(3, gamma_ests(3), 'xk', 'markersize', 20)
plot(4, gamma_ests(4), 'xk', 'markersize', 20)
plot(5, gamma_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([-0.05 0.16])
ylabel('Mean ($\hat{\gamma}$)')

%% Plot of rn estimates at different periods
figure
plot(x, rn_fitted_inter, 'b-','Linewidth',3)
title({'Fitted regression models','on $\hat{r_{n}}$'})
hold on;
plot(x, rn_fitted_linear, 'r-','Linewidth',3)
plot(x, rn_fitted_quad, 'k-','Linewidth',7)
plot(1, rn_ests(1), 'xk', 'markersize', 20)
plot(2, rn_ests(2), 'xk', 'markersize', 20)
plot(3, rn_ests(3), 'xk', 'markersize', 20)
plot(4, rn_ests(4), 'xk', 'markersize', 20)
plot(5, rn_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0 0.06])
ylabel('Mean ($\hat{r_{n}}$)')

%% Plot of eta estimations at different periods
figure
plot(x, eta_fitted_inter, 'b-','Linewidth',3)
title({'Fitted regression models','on $\hat{\eta}$'})
hold on;
plot(x, eta_fitted_linear, 'r-','Linewidth',3)
plot(x, eta_fitted_quad, 'k-','Linewidth',7)
plot(1, eta_ests(1), 'xk', 'markersize', 20)
plot(2, eta_ests(2), 'xk', 'markersize', 20)
plot(3, eta_ests(3), 'xk', 'markersize', 20)
plot(4, eta_ests(4), 'xk', 'markersize', 20)
plot(5, eta_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([10 20])
ylabel('Mean ($\hat{\eta}$)')

%% Plot of dm estimations at different periods
figure
plot(x, dm_fitted_inter, 'b-','Linewidth',7)
title({'Fitted regression models','on $\hat{d_{m}}$'})
hold on;
plot(x, dm_fitted_linear, 'r-','Linewidth',3)
plot(x, dm_fitted_quad, 'k-','Linewidth',3)
plot(1, dm_ests(1), 'xk', 'markersize', 20)
plot(2, dm_ests(2), 'xk', 'markersize', 20)
plot(3, dm_ests(3), 'xk', 'markersize', 20)
plot(4, dm_ests(4), 'xk', 'markersize', 20)
plot(5, dm_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0.012 0.017])
ylabel('Mean ($\hat{d_{m}}$)')

%% Plot of alpha estimations at different periods
figure
plot(x, alpha_fitted_inter, 'b-','Linewidth',3)
title({'Fitted regression models','on $\hat{\alpha}$'})
hold on;
plot(x, alpha_fitted_linear, 'r-','Linewidth',3)
plot(x, alpha_fitted_quad, 'k-','Linewidth',7)
plot(1, alpha_ests(1), 'xk', 'markersize', 20)
plot(2, alpha_ests(2), 'xk', 'markersize', 20)
plot(3, alpha_ests(3), 'xk', 'markersize', 20)
plot(4, alpha_ests(4), 'xk', 'markersize', 20)
plot(5, alpha_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0.1 0.2])
ylabel('Mean ($\hat{\alpha}$)')

%% Plot of P_ext estimations at different periods
figure
plot(x, p_ext_fitted_inter, 'b-','Linewidth',7)
title({'Fitted regression models','on $\hat{P_{ext,}}$'})
hold on;
plot(x, p_ext_fitted_linear, 'r-','Linewidth',3)
plot(x, p_ext_fitted_quad, 'k-','Linewidth',3)
plot(1, p_ext_ests(1), 'xk', 'markersize', 20)
plot(2, p_ext_ests(2), 'xk', 'markersize', 20)
plot(3, p_ext_ests(3), 'xk', 'markersize', 20)
plot(4, p_ext_ests(4), 'xk', 'markersize', 20)
plot(5, p_ext_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0.04 0.07])
ylabel('Mean ($\hat{P_{ext.}}$)')

%% Plot of P_mit estimations at different periods
figure
plot(x, p_mit_fitted_inter, 'b-','Linewidth',3)
title({'Fitted regression models','on $\hat{P_{mit.}}$'})
hold on;
plot(x, p_mit_fitted_linear, 'r-','Linewidth',3)
plot(x, p_mit_fitted_quad, 'k-','Linewidth',7)
plot(1, p_mit_ests(1), 'xk', 'markersize', 20)
plot(2, p_mit_ests(2), 'xk', 'markersize', 20)
plot(3, p_mit_ests(3), 'xk', 'markersize', 20)
plot(4, p_mit_ests(4), 'xk', 'markersize', 20)
plot(5, p_mit_ests(5), 'xk', 'markersize', 20)
xlim([0.5 5.5])
xlabel('Periods of invasion')
ylim([0.25 0.9])
ylabel('Mean ($\hat{P_{mit.}}$)')