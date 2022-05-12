% Posterior_prior_comparisons.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of posterior probability densities 
% of the parameter estimates obtained using non-error calibrated
% ABC on the SCC invasion pattern dataset with time-dependent parameters
% introduced.

%% Environment settings
clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

%% Posterior parameters
post_r9 = readtable("Round 9 parameters.txt");
prior_r1 = readtable("Round 1 initial time varying parameters.txt");
post_reg_coefs_bounds = readtable("Time-dependent parameters bounds.txt");

%% Separate the final parameter samples 
dn_r1 = table2array(prior_r1(:,2));
gamma_r1 = table2array(prior_r1(:,3));
rn_r1 = table2array(prior_r1(:,4));
eta_r1 = table2array(prior_r1(:,5));
dm_r1 = table2array(prior_r1(:,6));
alpha_r1 = table2array(prior_r1(:,7));
r_init_r1 = table2array(prior_r1(:,8));
p_ext_r1 = table2array(prior_r1(:,9));
p_mit_r1 = table2array(prior_r1(:,10));

dn_quad_r1 = table2array(prior_r1(:,11));
dn_lin_r1 = table2array(prior_r1(:,12));
gamma_quad_r1 = table2array(prior_r1(:,13));
gamma_lin_r1 = table2array(prior_r1(:,14));
rn_quad_r1 = table2array(prior_r1(:,15));
rn_lin_r1 = table2array(prior_r1(:,16));
eta_quad_r1 = table2array(prior_r1(:,17));
eta_lin_r1 = table2array(prior_r1(:,18));
alpha_quad_r1 = table2array(prior_r1(:,19));
alpha_lin_r1 = table2array(prior_r1(:,20));
p_mit_quad_r1 = table2array(prior_r1(:,21));
p_mit_lin_r1 = table2array(prior_r1(:,22));

dn = table2array(post_r9(:,2));
gamma = table2array(post_r9(:,3));
rn = table2array(post_r9(:,4));
eta = table2array(post_r9(:,5));
dm = table2array(post_r9(:,6));
alpha = table2array(post_r9(:,7));
r_init = table2array(post_r9(:,8));
p_ext = table2array(post_r9(:,9));
p_mit = table2array(post_r9(:,10));

dn_quad = table2array(post_r9(:,11));
dn_lin = table2array(post_r9(:,12));
gamma_quad = table2array(post_r9(:,13));
gamma_lin = table2array(post_r9(:,14));
rn_quad = table2array(post_r9(:,15));
rn_lin = table2array(post_r9(:,16));
eta_quad = table2array(post_r9(:,17));
eta_lin = table2array(post_r9(:,18));
alpha_quad = table2array(post_r9(:,19));
alpha_lin = table2array(post_r9(:,20));
p_mit_quad = table2array(post_r9(:,21));
p_mit_lin = table2array(post_r9(:,22));

%% dn posteriors
[g,xii] = ksdensity(dn, 'bandwidth', 0.0005);
figure
yline(1/0.019931,'b-','Linewidth',3.5);
xlim([0.000069 0.015]);
hold on;
plot(xii,g,'k-')
plot(mean(dn),0,'xk','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$d_{n}$ values')
ylabel('Probability density')

%% gamma posteriors
[g,xii] = ksdensity(gamma, 'bandwidth', 0.002);
figure
yline(1/0.255,'b-','Linewidth',3.5);
xlim([0.005 0.05]);
hold on;
plot(xii,g,'k-')
plot(mean(gamma),0,'xk','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density','Posterior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$\gamma$ values')
ylabel('Probability density')

%% rn posteriors
[g,xii] = ksdensity(rn, 'bandwidth', 0.008);
figure
yline(1/0.0792,'b-','Linewidth',3.5);
xlim([0.0008 0.08]);
hold on;
plot(xii,g,'k-')
plot(mean(rn),0,'xk','markersize',20)
plot(mean(rn_r1),0,'xb','markersize', 20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean', 'Prior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$r_{n}$ values')
ylabel('Probability density')

%% eta posteriors
[g,xii] = ksdensity(eta, 'bandwidth', 0.6);
figure
yline(1/11,'b-','Linewidth',3.5);
xlim([7 18]);
hold on;
plot(xii,g,'k-')
plot(mean(eta),0,'xk','markersize',20)
plot(mean(eta_r1),0,'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean', 'Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$\eta$ values')
ylabel('Probability density')

%% dm posteriors
[g,xii] = ksdensity(dm, 'bandwidth', 0.0025);
figure
yline(1/0.0329,'b-','Linewidth',3.5);
xlim([0.0001 0.033]);
hold on;
plot(xii,g,'k-')
plot(mean(dm),0,'xk','markersize',20)
plot(mean(dm_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean', 'Prior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$d_{m}$ values')
ylabel('Probability density')

%% alpha posteriors
[g,xii] = ksdensity(alpha, 'bandwidth', 0.008);
figure
yline(1/0.11,'b-','Linewidth',3.5);
xlim([0.07 0.18]);
hold on;
plot(xii,g,'k-')
plot(mean(alpha),0,'xk','markersize',20)
plot(mean(alpha_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean', 'Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$\alpha$ values')
ylabel('Probability density')

%% r_init posteriors
[g,xii] = ksdensity(r_init, 'bandwidth', 0.3);
figure
yline(1/4,'b-','Linewidth',3.5);
xlim([1 5]);
hold on;
plot(xii,g,'k-')
plot(mean(r_init),0,'xk','markersize',20)
plot(mean(r_init_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$R_{init.}$ values')
ylabel('Probability density')

%% P_ext posteriors
[g,xii] = ksdensity(p_ext, 'bandwidth', 0.008);
figure
yline(1/0.09,'b-','Linewidth',3.5);
xlim([0.01 0.1]);
hold on;
plot(xii,g,'k-')
plot(mean(p_ext),0,'xk','markersize',20)
plot(mean(p_ext_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

%% P_mit posteriors
[g,xii] = ksdensity(p_mit, 'bandwidth', 0.04);
figure
yline(1/0.8,'b-','Linewidth',3.5);
xlim([0.2 1]);
hold on;
plot(xii,g,'k-')
plot(mean(p_mit),0,'xk','markersize',20)
plot(mean(p_mit_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$P_{mit.}$ values')
ylabel('Probability density')

%% dn_quad posteriors
dn_quad_lb = table2array(post_reg_coefs_bounds(1,2));
dn_quad_ub = table2array(post_reg_coefs_bounds(1,3));
[g,xii] = ksdensity(dn_quad, 'bandwidth', 0.00005);
figure
yline(1/(dn_quad_ub - dn_quad_lb),'b-','Linewidth',3.5);
xlim([dn_quad_lb dn_quad_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(dn_quad),0,'xk','markersize',20)
plot(mean(dn_quad_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$d_{n, quad}$ values')
ylabel('Probability density')

%% dn_lin posteriors
dn_lin_lb = table2array(post_reg_coefs_bounds(2,2));
dn_lin_ub = table2array(post_reg_coefs_bounds(2,3));
[g,xii] = ksdensity(dn_lin, 'bandwidth', 0.00015);
figure
yline(1/(dn_lin_ub - dn_lin_lb),'b-','Linewidth',3.5);
xlim([dn_lin_lb dn_lin_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(dn_lin),0,'xk','markersize',20)
plot(mean(dn_lin_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$d_{n, lin}$ values')
ylabel('Probability density')

%% gamma_quad posteriors
gamma_quad_lb = table2array(post_reg_coefs_bounds(3,2));
gamma_quad_ub = table2array(post_reg_coefs_bounds(3,3));
[g,xii] = ksdensity(gamma_quad, 'bandwidth', 0.001);
figure
yline(1/(gamma_quad_ub - gamma_quad_lb),'b-','Linewidth',3.5);
xlim([gamma_quad_lb -0.012]);
hold on;
plot(xii,g,'k-')
plot(mean(gamma_quad),0,'xk','markersize',20)
plot(mean(gamma_quad_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$\gamma_{ quad}$ values')
ylabel('Probability density')

%% gamma_lin posteriors
gamma_lin_lb = table2array(post_reg_coefs_bounds(4,2));
gamma_lin_ub = table2array(post_reg_coefs_bounds(4,3));
[g,xii] = ksdensity(gamma_lin, 'bandwidth', 0.015);
figure
yline(1/(gamma_lin_ub - gamma_lin_lb),'b-','Linewidth',3.5);
xlim([gamma_lin_lb gamma_lin_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(gamma_lin),0,'xk','markersize',20)
plot(mean(gamma_lin_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$\gamma_{ lin}$ values')
ylabel('Probability density')

%% rn_quad posteriors
rn_quad_lb = table2array(post_reg_coefs_bounds(5,2));
rn_quad_ub = table2array(post_reg_coefs_bounds(5,3));
[g,xii] = ksdensity(rn_quad, 'bandwidth', 0.0005);
figure
yline(1/(rn_quad_ub - rn_quad_lb),'b-','Linewidth',3.5);
xlim([rn_quad_lb rn_quad_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(rn_quad),0,'xk','markersize',20)
plot(mean(rn_quad_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$r_{n, quad}$ values')
ylabel('Probability density')

%% rn_lin posteriors
rn_lin_lb = table2array(post_reg_coefs_bounds(6,2));
rn_lin_ub = table2array(post_reg_coefs_bounds(6,3));
[g,xii] = ksdensity(rn_lin, 'bandwidth', 0.0025);
figure
yline(1/(rn_lin_ub - rn_lin_lb),'b-','Linewidth',3.5);
xlim([rn_lin_lb rn_lin_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(rn_lin),0,'xk','markersize',20)
plot(0.0262, 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$r_{n, lin}$ values')
ylabel('Probability density')

%% eta_quad posteriors 
eta_quad_lb = table2array(post_reg_coefs_bounds(7,2));
eta_quad_ub = table2array(post_reg_coefs_bounds(7,3));
[g,xii] = ksdensity(eta_quad, 'bandwidth', 0.07);
figure
yline(1/(eta_quad_ub - eta_quad_lb),'b-','Linewidth',3.5);
xlim([eta_quad_lb eta_quad_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(eta_quad),0,'xk','markersize',20)
plot(mean(eta_quad_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$\eta_{ quad}$ values')
ylabel('Probability density')

%% eta_lin posteriors
eta_lin_lb = table2array(post_reg_coefs_bounds(8,2));
eta_lin_ub = table2array(post_reg_coefs_bounds(8,3));
[g,xii] = ksdensity(eta_lin, 'bandwidth', 0.7);
figure
yline(1/(eta_lin_ub - eta_lin_lb),'b-','Linewidth',3.5);
xlim([eta_lin_lb eta_lin_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(eta_lin),0,'xk','markersize',20)
plot(mean(eta_lin_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$\eta_{ lin}$ values')
ylabel('Probability density')

%% alpha_quad posteriors
alpha_quad_lb = table2array(post_reg_coefs_bounds(9,2));
alpha_quad_ub = table2array(post_reg_coefs_bounds(9,3));
[g,xii] = ksdensity(alpha_quad, 'bandwidth', 0.0008);
figure
yline(1/(alpha_quad_ub - alpha_quad_lb),'b-','Linewidth',3.5);
xlim([alpha_quad_lb alpha_quad_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(alpha_quad),0,'xk','markersize',20)
plot(mean(alpha_quad_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$\alpha_{ quad}$ values')
ylabel('Probability density')

%% alpha_lin posteriors
alpha_lin_lb = table2array(post_reg_coefs_bounds(10,2));
alpha_lin_ub = table2array(post_reg_coefs_bounds(10,3));
[g,xii] = ksdensity(alpha_lin, 'bandwidth', 0.005);
figure
yline(1/(alpha_lin_ub - alpha_lin_lb),'b-','Linewidth',3.5);
xlim([alpha_lin_lb alpha_lin_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(alpha_lin),0,'xk','markersize',20)
plot(mean(alpha_lin_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$\alpha_{ lin}$ values')
ylabel('Probability density')

%% p_mit_quad posteriors
p_mit_quad_lb = table2array(post_reg_coefs_bounds(11,2));
p_mit_quad_ub = table2array(post_reg_coefs_bounds(11,3));
[g,xii] = ksdensity(p_mit_quad, 'bandwidth', 0.008);
figure
yline(1/(p_mit_quad_ub - p_mit_quad_lb),'b-','Linewidth',3.5);
xlim([-0.2 p_mit_quad_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(p_mit_quad),0,'xk','markersize',20)
plot(mean(p_mit_quad_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$P_{mit, quad}$ values')
ylabel('Probability density')

%% p_mit_lin posteriors
p_mit_lin_lb = table2array(post_reg_coefs_bounds(12,2));
p_mit_lin_ub = table2array(post_reg_coefs_bounds(12,3));
[g,xii] = ksdensity(p_mit_lin, 'bandwidth', 0.05);
figure
yline(1/(p_mit_lin_ub - p_mit_lin_lb),'b-','Linewidth',3.5);
xlim([p_mit_lin_lb p_mit_lin_ub]);
hold on;
plot(xii,g,'k-')
plot(mean(p_mit_lin),0,'xk','markersize',20)
plot(mean(p_mit_lin_r1), 0, 'xb','markersize',20)
hold off;
lgd = legend({'Posterior density','Prior density', 'Posterior mean','Prior mean'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$P_{mit, lin}$ values')
ylabel('Probability density')