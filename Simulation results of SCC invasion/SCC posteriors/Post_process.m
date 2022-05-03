% Post_process.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of posterior densities of the
% parameter samples obtained by applying ABC scheme on the reference 
% SCC invasion dataset.

% Environment settings
clc
clear all
close all
 
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',4)
set(0,'defaultTextInterpreter','latex')

%% Read in all the final parameter estimates
paras_d3_all = readtable("Round 11 parameters run 3 day 3.txt");
paras_d6_all = readtable("Round 6 parameters run 2 day 6.txt");
paras_d9_all = readtable("Round 6 parameters run 1 day 9.txt");
paras_d12_all = readtable("Round 6 parameters run 1 day 12.txt");
paras_d14_all = readtable("Round 5 parameters run 2 day 14.txt");


%% dn posteriors
dn_d3_post = table2array(paras_d3_all(:,2));
dn_d6_post = table2array(paras_d6_all(:,2));
dn_d9_post = table2array(paras_d9_all(:,2));
dn_d12_post = table2array(paras_d12_all(:,2));
dn_d14_post = table2array(paras_d14_all(:,2));

[g,xii] = ksdensity(dn_d3_post,'Bandwidth',0.000006);
[h,xiii] = ksdensity(dn_d6_post,'Bandwidth',0.00006);
[j,xiv] = ksdensity(dn_d9_post,'Bandwidth',0.0004);
[z,xv] = ksdensity(dn_d12_post,'Bandwidth',0.0004);
[m,xvi] = ksdensity(dn_d14_post,'Bandwidth',0.0007);

figure
plot(xii, g, 'k-')
xlim([0.000069, 0.0003]);
ylim([0, 32000])
hold on;
plot(mean(dn_d3_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{d_{n}}$)')
ylabel('Probability density')
title({'Post-day 3', '$d_{n}$ posterior'})

figure
plot(xiii, h, 'm-')
xlim([0.000069, 0.0015]);
ylim([0, 3200])
hold on;
plot(mean(dn_d6_post), 0, 'xm', 'markersize', 20)
hold off;
xlabel('mean ($\hat{d_{n}}$)')
ylabel('Probability density')
title({'Post-day 6', '$d_{n}$ posterior'})

figure
plot(xiv, j, 'b-')
xlim([0.000069, 0.005]);
ylim([0, 700])
hold on;
plot(mean(dn_d9_post), 0, 'xb', 'markersize', 20)
hold off;
xlabel('mean ($\hat{d_{n}}$)')
ylabel('Probability density')
title({'Post-day 9', '$d_{n}$ posterior'})

figure
plot(xv, z, 'r-')
xlim([0.000069, 0.005]);
ylim([0, 700])
hold on;
plot(mean(dn_d12_post), 0, 'xr', 'markersize', 20)
hold off;
xlabel('mean ($\hat{d_{n}}$)')
ylabel('Probability density')
title({'Post-day 12', '$d_{n}$ posterior'})

figure
plot(xvi, m, 'g-')
xlim([0.000069, 0.02]);
ylim([0, 250])
hold on;
plot(mean(dn_d14_post), 0, 'xg', 'markersize', 20)
hold off;
xlabel('mean ($\hat{d_{n}}$)')
ylabel('Probability density')
title({'Post-day 14', '$d_{n}$ posterior'})

%% gamma posteriors
gamma_d3_post = table2array(paras_d3_all(:,3));
gamma_d6_post = table2array(paras_d6_all(:,3));
gamma_d9_post = table2array(paras_d9_all(:,3));
gamma_d12_post = table2array(paras_d12_all(:,3));
gamma_d14_post = table2array(paras_d14_all(:,3));

[g,xii] = ksdensity(gamma_d3_post,'Bandwidth',0.0002);
[h,xiii] = ksdensity(gamma_d6_post,'Bandwidth',0.03);
[j,xiv] = ksdensity(gamma_d9_post,'Bandwidth',0.04);
[z,xv] = ksdensity(gamma_d12_post,'Bandwidth',0.04);
[m,xvi] = ksdensity(gamma_d14_post,'Bandwidth',0.03);

figure
plot(xii, g, 'k-')
xlim([0.005, 0.007]);
ylim([0, 3500])
hold on;
plot(mean(gamma_d3_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{\gamma}$)')
ylabel('Probability density')
title({'Post-day 3', '$\gamma$ posterior'})

%figure
%plot(xiii, h, 'k-')
%xlim([0.005, 0.26]);
%ylim([0, 5])
%hold on;
%plot(mean(gamma_d6_post), 0, 'xk', 'markersize', 20)
%hold off;
%xlabel('mean ($\hat{\gamma}$)')
%ylabel('Probability density')
%title({'Post-day 6', '$\gamma$ posterior'})

%figure
%plot(xiv, j, 'k-')
%xlim([0.005 0.26]);
%ylim([0, 5])
%hold on;
%plot(mean(gamma_d9_post), 0, 'xk', 'markersize', 20)
%hold off;
%xlabel('mean ($\hat{\gamma}$)')
%ylabel('Probability density')
%title({'Post-day 9', '$\gamma$ posterior'})

%figure
%plot(xv, z, 'k-')
%xlim([0.005 0.26]);
%ylim([0, 6])
%hold on;
%plot(mean(gamma_d12_post), 0, 'xk', 'markersize', 20)
%hold off;
%xlabel('mean ($\hat{d_{n}}$)')
%ylabel('Probability density')
%title({'Post-day 12', '$\gamma$ posterior'})

%figure
%plot(xvi, m, 'k-')
%xlim([0.005, 0.26]);
%ylim([0, 5.5])
%hold on;
%plot(mean(gamma_d14_post), 0, 'xk', 'markersize', 20)
%hold off;
%xlabel('mean ($\hat{\gamma}$)')
%ylabel('Probability density')
%title({'Post-day 14', '$\gamma$ posterior'})

figure
plot(xiii, h, 'm-',xiv, j, 'b-', xv, z, 'r-', xvi, m , 'g-')
xlim([0.005 0.26])
ylim([0 7.5])
hold on;
plot(mean(gamma_d6_post), 0, 'xm', 'markersize', 20)
plot(mean(gamma_d9_post), 0, 'xb', 'markersize', 20)
plot(mean(gamma_d12_post), 0, 'xr', 'markersize', 20)
plot(mean(gamma_d14_post), 0, 'xg', 'markersize', 20)
hold off;
xlabel('mean ($\hat{\gamma}$)')
ylabel('Probability density')
lgd = legend({'$\gamma$ posterior (Post-day 6)','$\gamma$ posterior (Post-day 9)','$\gamma$ posterior (Post-day 12)','$\gamma$ posterior (Post-day 14)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
title('$\gamma$ posteriors')

%% rn posteriors 
rn_d3_post = table2array(paras_d3_all(:,4));
rn_d6_post = table2array(paras_d6_all(:,4));
rn_d9_post = table2array(paras_d9_all(:,4));
rn_d12_post = table2array(paras_d12_all(:,4));
rn_d14_post = table2array(paras_d14_all(:,4));

[g,xii] = ksdensity(rn_d3_post,'Bandwidth',0.001);
[h,xiii] = ksdensity(rn_d6_post,'Bandwidth',0.01);
[j,xiv] = ksdensity(rn_d9_post,'Bandwidth',0.015);
[z,xv] = ksdensity(rn_d12_post,'Bandwidth',0.012);
[m,xvi] = ksdensity(rn_d14_post,'Bandwidth',0.01);

figure
plot(xii, g, 'k-')
xlim([0.015, 0.03]);
ylim([0, 250])
hold on;
plot(mean(rn_d3_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{r_{n}}$)')
ylabel('Probability density')
title({'Post-day 3', '$r_{n}$ posterior'})

figure
plot(xiii, h, 'm-',xiv, j, 'b-', xv, z, 'r-', xvi, m , 'g-')
xlim([0.0008 0.08])
ylim([0 25])
hold on;
plot(mean(rn_d6_post), 0, 'xm', 'markersize', 20)
plot(mean(rn_d9_post), 0, 'xb', 'markersize', 20)
plot(mean(rn_d12_post), 0, 'xr', 'markersize', 20)
plot(mean(rn_d14_post), 0, 'xg', 'markersize', 20)
hold off;
xlabel('mean ($\hat{r_{n}}$)')
ylabel('Probability density')
lgd = legend({'$r_{n}$ posterior (Post-day 6)','$r_{n}$ posterior (Post-day 9)','$r_{n}$ posterior (Post-day 12)','$r_{n}$ posterior (Post-day 14)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
title('$r_{n}$ posteriors')

%% eta posteriors
eta_d3_post = table2array(paras_d3_all(:,5));
eta_d6_post = table2array(paras_d6_all(:,5));
eta_d9_post = table2array(paras_d9_all(:,5));
eta_d12_post = table2array(paras_d12_all(:,5));
eta_d14_post = table2array(paras_d14_all(:,5));

[g,xii] = ksdensity(eta_d3_post,'Bandwidth',0.12);
[h,xiii] = ksdensity(eta_d6_post,'Bandwidth',1.5);
[j,xiv] = ksdensity(eta_d9_post,'Bandwidth',1.5);
[z,xv] = ksdensity(eta_d12_post,'Bandwidth',1.5);
[m,xvi] = ksdensity(eta_d14_post,'Bandwidth',1.5);

figure
plot(xii, g, 'k-')
xlim([17.3, 18]);
ylim([0, 3])
hold on;
plot(mean(eta_d3_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{\eta}$)')
ylabel('Probability density')
title({'Post-day 3', '$\eta$ posterior'})

figure
plot(xiii, h, 'm-',xiv, j, 'b-', xv, z, 'r-', xvi, m , 'g-')
xlim([7 18])
ylim([0 0.2])
hold on;
plot(mean(eta_d6_post), 0, 'xm', 'markersize', 20)
plot(mean(eta_d9_post), 0, 'xb', 'markersize', 20)
plot(mean(eta_d12_post), 0, 'xr', 'markersize', 20)
plot(mean(eta_d14_post), 0, 'xg', 'markersize', 20)
hold off;
xlabel('mean ($\hat{\eta}$)')
ylabel('Probability density')
lgd = legend({'$\eta$ posterior (Post-day 6)','$\eta$ posterior (Post-day 9)','$\eta$ posterior (Post-day 12)','$\eta$ posterior (Post-day 14)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
title('$\eta$ posteriors')

%% dm posteriors
dm_d3_post = table2array(paras_d3_all(:,6));
dm_d6_post = table2array(paras_d6_all(:,6));
dm_d9_post = table2array(paras_d9_all(:,6));
dm_d12_post = table2array(paras_d12_all(:,6));
dm_d14_post = table2array(paras_d14_all(:,6));

[g,xii] = ksdensity(dm_d3_post,'Bandwidth',0.005);
[h,xiii] = ksdensity(dm_d6_post,'Bandwidth',0.003);
[j,xiv] = ksdensity(dm_d9_post,'Bandwidth',0.005);
[z,xv] = ksdensity(dm_d12_post,'Bandwidth',0.005);
[m,xvi] = ksdensity(dm_d14_post,'Bandwidth',0.005);

figure
plot(xii, g, 'k-')
xlim([0.0001 0.033]);
%ylim([0, 600])
hold on;
plot(mean(dm_d3_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{d_{m}}$)')
ylabel('Probability density')
title({'Post-day 3', '$d_{m}$ posterior'})

figure
plot(xiii, h, 'm-',xiv, j, 'b-', xv, z, 'r-', xvi, m , 'g-')
xlim([0.0001 0.033])
ylim([0 80])
hold on;
plot(mean(dm_d6_post), 0, 'xm', 'markersize', 20)
plot(mean(dm_d9_post), 0, 'xb', 'markersize', 20)
plot(mean(dm_d12_post), 0, 'xr', 'markersize', 20)
plot(mean(dm_d14_post), 0, 'xg', 'markersize', 20)
hold off;
xlabel('mean ($\hat{d_{m}}$)')
ylabel('Probability density')
lgd = legend({'$d_{m}$ posterior (Post-day 6)','$d_{m}$ posterior (Post-day 9)','$d_{m}$ posterior (Post-day 12)','$d_{m}$ posterior (Post-day 14)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
title('$d_{m}$ posteriors')

%% alpha posteriors
alpha_d3_post = table2array(paras_d3_all(:,7));
alpha_d6_post = table2array(paras_d6_all(:,7));
alpha_d9_post = table2array(paras_d9_all(:,7));
alpha_d12_post = table2array(paras_d12_all(:,7));
alpha_d14_post = table2array(paras_d14_all(:,7));

[g,xii] = ksdensity(alpha_d3_post,'Bandwidth',0.005);
[h,xiii] = ksdensity(alpha_d6_post,'Bandwidth',0.02);
[j,xiv] = ksdensity(alpha_d9_post,'Bandwidth',0.012);
[z,xv] = ksdensity(alpha_d12_post,'Bandwidth',0.01);
[m,xvi] = ksdensity(alpha_d14_post,'Bandwidth',0.01);

figure
plot(xii, g, 'k-')
xlim([0.16 0.18]);
ylim([0, 130])
hold on;
plot(mean(alpha_d3_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{\alpha}$)')
ylabel('Probability density')
title({'Post-day 3', '$\alpha$ posterior'})

figure
plot(xiii, h, 'm-',xiv, j, 'b-', xv, z, 'r-', xvi, m , 'g-')
xlim([0.07 0.18])
ylim([0 20])
hold on;
plot(mean(alpha_d6_post), 0, 'xm', 'markersize', 20)
plot(mean(alpha_d9_post), 0, 'xb', 'markersize', 20)
plot(mean(alpha_d12_post), 0, 'xr', 'markersize', 20)
plot(mean(alpha_d14_post), 0, 'xg', 'markersize', 20)
hold off;
xlabel('mean ($\hat{\alpha}$)')
ylabel('Probability density')
lgd = legend({'$\alpha$ posterior (Post-day 6)','$\alpha$ posterior (Post-day 9)','$\alpha$ posterior (Post-day 12)','$\alpha$ posterior (Post-day 14)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
title('$\alpha$ posteriors')

%% r.init posterior
r_init_post = table2array(paras_d3_all(:,8));
[g,xii] = ksdensity(r_init_post,'Bandwidth',0.1);

figure
plot(xii, g, 'k-')
xlim([1 5]);
ylim([0 5])
hold on;
plot(mean(r_init_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{R_{init.}}$)')
ylabel('Probability density')
title({'$R_{init.}$ posterior'})

%% p.ext posterior
p_ext_d3_post = table2array(paras_d3_all(:,9));
p_ext_d6_post = table2array(paras_d6_all(:,8));
p_ext_d9_post = table2array(paras_d9_all(:,8));
p_ext_d12_post = table2array(paras_d12_all(:,8));
p_ext_d14_post = table2array(paras_d14_all(:,8));

[g,xii] = ksdensity(p_ext_d3_post,'Bandwidth',0.005);
[h,xiii] = ksdensity(p_ext_d6_post,'Bandwidth',0.013);
[j,xiv] = ksdensity(p_ext_d9_post,'Bandwidth',0.013);
[z,xv] = ksdensity(p_ext_d12_post,'Bandwidth',0.013);
[m,xvi] = ksdensity(p_ext_d14_post,'Bandwidth',0.01);

figure
plot(xii, g, 'k-')
xlim([0.01 0.1]);
%ylim([0, 200])
hold on;
plot(mean(p_ext_d3_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{P_{ext.}}$)')
ylabel('Probability density')
title({'Post-day 3', '$P_{ext.}$ posterior'})

figure
plot(xiii, h, 'm-',xiv, j, 'b-', xv, z, 'r-', xvi, m , 'g-')
xlim([0.01 0.1])
ylim([0 23])
hold on;
plot(mean(p_ext_d6_post), 0, 'xm', 'markersize', 20)
plot(mean(p_ext_d9_post), 0, 'xb', 'markersize', 20)
plot(mean(p_ext_d12_post), 0, 'xr', 'markersize', 20)
plot(mean(p_ext_d14_post), 0, 'xg', 'markersize', 20)
hold off;
xlabel('mean ($\hat{P_{ext.}}$)')
ylabel('Probability density')
lgd = legend({'$P_{ext.}$ posterior (Post-day 6)','$P_{ext.}$ posterior (Post-day 9)','$P_{ext.}$ posterior (Post-day 12)','$P_{ext.}$ posterior (Post-day 14)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
title('$P_{ext.}$ posteriors')

%% p.mit. posteriors
p_mit_d3_post = table2array(paras_d3_all(:,10));
p_mit_d6_post = table2array(paras_d6_all(:,9));
p_mit_d9_post = table2array(paras_d9_all(:,9));
p_mit_d12_post = table2array(paras_d12_all(:,9));
p_mit_d14_post = table2array(paras_d14_all(:,9));

[g,xii] = ksdensity(p_mit_d3_post,'Bandwidth',0.03);
[h,xiii] = ksdensity(p_mit_d6_post,'Bandwidth',0.05);
[j,xiv] = ksdensity(p_mit_d9_post,'Bandwidth',0.06);
[z,xv] = ksdensity(p_mit_d12_post,'Bandwidth',0.03);
[m,xvi] = ksdensity(p_mit_d14_post,'Bandwidth',0.03);

figure
plot(xii, g, 'k-')
xlim([0.2 1]);
ylim([0 13])
hold on;
plot(mean(p_mit_d3_post), 0, 'xk', 'markersize', 20)
hold off;
xlabel('mean ($\hat{P_{mit.}}$)')
ylabel('Probability density')
title({'Post-day 3', '$P_{mit.}$ posterior'})

figure
plot(xiii, h, 'm-',xiv, j, 'b-', xv, z, 'r-', xvi, m , 'g-')
xlim([0.2 1])
ylim([0 18])
hold on;
plot(mean(p_mit_d6_post), 0, 'xm', 'markersize', 20)
plot(mean(p_mit_d9_post), 0, 'xb', 'markersize', 20)
plot(mean(p_mit_d12_post), 0, 'xr', 'markersize', 20)
plot(mean(p_mit_d14_post), 0, 'xg', 'markersize', 20)
hold off;
xlabel('mean ($\hat{P_{mit.}}$)')
ylabel('Probability density')
lgd = legend({'$P_{mit.}$ posterior (Post-day 6)','$P_{mit.}$ posterior (Post-day 9)','$P_{mit.}$ posterior (Post-day 12)','$P_{mit.}$ posterior (Post-day 14)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
title('$P_{mit.}$ posteriors')
