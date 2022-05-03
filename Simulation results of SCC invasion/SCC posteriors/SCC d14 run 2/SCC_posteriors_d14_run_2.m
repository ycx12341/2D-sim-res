% SCC_posteriors_d14_run2.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of probability densities of the
% parameter estimations at the end of each round of the second run of 
% applying ABC scheme on the post-day 14 pattern of the SCC reference 
% dataset.

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
paras_init_d14 = readtable("Round 1 initial parameters.txt");
paras_r2_d14 = readtable("Round 2 parameters.txt");
paras_r3_d14 = readtable("Round 3 parameters.txt");
paras_r4_d14 = readtable("Round 4 parameters.txt");
paras_r5_d14 = readtable("Round 5 parameters.txt");

%% Density plot of dn
dn_init = table2array(paras_init_d14(:,2));
dn_r2 = table2array(paras_r2_d14(:,2));
dn_r3 = table2array(paras_r3_d14(:,2));
dn_r4 = table2array(paras_r4_d14(:,2));
dn_r5 = table2array(paras_r5_d14(:,2));

[g,xii] = ksdensity(dn_r2,'Bandwidth',0.0007);
[h,xiii] = ksdensity(dn_r3,'Bandwidth',0.0007);
[j,xiv] = ksdensity(dn_r4,'Bandwidth',0.0007);
[z,xv] = ksdensity(dn_r5,'Bandwidth',0.0007);

figure
yline(1/(0.02 - 0.000069),'b','Linewidth',3.5);
xlim([0.000069 0.02]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-')
plot(mean(dn_init),0,'xb','markersize',20)
plot(mean(dn_r2),0,'xg','markersize',20)
plot(mean(dn_r3),0,'xr','markersize',20)
plot(mean(dn_r4),0,'xc','markersize',20)
plot(mean(dn_r5),0,'xm','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n}$)','Post-round 1 density($d_{n}$)','Post-round 2 density($d_{n}$)','Post-round 3 density($d_{n}$)','Post-round 4 density($d_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$d_{n}$ values')
ylabel('Probability density')

%% Density plot of gamma
gamma_init = table2array(paras_init_d14(:,3));
gamma_r2 = table2array(paras_r2_d14(:,3));
gamma_r3 = table2array(paras_r3_d14(:,3));
gamma_r4 = table2array(paras_r4_d14(:,3));
gamma_r5 = table2array(paras_r5_d14(:,3));

[g,xii] = ksdensity(gamma_r2,'Bandwidth',0.03);
[h,xiii] = ksdensity(gamma_r3,'Bandwidth',0.03);
[j,xiv] = ksdensity(gamma_r4,'Bandwidth',0.03);
[z,xv] = ksdensity(gamma_r5,'Bandwidth',0.03);

figure
yline(1/(0.26 - 0.005),'b','Linewidth',3.5);
xlim([0.005 0.26]);
%ylim([0 6])
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-')
plot(mean(gamma_init),0,'xb','markersize',20)
plot(mean(gamma_r2),0,'xg','markersize',20)
plot(mean(gamma_r3),0,'xr','markersize',20)
plot(mean(gamma_r4),0,'xc','markersize',20)
plot(mean(gamma_r5),0,'xm','markersize',20)
hold off;
lgd = legend({'Initial density($\gamma$)','Post-round 1 density($\gamma$)','Post-round 2 density($\gamma$)','Post-round 3 density($\gamma$)','Post-round 4 density($\gamma$)'},'Interpreter','latex','Location','southwest','Orientation','vertical','Fontsize',10);
xlabel('$\gamma$ values')
ylabel('Probability density')

%% Density plot of rn
rn_init = table2array(paras_init_d14(:,4));
rn_r2 = table2array(paras_r2_d14(:,4));
rn_r3 = table2array(paras_r3_d14(:,4));
rn_r4 = table2array(paras_r4_d14(:,4));
rn_r5 = table2array(paras_r5_d14(:,4));

[g,xii] = ksdensity(rn_r2,'Bandwidth',0.01);
[h,xiii] = ksdensity(rn_r3,'Bandwidth',0.01);
[j,xiv] = ksdensity(rn_r4,'Bandwidth',0.01);
[z,xv] = ksdensity(rn_r5,'Bandwidth',0.01);

figure
yline(1/(0.08 - 0.0008),'b','Linewidth',3.5);
xlim([0.0008 0.08]);
ylim([0 17])
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-')
plot(mean(rn_init),0,'xb','markersize',20)
plot(mean(rn_r2),0,'xg','markersize',20)
plot(mean(rn_r3),0,'xr','markersize',20)
plot(mean(rn_r4),0,'xc','markersize',20)
plot(mean(rn_r5),0,'xm','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n}$)','Post-round 1 density($r_{n}$)','Post-round 2 density($r_{n}$)','Post-round 3 density($r_{n}$)','Post-round 4 density($r_{n}$)'},'Interpreter','latex','Location','southwest','Orientation','vertical','Fontsize',10);
xlabel('$r_{n}$ values')
ylabel('Probability density')

%% Density plots of eta
eta_init = table2array(paras_init_d14(:,5));
eta_r2 = table2array(paras_r2_d14(:,5));
eta_r3 = table2array(paras_r3_d14(:,5));
eta_r4 = table2array(paras_r4_d14(:,5));
eta_r5 = table2array(paras_r5_d14(:,5));

[g,xii] = ksdensity(eta_r2,'Bandwidth',1.5);
[h,xiii] = ksdensity(eta_r3,'Bandwidth',1.5);
[j,xiv] = ksdensity(eta_r4,'Bandwidth',1.5);
[z,xv] = ksdensity(eta_r5,'Bandwidth',1.5);

figure
yline(1/(18 - 7),'b','Linewidth',3.5);
xlim([7 18]);
ylim([0 0.15])
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-')
plot(mean(eta_init),0,'xb','markersize',20)
plot(mean(eta_r2),0,'xg','markersize',20)
plot(mean(eta_r3),0,'xr','markersize',20)
plot(mean(eta_r4),0,'xc','markersize',20)
plot(mean(eta_r5),0,'xm','markersize',20)
hold off;
lgd = legend({'Initial density($\eta$)','Post-round 1 density($\eta$)','Post-round 2 density($\eta$)','Post-round 3 density($\eta$)','Post-round 4 density($\eta$)'},'Interpreter','latex','Location','southwest','Orientation','vertical','Fontsize',10);
xlabel('$\eta$ values')
ylabel('Probability density')

%% Density plots of dm
dm_init = table2array(paras_init_d14(:,6));
dm_r2 = table2array(paras_r2_d14(:,6));
dm_r3 = table2array(paras_r3_d14(:,6));
dm_r4 = table2array(paras_r4_d14(:,6));
dm_r5 = table2array(paras_r5_d14(:,6));

[g,xii] = ksdensity(dm_r2,'Bandwidth',0.005);
[h,xiii] = ksdensity(dm_r3,'Bandwidth',0.005);
[j,xiv] = ksdensity(dm_r4,'Bandwidth',0.005);
[z,xv] = ksdensity(dm_r5,'Bandwidth',0.005);

figure
yline(1/(0.033 - 0.0001),'b','Linewidth',3.5);
xlim([0.0001 0.033]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-')
plot(mean(dm_init),0,'xb','markersize',20)
plot(mean(dm_r2),0,'xg','markersize',20)
plot(mean(dm_r3),0,'xr','markersize',20)
plot(mean(dm_r4),0,'xc','markersize',20)
plot(mean(dm_r5),0,'xm','markersize',20)
hold off;
lgd = legend({'Initial density($d_{m}$)','Post-round 1 density($d_{m}$)','Post-round 2 density($d_{m}$)','Post-round 3 density($d_{m}$)','Post-round 4 density($d_{m}$)'},'Interpreter','latex','Location','southwest','Orientation','vertical','Fontsize',10);
xlabel('$d_{m}$ values')
ylabel('Probability density')

%% Density plots of alpha
alpha_init = table2array(paras_init_d14(:,7));
alpha_r2 = table2array(paras_r2_d14(:,7));
alpha_r3 = table2array(paras_r3_d14(:,7));
alpha_r4 = table2array(paras_r4_d14(:,7));
alpha_r5 = table2array(paras_r5_d14(:,7));

[g,xii] = ksdensity(alpha_r2,'Bandwidth',0.01);
[h,xiii] = ksdensity(alpha_r3,'Bandwidth',0.01);
[j,xiv] = ksdensity(alpha_r4,'Bandwidth',0.01);
[z,xv] = ksdensity(alpha_r5,'Bandwidth',0.01);

figure
yline(1/(0.18 - 0.07),'b','Linewidth',3.5);
xlim([0.07 0.18]);
ylim([0 13]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-')
plot(mean(alpha_init),0,'xb','markersize',20)
plot(mean(alpha_r2),0,'xg','markersize',20)
plot(mean(alpha_r3),0,'xr','markersize',20)
plot(mean(alpha_r4),0,'xc','markersize',20)
plot(mean(alpha_r5),0,'xm','markersize',20)
hold off;
lgd = legend({'Initial density($\alpha$)','Post-round 1 density($\alpha$)','Post-round 2 density($\alpha$)','Post-round 3 density($\alpha$)','Post-round 4 density($\alpha$)'},'Interpreter','latex','Location','southeast','Orientation','vertical','Fontsize',10);
xlabel('$\alpha$ values')
ylabel('Probability density')

%% Density plots of p.ext
p_ext_init = table2array(paras_init_d14(:,8));
p_ext_r2 = table2array(paras_r2_d14(:,8));
p_ext_r3 = table2array(paras_r3_d14(:,8));
p_ext_r4 = table2array(paras_r4_d14(:,8));
p_ext_r5 = table2array(paras_r5_d14(:,8));

[g,xii] = ksdensity(p_ext_r2,'Bandwidth',0.01);
[h,xiii] = ksdensity(p_ext_r3,'Bandwidth',0.01);
[j,xiv] = ksdensity(p_ext_r4,'Bandwidth',0.01);
[z,xv] = ksdensity(p_ext_r5,'Bandwidth',0.01);

figure
yline(1/(0.1 - 0.01),'b','Linewidth',3.5);
xlim([0.01 0.1]);
%ylim([0 13]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-')
plot(mean(p_ext_init),0,'xb','markersize',20)
plot(mean(p_ext_r2),0,'xg','markersize',20)
plot(mean(p_ext_r3),0,'xr','markersize',20)
plot(mean(p_ext_r4),0,'xc','markersize',20)
plot(mean(p_ext_r5),0,'xm','markersize',20)
hold off;
lgd = legend({'Initial density($P_{ext.}$)','Post-round 1 density($P_{ext.}$)','Post-round 2 density($P_{ext.}$)','Post-round 3 density($P_{ext.}$)','Post-round 4 density($P_{ext.}$)'},'Interpreter','latex','Location','southeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

%% Density plots of p.mit
p_mit_init = table2array(paras_init_d14(:,9));
p_mit_r2 = table2array(paras_r2_d14(:,9));
p_mit_r3 = table2array(paras_r3_d14(:,9));
p_mit_r4 = table2array(paras_r4_d14(:,9));
p_mit_r5 = table2array(paras_r5_d14(:,9));

[g,xii] = ksdensity(p_mit_r2,'Bandwidth',0.03);
[h,xiii] = ksdensity(p_mit_r3,'Bandwidth',0.03);
[j,xiv] = ksdensity(p_mit_r4,'Bandwidth',0.03);
[z,xv] = ksdensity(p_mit_r5,'Bandwidth',0.03);

figure
yline(1/(1 - 0.2),'b','Linewidth',3.5);
xlim([0.2 1]);
%ylim([0 13]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-')
plot(mean(p_mit_init),0,'xb','markersize',20)
plot(mean(p_mit_r2),0,'xg','markersize',20)
plot(mean(p_mit_r3),0,'xr','markersize',20)
plot(mean(p_mit_r4),0,'xc','markersize',20)
plot(mean(p_mit_r5),0,'xm','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit.}$)','Post-round 1 density($P_{mit.}$)','Post-round 2 density($P_{mit.}$)','Post-round 3 density($P_{mit.}$)','Post-round 4 density($P_{mit.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit.}$ values')
ylabel('Probability density')
