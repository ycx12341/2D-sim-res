% SCC_posteriors_d3_run3.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of probability densities of the
% parameter estimations at the end of each round of the third run of 
% applying ABC scheme on the SCC invasion patterns with time-dependent
% parameter values introduced to the model.

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

%% Read in all the final parameter estimates
paras_init_tv = readtable("Round 1 initial time varying parameters.txt");
paras_r2_tv = readtable("Round 2 parameters.txt");
paras_r3_tv = readtable("Round 3 parameters.txt");
paras_r4_tv = readtable("Round 4 parameters.txt");
paras_r5_tv = readtable("Round 5 parameters.txt");
paras_r6_tv = readtable("Round 6 parameters.txt");
paras_r7_tv = readtable("Round 7 parameters.txt");
paras_r8_tv = readtable("Round 8 parameters.txt");
paras_r9_tv = readtable("Round 9 parameters.txt");

%% Regression coefficients bounds
post_reg_coefs_bounds = readtable("Time-dependent parameters bounds.txt");

%% Density plot of dn
dn_init = table2array(paras_init_tv(:,2));
dn_r2 = table2array(paras_r2_tv(:,2));
dn_r3 = table2array(paras_r3_tv(:,2));
dn_r4 = table2array(paras_r4_tv(:,2));
dn_r5 = table2array(paras_r5_tv(:,2));
dn_r6 = table2array(paras_r6_tv(:,2));
dn_r7 = table2array(paras_r7_tv(:,2));
dn_r8 = table2array(paras_r8_tv(:,2));
dn_r9 = table2array(paras_r9_tv(:,2));

[g,xii] = ksdensity(dn_r2,'Bandwidth',0.0005);
[h,xiii] = ksdensity(dn_r3,'Bandwidth',0.0005);
[j,xiv] = ksdensity(dn_r4,'Bandwidth',0.0005);
[z,xv] = ksdensity(dn_r5,'Bandwidth',0.0005);
[m,xvi] = ksdensity(dn_r6,'Bandwidth',0.0005);

figure
yline(1/(0.02 - 0.000069),'b','Linewidth',3.5);
xlim([0.000069 0.015]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(dn_init),0,'xb','markersize',20)
plot(mean(dn_r2),0,'xg','markersize',20)
plot(mean(dn_r3),0,'xr','markersize',20)
plot(mean(dn_r4),0,'xc','markersize',20)
plot(mean(dn_r5),0,'xm','markersize',20)
plot(mean(dn_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n}$)','Post-round 1 density($d_{n}$)','Post-round 2 density($d_{n}$)','Post-round 3 density($d_{n}$)','Post-round 4 density($d_{n}$)', 'Post-round 5 density($d_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$d_{n}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(dn_r7,'Bandwidth',0.0005);
[l,xviii] = ksdensity(dn_r8,'Bandwidth',0.0005);
[q,xix] = ksdensity(dn_r9,'Bandwidth',0.0005);

figure
yline(1/(0.02 - 0.000069),'b','Linewidth',3.5);
xlim([0.000069 0.015]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(dn_r7),0,'xg','markersize',20)
plot(mean(dn_r8),0,'xr','markersize',20)
plot(mean(dn_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n}$)','Post-round 6 density($d_{n}$)','Post-round 7 density($d_{n}$)','Post-round 8 density($d_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$d_{n}$ values')
ylabel('Probability density')

%% Density plot of gamma
gamma_init = table2array(paras_init_tv(:,3));
gamma_r2 = table2array(paras_r2_tv(:,3));
gamma_r3 = table2array(paras_r3_tv(:,3));
gamma_r4 = table2array(paras_r4_tv(:,3));
gamma_r5 = table2array(paras_r5_tv(:,3));
gamma_r6 = table2array(paras_r6_tv(:,3));
gamma_r7 = table2array(paras_r7_tv(:,3));
gamma_r8 = table2array(paras_r8_tv(:,3));
gamma_r9 = table2array(paras_r9_tv(:,3));

[g,xii] = ksdensity(gamma_r2,'Bandwidth',0.002);
[h,xiii] = ksdensity(gamma_r3,'Bandwidth',0.002);
[j,xiv] = ksdensity(gamma_r4,'Bandwidth',0.002);
[z,xv] = ksdensity(gamma_r5,'Bandwidth',0.002);
[m,xvi] = ksdensity(gamma_r6,'Bandwidth',0.002);

figure
yline(1/(0.26 - 0.005),'b','Linewidth',3.5);
xlim([0.005 0.05]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(gamma_init),0,'xb','markersize',20)
plot(mean(gamma_r2),0,'xg','markersize',20)
plot(mean(gamma_r3),0,'xr','markersize',20)
plot(mean(gamma_r4),0,'xc','markersize',20)
plot(mean(gamma_r5),0,'xm','markersize',20)
plot(mean(gamma_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\gamma$)','Post-round 1 density($\gamma$)','Post-round 2 density($\gamma$)','Post-round 3 density($\gamma$)','Post-round 4 density($\gamma$)', 'Post-round 5 density($\gamma$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$\gamma$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(gamma_r7,'Bandwidth',0.002);
[l,xviii] = ksdensity(gamma_r8,'Bandwidth',0.002);
[q,xix] = ksdensity(gamma_r9,'Bandwidth',0.002);

figure
yline(1/(0.26 - 0.005),'b','Linewidth',3.5);
xlim([0.005 0.05]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(gamma_init),0,'xb','markersize',20)
plot(mean(gamma_r7),0,'xg','markersize',20)
plot(mean(gamma_r8),0,'xr','markersize',20)
plot(mean(gamma_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\gamma$)','Post-round 6 density($\gamma$)','Post-round 7 density($\gamma$)','Post-round 8 density($\gamma$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$\gamma$ values')
ylabel('Probability density')

%% Density plot of rn
rn_init = table2array(paras_init_tv(:,4));
rn_r2 = table2array(paras_r2_tv(:,4));
rn_r3 = table2array(paras_r3_tv(:,4));
rn_r4 = table2array(paras_r4_tv(:,4));
rn_r5 = table2array(paras_r5_tv(:,4));
rn_r6 = table2array(paras_r6_tv(:,4));
rn_r7 = table2array(paras_r7_tv(:,4));
rn_r8 = table2array(paras_r8_tv(:,4));
rn_r9 = table2array(paras_r9_tv(:,4));

[g,xii] = ksdensity(rn_r2,'Bandwidth',0.008);
[h,xiii] = ksdensity(rn_r3,'Bandwidth',0.008);
[j,xiv] = ksdensity(rn_r4,'Bandwidth',0.008);
[z,xv] = ksdensity(rn_r5,'Bandwidth',0.008);
[m,xvi] = ksdensity(rn_r6,'Bandwidth',0.008);

figure
yline(1/(0.08 - 0.0008),'b','Linewidth',3.5);
xlim([0.0008 0.08]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(rn_init),0,'xb','markersize',20)
plot(mean(rn_r2),0,'xg','markersize',20)
plot(mean(rn_r3),0,'xr','markersize',20)
plot(mean(rn_r4),0,'xc','markersize',20)
plot(mean(rn_r5),0,'xm','markersize',20)
plot(mean(rn_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n}$)','Post-round 1 density($r_{n}$)','Post-round 2 density($r_{n}$)','Post-round 3 density($r_{n}$)','Post-round 4 density($r_{n}$)', 'Post-round 5 density($r_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$r_{n}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(rn_r7,'Bandwidth',0.008);
[l,xviii] = ksdensity(rn_r8,'Bandwidth',0.008);
[q,xix] = ksdensity(rn_r9,'Bandwidth',0.008);

figure
yline(1/(0.08 - 0.0008),'b','Linewidth',3.5);
xlim([0.0008 0.08]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(rn_init),0,'xb','markersize',20)
plot(mean(rn_r7),0,'xg','markersize',20)
plot(mean(rn_r8),0,'xr','markersize',20)
plot(mean(rn_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n}$)','Post-round 6 density($r_{n}$)','Post-round 7 density($r_{n}$)','Post-round 8 density($r_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$r_{n}$ values')
ylabel('Probability density')

%% Density plots of eta
eta_init = table2array(paras_init_tv(:,5));
eta_r2 = table2array(paras_r2_tv(:,5));
eta_r3 = table2array(paras_r3_tv(:,5));
eta_r4 = table2array(paras_r4_tv(:,5));
eta_r5 = table2array(paras_r5_tv(:,5));
eta_r6 = table2array(paras_r6_tv(:,5));
eta_r7 = table2array(paras_r7_tv(:,5));
eta_r8 = table2array(paras_r8_tv(:,5));
eta_r9 = table2array(paras_r9_tv(:,5));

[g,xii] = ksdensity(eta_r2,'Bandwidth',0.6);
[h,xiii] = ksdensity(eta_r3,'Bandwidth',0.6);
[j,xiv] = ksdensity(eta_r4,'Bandwidth',0.6);
[z,xv] = ksdensity(eta_r5,'Bandwidth',0.6);
[m,xvi] = ksdensity(eta_r6,'Bandwidth',0.6);

figure
yline(1/(18 - 7),'b','Linewidth',3.5);
xlim([7 18]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(eta_init),0,'xb','markersize',20)
plot(mean(eta_r2),0,'xg','markersize',20)
plot(mean(eta_r3),0,'xr','markersize',20)
plot(mean(eta_r4),0,'xc','markersize',20)
plot(mean(eta_r5),0,'xm','markersize',20)
plot(mean(eta_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\eta$)','Post-round 1 density($\eta$)','Post-round 2 density($\eta$)','Post-round 3 density($\eta$)','Post-round 4 density($\eta$)', 'Post-round 5 density($\eta$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\eta$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(eta_r7,'Bandwidth',0.6);
[l,xviii] = ksdensity(eta_r8,'Bandwidth',0.6);
[q,xix] = ksdensity(eta_r9,'Bandwidth',0.6);

figure
yline(1/(18 - 7),'b','Linewidth',3.5);
xlim([7 18]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(eta_init),0,'xb','markersize',20)
plot(mean(eta_r7),0,'xg','markersize',20)
plot(mean(eta_r8),0,'xr','markersize',20)
plot(mean(eta_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\eta$)','Post-round 6 density($\eta$)','Post-round 7 density($\eta$)','Post-round 8 density($\eta$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\eta$ values')
ylabel('Probability density')

%% Density plots of dm
dm_init = table2array(paras_init_tv(:,6));
dm_r2 = table2array(paras_r2_tv(:,6));
dm_r3 = table2array(paras_r3_tv(:,6));
dm_r4 = table2array(paras_r4_tv(:,6));
dm_r5 = table2array(paras_r5_tv(:,6));
dm_r6 = table2array(paras_r6_tv(:,6));
dm_r7 = table2array(paras_r7_tv(:,6));
dm_r8 = table2array(paras_r8_tv(:,6));
dm_r9 = table2array(paras_r9_tv(:,6));

[g,xii] = ksdensity(dm_r2,'Bandwidth',0.0025);
[h,xiii] = ksdensity(dm_r3,'Bandwidth',0.0025);
[j,xiv] = ksdensity(dm_r4,'Bandwidth',0.0025);
[z,xv] = ksdensity(dm_r5,'Bandwidth',0.0025);
[m,xvi] = ksdensity(dm_r6,'Bandwidth',0.0025);

figure
yline(1/(0.033 - 0.0001),'b','Linewidth',3.5);
xlim([0.0001 0.033]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(dm_init),0,'xb','markersize',20)
plot(mean(dm_r2),0,'xg','markersize',20)
plot(mean(dm_r3),0,'xr','markersize',20)
plot(mean(dm_r4),0,'xc','markersize',20)
plot(mean(dm_r5),0,'xm','markersize',20)
plot(mean(dm_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($d_{m}$)','Post-round 1 density($d_{m}$)','Post-round 2 density($d_{m}$)','Post-round 3 density($d_{m}$)','Post-round 4 density($d_{m}$)', 'Post-round 5 density($d_{m}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',7);
xlabel('$d_{m}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(dm_r7,'Bandwidth',0.0025);
[l,xviii] = ksdensity(dm_r8,'Bandwidth',0.0025);
[q,xix] = ksdensity(dm_r9,'Bandwidth',0.0025);

figure
yline(1/(0.033 - 0.0001),'b','Linewidth',3.5);
xlim([0.0001 0.033]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(dm_init),0,'xb','markersize',20)
plot(mean(dm_r7),0,'xg','markersize',20)
plot(mean(dm_r8),0,'xr','markersize',20)
plot(mean(dm_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($d_{m}$)','Post-round 6 density($d_{m}$)','Post-round 7 density($d_{m}$)','Post-round 8 density($d_{m}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$d_{m}$ values')
ylabel('Probability density')

%% Density plots of alpha
alpha_init = table2array(paras_init_tv(:,7));
alpha_r2 = table2array(paras_r2_tv(:,7));
alpha_r3 = table2array(paras_r3_tv(:,7));
alpha_r4 = table2array(paras_r4_tv(:,7));
alpha_r5 = table2array(paras_r5_tv(:,7));
alpha_r6 = table2array(paras_r6_tv(:,7));
alpha_r7 = table2array(paras_r7_tv(:,7));
alpha_r8 = table2array(paras_r8_tv(:,7));
alpha_r9 = table2array(paras_r9_tv(:,7));

[g,xii] = ksdensity(alpha_r2,'Bandwidth',0.008);
[h,xiii] = ksdensity(alpha_r3,'Bandwidth',0.008);
[j,xiv] = ksdensity(alpha_r4,'Bandwidth',0.008);
[z,xv] = ksdensity(alpha_r5,'Bandwidth',0.008);
[m,xvi] = ksdensity(alpha_r6,'Bandwidth',0.008);

figure
yline(1/(0.18 - 0.07),'b','Linewidth',3.5);
xlim([0.07 0.18]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(alpha_init),0,'xb','markersize',20)
plot(mean(alpha_r2),0,'xg','markersize',20)
plot(mean(alpha_r3),0,'xr','markersize',20)
plot(mean(alpha_r4),0,'xc','markersize',20)
plot(mean(alpha_r5),0,'xm','markersize',20)
plot(mean(alpha_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\alpha$)','Post-round 1 density($\alpha$)','Post-round 2 density($\alpha$)','Post-round 3 density($\alpha$)','Post-round 4 density($\alpha$)', 'Post-round 5 density($\alpha$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\alpha$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(alpha_r7,'Bandwidth',0.008);
[l,xviii] = ksdensity(alpha_r8,'Bandwidth',0.008);
[q,xix] = ksdensity(alpha_r9,'Bandwidth',0.008);

figure
yline(1/(0.18 - 0.07),'b','Linewidth',3.5);
xlim([0.07 0.18]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(alpha_init),0,'xb','markersize',20)
plot(mean(alpha_r7),0,'xg','markersize',20)
plot(mean(alpha_r8),0,'xr','markersize',20)
plot(mean(alpha_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\alpha$)','Post-round 6 density($\alpha$)','Post-round 7 density($\alpha$)','Post-round 8 density($\alpha$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\alpha$ values')
ylabel('Probability density')

%% Density plots of r.init
r_init_init = table2array(paras_init_tv(:,8));
r_init_r2 = table2array(paras_r2_tv(:,8));
r_init_r3 = table2array(paras_r3_tv(:,8));
r_init_r4 = table2array(paras_r4_tv(:,8));
r_init_r5 = table2array(paras_r5_tv(:,8));
r_init_r6 = table2array(paras_r6_tv(:,8));
r_init_r7 = table2array(paras_r7_tv(:,8));
r_init_r8 = table2array(paras_r8_tv(:,8));
r_init_r9 = table2array(paras_r9_tv(:,8));

[g,xii] = ksdensity(r_init_r2,'Bandwidth',0.3);
[h,xiii] = ksdensity(r_init_r3,'Bandwidth',0.3);
[j,xiv] = ksdensity(r_init_r4,'Bandwidth',0.3);
[z,xv] = ksdensity(r_init_r5,'Bandwidth',0.3);
[m,xvi] = ksdensity(r_init_r6,'Bandwidth',0.3);

figure
yline(1/(5 - 1),'b','Linewidth',3.5);
xlim([1 5]);
%ylim([0 13]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(r_init_init),0,'xb','markersize',20)
plot(mean(r_init_r2),0,'xg','markersize',20)
plot(mean(r_init_r3),0,'xr','markersize',20)
plot(mean(r_init_r4),0,'xc','markersize',20)
plot(mean(r_init_r5),0,'xm','markersize',20)
plot(mean(r_init_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($R_{init.}$)','Post-round 1 density($R_{init.}$)','Post-round 2 density($R_{init.}$)','Post-round 3 density($R_{init.}$)','Post-round 4 density($R_{init.}$)', 'Post-round 5 density($R_{init.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$R_{init.}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(r_init_r7,'Bandwidth',0.3);
[l,xviii] = ksdensity(r_init_r8,'Bandwidth',0.3);
[q,xix] = ksdensity(r_init_r9,'Bandwidth',0.3);

figure
yline(1/(5 - 1),'b','Linewidth',3.5);
xlim([1 5]);
%ylim([0 13]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(r_init_init),0,'xb','markersize',20)
plot(mean(r_init_r7),0,'xg','markersize',20)
plot(mean(r_init_r8),0,'xr','markersize',20)
plot(mean(r_init_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($R_{init.}$)','Post-round 6 density($R_{init.}$)','Post-round 7 density($R_{init.}$)','Post-round 8 density($R_{init.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$R_{init.}$ values')
ylabel('Probability density')

%% Density plots of p.ext
p_ext_init = table2array(paras_init_tv(:,9));
p_ext_r2 = table2array(paras_r2_tv(:,9));
p_ext_r3 = table2array(paras_r3_tv(:,9));
p_ext_r4 = table2array(paras_r4_tv(:,9));
p_ext_r5 = table2array(paras_r5_tv(:,9));
p_ext_r6 = table2array(paras_r6_tv(:,9));
p_ext_r7 = table2array(paras_r7_tv(:,9));
p_ext_r8 = table2array(paras_r8_tv(:,9));
p_ext_r9 = table2array(paras_r9_tv(:,9));

[g,xii] = ksdensity(p_ext_r2,'Bandwidth',0.008);
[h,xiii] = ksdensity(p_ext_r3,'Bandwidth',0.008);
[j,xiv] = ksdensity(p_ext_r4,'Bandwidth',0.008);
[z,xv] = ksdensity(p_ext_r5,'Bandwidth',0.008);
[m,xvi] = ksdensity(p_ext_r6,'Bandwidth',0.008);

figure
yline(1/(0.1 - 0.01),'b','Linewidth',3.5);
xlim([0.01 0.1]);
%ylim([0 13]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(p_ext_init),0,'xb','markersize',20)
plot(mean(p_ext_r2),0,'xg','markersize',20)
plot(mean(p_ext_r3),0,'xr','markersize',20)
plot(mean(p_ext_r4),0,'xc','markersize',20)
plot(mean(p_ext_r5),0,'xm','markersize',20)
plot(mean(p_ext_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($P_{ext.}$)','Post-round 1 density($P_{ext.}$)','Post-round 2 density($P_{ext.}$)','Post-round 3 density($P_{ext.}$)','Post-round 4 density($P_{ext.}$)', 'Post-round 5 density($P_{ext.}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(p_ext_r7,'Bandwidth',0.008);
[l,xviii] = ksdensity(p_ext_r8,'Bandwidth',0.008);
[q,xix] = ksdensity(p_ext_r9,'Bandwidth',0.008);

figure
yline(1/(0.1 - 0.01),'b','Linewidth',3.5);
xlim([0.01 0.1]);
%ylim([0 13]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(p_ext_init),0,'xb','markersize',20)
plot(mean(p_ext_r7),0,'xg','markersize',20)
plot(mean(p_ext_r8),0,'xr','markersize',20)
plot(mean(p_ext_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($P_{ext.}$)','Post-round 6 density($P_{ext.}$)','Post-round 7 density($P_{ext.}$)','Post-round 8 density($P_{ext.}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

%% Density plots of p.mit
p_mit_init = table2array(paras_init_tv(:,10));
p_mit_r2 = table2array(paras_r2_tv(:,10));
p_mit_r3 = table2array(paras_r3_tv(:,10));
p_mit_r4 = table2array(paras_r4_tv(:,10));
p_mit_r5 = table2array(paras_r5_tv(:,10));
p_mit_r6 = table2array(paras_r6_tv(:,10));
p_mit_r7 = table2array(paras_r7_tv(:,10));
p_mit_r8 = table2array(paras_r8_tv(:,10));
p_mit_r9 = table2array(paras_r9_tv(:,10));

[g,xii] = ksdensity(p_mit_r2,'Bandwidth',0.04);
[h,xiii] = ksdensity(p_mit_r3,'Bandwidth',0.04);
[j,xiv] = ksdensity(p_mit_r4,'Bandwidth',0.04);
[z,xv] = ksdensity(p_mit_r5,'Bandwidth',0.04);
[m,xvi] = ksdensity(p_mit_r6,'Bandwidth',0.04);

figure
yline(1/(1 - 0.2),'b','Linewidth',3.5);
xlim([0.2 1]);
%ylim([0 13]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(p_mit_init),0,'xb','markersize',20)
plot(mean(p_mit_r2),0,'xg','markersize',20)
plot(mean(p_mit_r3),0,'xr','markersize',20)
plot(mean(p_mit_r4),0,'xc','markersize',20)
plot(mean(p_mit_r5),0,'xm','markersize',20)
plot(mean(p_mit_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit.}$)','Post-round 1 density($P_{mit.}$)','Post-round 2 density($P_{mit.}$)','Post-round 3 density($P_{mit.}$)','Post-round 4 density($P_{mit.}$)', 'Post-round 5 density($P_{mit.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit.}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(p_mit_r7,'Bandwidth',0.04);
[l,xviii] = ksdensity(p_mit_r8,'Bandwidth',0.04);
[q,xix] = ksdensity(p_mit_r9,'Bandwidth',0.04);

figure
yline(1/(1 - 0.2),'b','Linewidth',3.5);
xlim([0.2 1]);
%ylim([0 13]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(p_mit_init),0,'xb','markersize',20)
plot(mean(p_mit_r7),0,'xg','markersize',20)
plot(mean(p_mit_r8),0,'xr','markersize',20)
plot(mean(p_mit_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit.}$)','Post-round 6 density($P_{mit.}$)','Post-round 7 density($P_{mit.}$)','Post-round 8 density($P_{mit.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit.}$ values')
ylabel('Probability density')

%% Density plot of dn_quad
dn_quad_init = table2array(paras_init_tv(:,11));
dn_quad_r2 = table2array(paras_r2_tv(:,11));
dn_quad_r3 = table2array(paras_r3_tv(:,11));
dn_quad_r4 = table2array(paras_r4_tv(:,11));
dn_quad_r5 = table2array(paras_r5_tv(:,11));
dn_quad_r6 = table2array(paras_r6_tv(:,11));
dn_quad_r7 = table2array(paras_r7_tv(:,11));
dn_quad_r8 = table2array(paras_r8_tv(:,11));
dn_quad_r9 = table2array(paras_r9_tv(:,11));

dn_quad_lb = table2array(post_reg_coefs_bounds(1,2));
dn_quad_ub = table2array(post_reg_coefs_bounds(1,3));

[g,xii] = ksdensity(dn_quad_r2,'Bandwidth',0.00005);
[h,xiii] = ksdensity(dn_quad_r3,'Bandwidth',0.00005);
[j,xiv] = ksdensity(dn_quad_r4,'Bandwidth',0.00005);
[z,xv] = ksdensity(dn_quad_r5,'Bandwidth',0.00005);
[m,xvi] = ksdensity(dn_quad_r6,'Bandwidth',0.00005);

figure
yline(1/(dn_quad_ub-dn_quad_lb),'b','Linewidth',3.5);
xlim([dn_quad_lb dn_quad_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(dn_quad_init),0,'xb','markersize',20)
plot(mean(dn_quad_r2),0,'xg','markersize',20)
plot(mean(dn_quad_r3),0,'xr','markersize',20)
plot(mean(dn_quad_r4),0,'xc','markersize',20)
plot(mean(dn_quad_r5),0,'xm','markersize',20)
plot(mean(dn_quad_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n, quad}$)','Post-round 1 density($d_{n, quad}$)','Post-round 2 density($d_{n, quad}$)','Post-round 3 density($d_{n, quad}$)','Post-round 4 density($d_{n, quad}$)', 'Post-round 5 density($d_{n, quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',8);
xlabel('$d_{n, quad}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(dn_quad_r7,'Bandwidth',0.00005);
[l,xviii] = ksdensity(dn_quad_r8,'Bandwidth',0.00005);
[q,xix] = ksdensity(dn_quad_r9,'Bandwidth',0.00005);

figure
yline(1/(dn_quad_ub-dn_quad_lb),'b','Linewidth',3.5);
xlim([dn_quad_lb dn_quad_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(dn_quad_r7),0,'xg','markersize',20)
plot(mean(dn_quad_r8),0,'xr','markersize',20)
plot(mean(dn_quad_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n, quad}$)','Post-round 6 density($d_{n, quad}$)','Post-round 7 density($d_{n, quad}$)','Post-round 8 density($d_{n, quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',8);
xlabel('$d_{n, quad}$ values')
ylabel('Probability density')

%% Density plot of dn_lin
dn_lin_init = table2array(paras_init_tv(:,12));
dn_lin_r2 = table2array(paras_r2_tv(:,12));
dn_lin_r3 = table2array(paras_r3_tv(:,12));
dn_lin_r4 = table2array(paras_r4_tv(:,12));
dn_lin_r5 = table2array(paras_r5_tv(:,12));
dn_lin_r6 = table2array(paras_r6_tv(:,12));
dn_lin_r7 = table2array(paras_r7_tv(:,12));
dn_lin_r8 = table2array(paras_r8_tv(:,12));
dn_lin_r9 = table2array(paras_r9_tv(:,12));

dn_lin_lb = table2array(post_reg_coefs_bounds(2,2));
dn_lin_ub = table2array(post_reg_coefs_bounds(2,3));

[g,xii] = ksdensity(dn_lin_r2,'Bandwidth',0.00015);
[h,xiii] = ksdensity(dn_lin_r3,'Bandwidth',0.00015);
[j,xiv] = ksdensity(dn_lin_r4,'Bandwidth',0.00015);
[z,xv] = ksdensity(dn_lin_r5,'Bandwidth',0.00015);
[m,xvi] = ksdensity(dn_lin_r6,'Bandwidth',0.00015);

figure
yline(1/(dn_lin_ub-dn_lin_lb),'b','Linewidth',3.5);
xlim([dn_lin_lb dn_lin_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(dn_lin_init),0,'xb','markersize',20)
plot(mean(dn_lin_r2),0,'xg','markersize',20)
plot(mean(dn_lin_r3),0,'xr','markersize',20)
plot(mean(dn_lin_r4),0,'xc','markersize',20)
plot(mean(dn_lin_r5),0,'xm','markersize',20)
plot(mean(dn_lin_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n, lin}$)','Post-round 1 density($d_{n, lin}$)','Post-round 2 density($d_{n, lin}$)','Post-round 3 density($d_{n, lin}$)','Post-round 4 density($d_{n, lin}$)', 'Post-round 5 density($d_{n, lin}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',8);
xlabel('$d_{n, lin}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(dn_lin_r7,'Bandwidth',0.00015);
[l,xviii] = ksdensity(dn_lin_r8,'Bandwidth',0.00015);
[q,xix] = ksdensity(dn_lin_r9,'Bandwidth',0.00015);

figure
yline(1/(dn_lin_ub-dn_lin_lb),'b','Linewidth',3.5);
xlim([dn_lin_lb dn_lin_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(dn_lin_r7),0,'xg','markersize',20)
plot(mean(dn_lin_r8),0,'xr','markersize',20)
plot(mean(dn_lin_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n, lin}$)','Post-round 6 density($d_{n, lin}$)','Post-round 7 density($d_{n, lin}$)','Post-round 8 density($d_{n, lin}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',8);
xlabel('$d_{n, lin}$ values')
ylabel('Probability density')

%% Density plot of gamma_quad
gamma_quad_init = table2array(paras_init_tv(:,13));
gamma_quad_r2 = table2array(paras_r2_tv(:,13));
gamma_quad_r3 = table2array(paras_r3_tv(:,13));
gamma_quad_r4 = table2array(paras_r4_tv(:,13));
gamma_quad_r5 = table2array(paras_r5_tv(:,13));
gamma_quad_r6 = table2array(paras_r6_tv(:,13));
gamma_quad_r7 = table2array(paras_r7_tv(:,13));
gamma_quad_r8 = table2array(paras_r8_tv(:,13));
gamma_quad_r9 = table2array(paras_r9_tv(:,13));

gamma_quad_lb = table2array(post_reg_coefs_bounds(3,2));
gamma_quad_ub = table2array(post_reg_coefs_bounds(3,3));

[g,xii] = ksdensity(gamma_quad_r2,'Bandwidth',0.001);
[h,xiii] = ksdensity(gamma_quad_r3,'Bandwidth',0.001);
[j,xiv] = ksdensity(gamma_quad_r4,'Bandwidth',0.001);
[z,xv] = ksdensity(gamma_quad_r5,'Bandwidth',0.001);
[m,xvi] = ksdensity(gamma_quad_r6,'Bandwidth',0.001);

figure
yline(1/(gamma_quad_ub-gamma_quad_lb),'b','Linewidth',3.5);
xlim([gamma_quad_lb -0.012]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(gamma_quad_init),0,'xb','markersize',20)
plot(mean(gamma_quad_r2),0,'xg','markersize',20)
plot(mean(gamma_quad_r3),0,'xr','markersize',20)
plot(mean(gamma_quad_r4),0,'xc','markersize',20)
plot(mean(gamma_quad_r5),0,'xm','markersize',20)
plot(mean(gamma_quad_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\gamma_{quad}$)','Post-round 1 density($\gamma_{quad}$)','Post-round 2 density($\gamma_{quad}$)','Post-round 3 density($\gamma_{quad}$)','Post-round 4 density($\gamma_{quad}$)', 'Post-round 5 density($\gamma_{quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\gamma_{quad}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(gamma_quad_r7,'Bandwidth',0.001);
[l,xviii] = ksdensity(gamma_quad_r8,'Bandwidth',0.001);
[q,xix] = ksdensity(gamma_quad_r9,'Bandwidth',0.001);

figure
yline(1/(gamma_quad_ub-gamma_quad_lb),'b','Linewidth',3.5);
xlim([gamma_quad_lb -0.012]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(gamma_quad_r7),0,'xg','markersize',20)
plot(mean(gamma_quad_r8),0,'xr','markersize',20)
plot(mean(gamma_quad_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\gamma_{quad}$)','Post-round 6 density($\gamma_{quad}$)','Post-round 7 density($\gamma_{quad}$)','Post-round 8 density($\gamma_{quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',8);
xlabel('$\gamma_{quad}$ values')
ylabel('Probability density')

%% Density plot of gamma_lin
gamma_lin_init = table2array(paras_init_tv(:,14));
gamma_lin_r2 = table2array(paras_r2_tv(:,14));
gamma_lin_r3 = table2array(paras_r3_tv(:,14));
gamma_lin_r4 = table2array(paras_r4_tv(:,14));
gamma_lin_r5 = table2array(paras_r5_tv(:,14));
gamma_lin_r6 = table2array(paras_r6_tv(:,14));
gamma_lin_r7 = table2array(paras_r7_tv(:,14));
gamma_lin_r8 = table2array(paras_r8_tv(:,14));
gamma_lin_r9 = table2array(paras_r9_tv(:,14));

gamma_lin_lb = table2array(post_reg_coefs_bounds(4,2));
gamma_lin_ub = table2array(post_reg_coefs_bounds(4,3));

[g,xii] = ksdensity(gamma_lin_r2,'Bandwidth',0.015);
[h,xiii] = ksdensity(gamma_lin_r3,'Bandwidth',0.015);
[j,xiv] = ksdensity(gamma_lin_r4,'Bandwidth',0.015);
[z,xv] = ksdensity(gamma_lin_r5,'Bandwidth',0.015);
[m,xvi] = ksdensity(gamma_lin_r6,'Bandwidth',0.015);

figure
yline(1/(gamma_lin_ub-gamma_lin_lb),'b','Linewidth',3.5);
xlim([gamma_lin_lb gamma_lin_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(gamma_lin_init),0,'xb','markersize',20)
plot(mean(gamma_lin_r2),0,'xg','markersize',20)
plot(mean(gamma_lin_r3),0,'xr','markersize',20)
plot(mean(gamma_lin_r4),0,'xc','markersize',20)
plot(mean(gamma_lin_r5),0,'xm','markersize',20)
plot(mean(gamma_lin_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\gamma_{lin}$)','Post-round 1 density($\gamma_{lin}$)','Post-round 2 density($\gamma_{lin}$)','Post-round 3 density($\gamma_{lin}$)','Post-round 4 density($\gamma_{lin}$)', 'Post-round 5 density($\gamma_{lin}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$\gamma_{lin}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(gamma_lin_r7,'Bandwidth',0.015);
[l,xviii] = ksdensity(gamma_lin_r8,'Bandwidth',0.015);
[q,xix] = ksdensity(gamma_lin_r9,'Bandwidth',0.015);

figure
yline(1/(gamma_lin_ub-gamma_lin_lb),'b','Linewidth',3.5);
xlim([gamma_lin_lb gamma_lin_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(gamma_lin_r7),0,'xg','markersize',20)
plot(mean(gamma_lin_r8),0,'xr','markersize',20)
plot(mean(gamma_lin_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\gamma_{lin}$)','Post-round 6 density($\gamma_{lin}$)','Post-round 7 density($\gamma_{lin}$)','Post-round 8 density($\gamma_{lin}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$\gamma_{lin}$ values')
ylabel('Probability density')

%% Density plot of rn_quad
rn_quad_init = table2array(paras_init_tv(:,15));
rn_quad_r2 = table2array(paras_r2_tv(:,15));
rn_quad_r3 = table2array(paras_r3_tv(:,15));
rn_quad_r4 = table2array(paras_r4_tv(:,15));
rn_quad_r5 = table2array(paras_r5_tv(:,15));
rn_quad_r6 = table2array(paras_r6_tv(:,15));
rn_quad_r7 = table2array(paras_r7_tv(:,15));
rn_quad_r8 = table2array(paras_r8_tv(:,15));
rn_quad_r9 = table2array(paras_r9_tv(:,15));

rn_quad_lb = table2array(post_reg_coefs_bounds(5,2));
rn_quad_ub = table2array(post_reg_coefs_bounds(5,3));

[g,xii] = ksdensity(rn_quad_r2,'Bandwidth',0.0005);
[h,xiii] = ksdensity(rn_quad_r3,'Bandwidth',0.0005);
[j,xiv] = ksdensity(rn_quad_r4,'Bandwidth',0.0005);
[z,xv] = ksdensity(rn_quad_r5,'Bandwidth',0.0005);
[m,xvi] = ksdensity(rn_quad_r6,'Bandwidth',0.0005);

figure
yline(1/(rn_quad_ub-rn_quad_lb),'b','Linewidth',3.5);
xlim([rn_quad_lb rn_quad_ub]);
ylim([0 600]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(rn_quad_init),0,'xb','markersize',20)
plot(mean(rn_quad_r2),0,'xg','markersize',20)
plot(mean(rn_quad_r3),0,'xr','markersize',20)
plot(mean(rn_quad_r4),0,'xc','markersize',20)
plot(mean(rn_quad_r5),0,'xm','markersize',20)
plot(mean(rn_quad_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n, quad}$)','Post-round 1 density($r_{n, quad}$)','Post-round 2 density($r_{n, quad}$)','Post-round 3 density($r_{n, quad}$)','Post-round 4 density($r_{n, quad}$)', 'Post-round 5 density($r_{n, quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',8);
xlabel('$r_{n, quad}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(rn_quad_r7,'Bandwidth',0.0005);
[l,xviii] = ksdensity(rn_quad_r8,'Bandwidth',0.0005);
[q,xix] = ksdensity(rn_quad_r9,'Bandwidth',0.0005);

figure
yline(1/(rn_quad_ub-rn_quad_lb),'b','Linewidth',3.5);
xlim([rn_quad_lb rn_quad_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(rn_quad_r7),0,'xg','markersize',20)
plot(mean(rn_quad_r8),0,'xr','markersize',20)
plot(mean(rn_quad_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n, quad}$)','Post-round 6 density($r_{n, quad}$)','Post-round 7 density($r_{n, quad}$)','Post-round 8 density($r_{n, quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',7);
xlabel('$r_{n, quad}$ values')
ylabel('Probability density')

%% Density plot of rn_lin
rn_lin_init = table2array(paras_init_tv(:,16));
rn_lin_r2 = table2array(paras_r2_tv(:,16));
rn_lin_r3 = table2array(paras_r3_tv(:,16));
rn_lin_r4 = table2array(paras_r4_tv(:,16));
rn_lin_r5 = table2array(paras_r5_tv(:,16));
rn_lin_r6 = table2array(paras_r6_tv(:,16));
rn_lin_r7 = table2array(paras_r7_tv(:,16));
rn_lin_r8 = table2array(paras_r8_tv(:,16));
rn_lin_r9 = table2array(paras_r9_tv(:,16));

rn_lin_lb = table2array(post_reg_coefs_bounds(6,2));
rn_lin_ub = table2array(post_reg_coefs_bounds(6,3));

[g,xii] = ksdensity(rn_lin_r2,'Bandwidth',0.0025);
[h,xiii] = ksdensity(rn_lin_r3,'Bandwidth',0.0025);
[j,xiv] = ksdensity(rn_lin_r4,'Bandwidth',0.0025);
[z,xv] = ksdensity(rn_lin_r5,'Bandwidth',0.0025);
[m,xvi] = ksdensity(rn_lin_r6,'Bandwidth',0.0025);

figure
yline(1/(rn_lin_ub-rn_lin_lb),'b','Linewidth',3.5);
xlim([rn_lin_lb rn_lin_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(rn_lin_init),0,'xb','markersize',20)
plot(mean(rn_lin_r2),0,'xg','markersize',20)
plot(mean(rn_lin_r3),0,'xr','markersize',20)
plot(mean(rn_lin_r4),0,'xc','markersize',20)
plot(mean(rn_lin_r5),0,'xm','markersize',20)
plot(mean(rn_lin_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n, lin}$)','Post-round 1 density($r_{n, lin}$)','Post-round 2 density($r_{n, lin}$)','Post-round 3 density($r_{n, lin}$)','Post-round 4 density($r_{n, lin}$)', 'Post-round 5 density($r_{n, lin}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$r_{n, lin}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(rn_lin_r7,'Bandwidth',0.0025);
[l,xviii] = ksdensity(rn_lin_r8,'Bandwidth',0.0025);
[q,xix] = ksdensity(rn_lin_r9,'Bandwidth',0.0025);

figure
yline(1/(rn_lin_ub-rn_lin_lb),'b','Linewidth',3.5);
xlim([rn_lin_lb rn_lin_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(rn_lin_r7),0,'xg','markersize',20)
plot(mean(rn_lin_r8),0,'xr','markersize',20)
plot(mean(rn_lin_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n, lin}$)','Post-round 6 density($r_{n, lin}$)','Post-round 7 density($r_{n, lin}$)','Post-round 8 density($r_{n, lin}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$r_{n, lin}$ values')
ylabel('Probability density')

%% Density plot of eta_quad
eta_quad_init = table2array(paras_init_tv(:,17));
eta_quad_r2 = table2array(paras_r2_tv(:,17));
eta_quad_r3 = table2array(paras_r3_tv(:,17));
eta_quad_r4 = table2array(paras_r4_tv(:,17));
eta_quad_r5 = table2array(paras_r5_tv(:,17));
eta_quad_r6 = table2array(paras_r6_tv(:,17));
eta_quad_r7 = table2array(paras_r7_tv(:,17));
eta_quad_r8 = table2array(paras_r8_tv(:,17));
eta_quad_r9 = table2array(paras_r9_tv(:,17));

eta_quad_lb = table2array(post_reg_coefs_bounds(7,2));
eta_quad_ub = table2array(post_reg_coefs_bounds(7,3));

[g,xii] = ksdensity(eta_quad_r2,'Bandwidth',0.07);
[h,xiii] = ksdensity(eta_quad_r3,'Bandwidth',0.07);
[j,xiv] = ksdensity(eta_quad_r4,'Bandwidth',0.07);
[z,xv] = ksdensity(eta_quad_r5,'Bandwidth',0.07);
[m,xvi] = ksdensity(eta_quad_r6,'Bandwidth',0.07);

figure
yline(1/(eta_quad_ub-eta_quad_lb),'b','Linewidth',3.5);
xlim([eta_quad_lb eta_quad_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(eta_quad_init),0,'xb','markersize',20)
plot(mean(eta_quad_r2),0,'xg','markersize',20)
plot(mean(eta_quad_r3),0,'xr','markersize',20)
plot(mean(eta_quad_r4),0,'xc','markersize',20)
plot(mean(eta_quad_r5),0,'xm','markersize',20)
plot(mean(eta_quad_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\eta_{quad}$)','Post-round 1 density($\eta_{quad}$)','Post-round 2 density($\eta_{quad}$)','Post-round 3 density($\eta_{quad}$)','Post-round 4 density($\eta_{quad}$)', 'Post-round 5 density($\eta_{quad}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$\eta_{quad}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(eta_quad_r7,'Bandwidth',0.07);
[l,xviii] = ksdensity(eta_quad_r8,'Bandwidth',0.07);
[q,xix] = ksdensity(eta_quad_r9,'Bandwidth',0.07);

figure
yline(1/(eta_quad_ub-eta_quad_lb),'b','Linewidth',3.5);
xlim([eta_quad_lb eta_quad_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(eta_quad_r7),0,'xg','markersize',20)
plot(mean(eta_quad_r8),0,'xr','markersize',20)
plot(mean(eta_quad_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\eta_{quad}$)','Post-round 6 density($\eta_{quad}$)','Post-round 7 density($\eta_{quad}$)','Post-round 8 density($\eta_{quad}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$\eta_{quad}$ values')
ylabel('Probability density')

%% Density plot of eta_lin
eta_lin_init = table2array(paras_init_tv(:,18));
eta_lin_r2 = table2array(paras_r2_tv(:,18));
eta_lin_r3 = table2array(paras_r3_tv(:,18));
eta_lin_r4 = table2array(paras_r4_tv(:,18));
eta_lin_r5 = table2array(paras_r5_tv(:,18));
eta_lin_r6 = table2array(paras_r6_tv(:,18));
eta_lin_r7 = table2array(paras_r7_tv(:,18));
eta_lin_r8 = table2array(paras_r8_tv(:,18));
eta_lin_r9 = table2array(paras_r9_tv(:,18));

eta_lin_lb = table2array(post_reg_coefs_bounds(8,2));
eta_lin_ub = table2array(post_reg_coefs_bounds(8,3));

[g,xii] = ksdensity(eta_lin_r2,'Bandwidth',0.7);
[h,xiii] = ksdensity(eta_lin_r3,'Bandwidth',0.7);
[j,xiv] = ksdensity(eta_lin_r4,'Bandwidth',0.7);
[z,xv] = ksdensity(eta_lin_r5,'Bandwidth',0.7);
[m,xvi] = ksdensity(eta_lin_r6,'Bandwidth',0.7);

figure
yline(1/(eta_lin_ub-eta_lin_lb),'b','Linewidth',3.5);
xlim([eta_lin_lb eta_lin_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(eta_lin_init),0,'xb','markersize',20)
plot(mean(eta_lin_r2),0,'xg','markersize',20)
plot(mean(eta_lin_r3),0,'xr','markersize',20)
plot(mean(eta_lin_r4),0,'xc','markersize',20)
plot(mean(eta_lin_r5),0,'xm','markersize',20)
plot(mean(eta_lin_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\eta_{lin}$)','Post-round 1 density($\eta_{lin}$)','Post-round 2 density($\eta_{lin}$)','Post-round 3 density($\eta_{lin}$)','Post-round 4 density($\eta_{lin}$)', 'Post-round 5 density($\eta_{lin}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\eta_{lin}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(eta_lin_r7,'Bandwidth',0.7);
[l,xviii] = ksdensity(eta_lin_r8,'Bandwidth',0.7);
[q,xix] = ksdensity(eta_lin_r9,'Bandwidth',0.7);

figure
yline(1/(eta_lin_ub-eta_lin_lb),'b','Linewidth',3.5);
xlim([eta_lin_lb eta_lin_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(eta_lin_r7),0,'xg','markersize',20)
plot(mean(eta_lin_r8),0,'xr','markersize',20)
plot(mean(eta_lin_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\eta_{lin}$)','Post-round 6 density($\eta_{quad}$)','Post-round 7 density($\eta_{quad}$)','Post-round 8 density($\eta_{quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\eta_{lin}$ values')
ylabel('Probability density')

%% Density plot of alpha_quad
alpha_quad_init = table2array(paras_init_tv(:,19));
alpha_quad_r2 = table2array(paras_r2_tv(:,19));
alpha_quad_r3 = table2array(paras_r3_tv(:,19));
alpha_quad_r4 = table2array(paras_r4_tv(:,19));
alpha_quad_r5 = table2array(paras_r5_tv(:,19));
alpha_quad_r6 = table2array(paras_r6_tv(:,19));
alpha_quad_r7 = table2array(paras_r7_tv(:,19));
alpha_quad_r8 = table2array(paras_r8_tv(:,19));
alpha_quad_r9 = table2array(paras_r9_tv(:,19));

alpha_quad_lb = table2array(post_reg_coefs_bounds(9,2));
alpha_quad_ub = table2array(post_reg_coefs_bounds(9,3));

[g,xii] = ksdensity(alpha_quad_r2,'Bandwidth',0.0008);
[h,xiii] = ksdensity(alpha_quad_r3,'Bandwidth',0.0008);
[j,xiv] = ksdensity(alpha_quad_r4,'Bandwidth',0.0008);
[z,xv] = ksdensity(alpha_quad_r5,'Bandwidth',0.0008);
[m,xvi] = ksdensity(alpha_quad_r6,'Bandwidth',0.0008);

figure
yline(1/(alpha_quad_ub-alpha_quad_lb),'b','Linewidth',3.5);
xlim([alpha_quad_lb alpha_quad_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(alpha_quad_init),0,'xb','markersize',20)
plot(mean(alpha_quad_r2),0,'xg','markersize',20)
plot(mean(alpha_quad_r3),0,'xr','markersize',20)
plot(mean(alpha_quad_r4),0,'xc','markersize',20)
plot(mean(alpha_quad_r5),0,'xm','markersize',20)
plot(mean(alpha_quad_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\alpha_{quad}$)','Post-round 1 density($\alpha_{quad}$)','Post-round 2 density($\alpha_{quad}$)','Post-round 3 density($\alpha_{quad}$)','Post-round 4 density($\alpha_{quad}$)', 'Post-round 5 density($\alpha_{quad}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$\alpha_{quad}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(alpha_quad_r7,'Bandwidth',0.0008);
[l,xviii] = ksdensity(alpha_quad_r8,'Bandwidth',0.0008);
[q,xix] = ksdensity(alpha_quad_r9,'Bandwidth',0.0008);

figure
yline(1/(alpha_quad_ub-alpha_quad_lb),'b','Linewidth',3.5);
xlim([alpha_quad_lb alpha_quad_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(alpha_quad_r7),0,'xg','markersize',20)
plot(mean(alpha_quad_r8),0,'xr','markersize',20)
plot(mean(alpha_quad_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\alpha_{quad}$)','Post-round 6 density($\alpha_{quad}$)','Post-round 7 density($\alpha_{quad}$)','Post-round 8 density($\alpha_{quad}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$\alpha_{quad}$ values')
ylabel('Probability density')

%% Density plot of alpha_lin
alpha_lin_init = table2array(paras_init_tv(:,20));
alpha_lin_r2 = table2array(paras_r2_tv(:,20));
alpha_lin_r3 = table2array(paras_r3_tv(:,20));
alpha_lin_r4 = table2array(paras_r4_tv(:,20));
alpha_lin_r5 = table2array(paras_r5_tv(:,20));
alpha_lin_r6 = table2array(paras_r6_tv(:,20));
alpha_lin_r7 = table2array(paras_r7_tv(:,20));
alpha_lin_r8 = table2array(paras_r8_tv(:,20));
alpha_lin_r9 = table2array(paras_r9_tv(:,20));

alpha_lin_lb = table2array(post_reg_coefs_bounds(10,2));
alpha_lin_ub = table2array(post_reg_coefs_bounds(10,3));

[g,xii] = ksdensity(alpha_lin_r2,'Bandwidth',0.005);
[h,xiii] = ksdensity(alpha_lin_r3,'Bandwidth',0.005);
[j,xiv] = ksdensity(alpha_lin_r4,'Bandwidth',0.005);
[z,xv] = ksdensity(alpha_lin_r5,'Bandwidth',0.005);
[m,xvi] = ksdensity(alpha_lin_r6,'Bandwidth',0.005);

figure
yline(1/(alpha_lin_ub-alpha_lin_lb),'b','Linewidth',3.5);
xlim([alpha_lin_lb alpha_lin_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(alpha_lin_init),0,'xb','markersize',20)
plot(mean(alpha_lin_r2),0,'xg','markersize',20)
plot(mean(alpha_lin_r3),0,'xr','markersize',20)
plot(mean(alpha_lin_r4),0,'xc','markersize',20)
plot(mean(alpha_lin_r5),0,'xm','markersize',20)
plot(mean(alpha_lin_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\alpha_{lin}$)','Post-round 1 density($\alpha_{lin}$)','Post-round 2 density($\alpha_{lin}$)','Post-round 3 density($\alpha_{lin}$)','Post-round 4 density($\alpha_{lin}$)', 'Post-round 5 density($\alpha_{lin}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\alpha_{lin}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(alpha_lin_r7,'Bandwidth',0.005);
[l,xviii] = ksdensity(alpha_lin_r8,'Bandwidth',0.005);
[q,xix] = ksdensity(alpha_lin_r9,'Bandwidth',0.005);

figure
yline(1/(alpha_lin_ub-alpha_lin_lb),'b','Linewidth',3.5);
xlim([alpha_lin_lb alpha_lin_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(alpha_lin_r7),0,'xg','markersize',20)
plot(mean(alpha_lin_r8),0,'xr','markersize',20)
plot(mean(alpha_lin_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($\alpha_{lin}$)','Post-round 6 density($\alpha_{lin}$)','Post-round 7 density($\alpha_{lin}$)','Post-round 8 density($\alpha_{lin}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\alpha_{lin}$ values')
ylabel('Probability density')

%% Density plot of p_mit_quad
p_mit_quad_init = table2array(paras_init_tv(:,21));
p_mit_quad_r2 = table2array(paras_r2_tv(:,21));
p_mit_quad_r3 = table2array(paras_r3_tv(:,21));
p_mit_quad_r4 = table2array(paras_r4_tv(:,21));
p_mit_quad_r5 = table2array(paras_r5_tv(:,21));
p_mit_quad_r6 = table2array(paras_r6_tv(:,21));
p_mit_quad_r7 = table2array(paras_r7_tv(:,21));
p_mit_quad_r8 = table2array(paras_r8_tv(:,21));
p_mit_quad_r9 = table2array(paras_r9_tv(:,21));

p_mit_quad_lb = table2array(post_reg_coefs_bounds(11,2));
p_mit_quad_ub = table2array(post_reg_coefs_bounds(11,3));

[g,xii] = ksdensity(p_mit_quad_r2,'Bandwidth',0.008);
[h,xiii] = ksdensity(p_mit_quad_r3,'Bandwidth',0.008);
[j,xiv] = ksdensity(p_mit_quad_r4,'Bandwidth',0.008);
[z,xv] = ksdensity(p_mit_quad_r5,'Bandwidth',0.008);
[m,xvi] = ksdensity(p_mit_quad_r6,'Bandwidth',0.008);

figure
yline(1/(p_mit_quad_ub-p_mit_quad_lb),'b','Linewidth',3.5);
xlim([-0.2 p_mit_quad_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(p_mit_quad_init),0,'xb','markersize',20)
plot(mean(p_mit_quad_r2),0,'xg','markersize',20)
plot(mean(p_mit_quad_r3),0,'xr','markersize',20)
plot(mean(p_mit_quad_r4),0,'xc','markersize',20)
plot(mean(p_mit_quad_r5),0,'xm','markersize',20)
plot(mean(p_mit_quad_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit, quad}$)','Post-round 1 density($P_{mit, quad}$)','Post-round 2 density($P_{mit, quad}$)','Post-round 3 density($P_{mit, quad}$)','Post-round 4 density($P_{mit, quad}$)', 'Post-round 5 density($P_{mit, quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit, quad}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(p_mit_quad_r7,'Bandwidth',0.008);
[l,xviii] = ksdensity(p_mit_quad_r8,'Bandwidth',0.008);
[q,xix] = ksdensity(p_mit_quad_r9,'Bandwidth',0.008);

figure
yline(1/(p_mit_quad_ub-p_mit_quad_lb),'b','Linewidth',3.5);
xlim([-0.2 p_mit_quad_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(p_mit_quad_r7),0,'xg','markersize',20)
plot(mean(p_mit_quad_r8),0,'xr','markersize',20)
plot(mean(p_mit_quad_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit, quad}$)','Post-round 6 density($P_{mit, quad}$)','Post-round 7 density($P_{mit, quad}$)','Post-round 8 density($P_{mit, quad}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit, quad}$ values')
ylabel('Probability density')

%% Density plot of p_mit_lin
p_mit_lin_init = table2array(paras_init_tv(:,22));
p_mit_lin_r2 = table2array(paras_r2_tv(:,22));
p_mit_lin_r3 = table2array(paras_r3_tv(:,22));
p_mit_lin_r4 = table2array(paras_r4_tv(:,22));
p_mit_lin_r5 = table2array(paras_r5_tv(:,22));
p_mit_lin_r6 = table2array(paras_r6_tv(:,22));
p_mit_lin_r7 = table2array(paras_r7_tv(:,22));
p_mit_lin_r8 = table2array(paras_r8_tv(:,22));
p_mit_lin_r9 = table2array(paras_r9_tv(:,22));

p_mit_lin_lb = table2array(post_reg_coefs_bounds(12,2));
p_mit_lin_ub = table2array(post_reg_coefs_bounds(12,3));

[g,xii] = ksdensity(p_mit_lin_r2,'Bandwidth',0.05);
[h,xiii] = ksdensity(p_mit_lin_r3,'Bandwidth',0.05);
[j,xiv] = ksdensity(p_mit_lin_r4,'Bandwidth',0.05);
[z,xv] = ksdensity(p_mit_lin_r5,'Bandwidth',0.05);
[m,xvi] = ksdensity(p_mit_lin_r6,'Bandwidth',0.05);

figure
yline(1/(p_mit_lin_ub-p_mit_lin_lb),'b','Linewidth',3.5);
xlim([p_mit_lin_lb p_mit_lin_ub]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(p_mit_lin_init),0,'xb','markersize',20)
plot(mean(p_mit_lin_r2),0,'xg','markersize',20)
plot(mean(p_mit_lin_r3),0,'xr','markersize',20)
plot(mean(p_mit_lin_r4),0,'xc','markersize',20)
plot(mean(p_mit_lin_r5),0,'xm','markersize',20)
plot(mean(p_mit_lin_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit, lin}$)','Post-round 1 density($P_{mit, lin}$)','Post-round 2 density($P_{mit, lin}$)','Post-round 3 density($P_{mit, lin}$)','Post-round 4 density($P_{mit, lin}$)', 'Post-round 5 density($P_{mit, lin}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit, lin}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(p_mit_lin_r7,'Bandwidth',0.05);
[l,xviii] = ksdensity(p_mit_lin_r8,'Bandwidth',0.05);
[q,xix] = ksdensity(p_mit_lin_r9,'Bandwidth',0.05);

figure
yline(1/(p_mit_lin_ub-p_mit_lin_lb),'b','Linewidth',3.5);
xlim([p_mit_lin_lb p_mit_lin_ub]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-')
plot(mean(p_mit_lin_r7),0,'xg','markersize',20)
plot(mean(p_mit_lin_r8),0,'xr','markersize',20)
plot(mean(p_mit_lin_r9),0,'xc','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit, lin}$)','Post-round 6 density($P_{mit, lin}$)','Post-round 7 density($P_{mit, lin}$)','Post-round 8 density($P_{mit, lin}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit, lin}$ values')
ylabel('Probability density')