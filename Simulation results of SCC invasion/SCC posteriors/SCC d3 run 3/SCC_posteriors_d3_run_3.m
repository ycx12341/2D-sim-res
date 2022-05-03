% SCC_posteriors_d3_run_3.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of probability densities of the
% parameter estimations at the end of each round of the second run of 
% applying ABC scheme on the post-day 3 pattern of the SCC reference 
% dataset.

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
paras_init_d3 = readtable("Round 1 initial parameters.txt");
paras_r2_d3 = readtable("Round 2 parameters.txt");
paras_r3_d3 = readtable("Round 3 parameters.txt");
paras_r4_d3 = readtable("Round 4 parameters.txt");
paras_r5_d3 = readtable("Round 5 parameters.txt");
paras_r6_d3 = readtable("Round 6 parameters.txt");
paras_r7_d3 = readtable("Round 7 parameters.txt");
paras_r8_d3 = readtable("Round 8 parameters.txt");
paras_r9_d3 = readtable("Round 9 parameters.txt");
paras_r10_d3 = readtable("Round 10 parameters.txt");
paras_r11_d3 = readtable("Round 11 parameters.txt");

%% Density plot of dn
dn_init = table2array(paras_init_d3(:,2));
dn_r2 = table2array(paras_r2_d3(:,2));
dn_r3 = table2array(paras_r3_d3(:,2));
dn_r4 = table2array(paras_r4_d3(:,2));
dn_r5 = table2array(paras_r5_d3(:,2));
dn_r6 = table2array(paras_r6_d3(:,2));
dn_r7 = table2array(paras_r7_d3(:,2));
dn_r8 = table2array(paras_r8_d3(:,2));
dn_r9 = table2array(paras_r9_d3(:,2));
dn_r10 = table2array(paras_r10_d3(:,2));
dn_r11 = table2array(paras_r11_d3(:,2));

[g,xii] = ksdensity(dn_r2,'Bandwidth',0.000006);
[h,xiii] = ksdensity(dn_r3,'Bandwidth',0.000006);
[j,xiv] = ksdensity(dn_r4,'Bandwidth',0.000006);
[z,xv] = ksdensity(dn_r5,'Bandwidth',0.000006);
[m,xvi] = ksdensity(dn_r6,'Bandwidth',0.000006);

figure
yline(1/(0.02 - 0.000069),'b','Linewidth',3.5);
xlim([0.000069 0.005]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(dn_init),0,'xb','markersize',20)
plot(mean(dn_r2),0,'xg','markersize',20)
plot(mean(dn_r3),0,'xr','markersize',20)
plot(mean(dn_r4),0,'xc','markersize',20)
plot(mean(dn_r5),0,'xm','markersize',20)
plot(mean(dn_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n}$)','Post-round 1 density($d_{n}$)','Post-round 2 density($d_{n}$)','Post-round 3 density($d_{n}$)','Post-round 4 density($d_{n}$)', 'Post-round 5 density($d_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$d_{n}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(dn_r7,'Bandwidth',0.000006);
[l,xviii] = ksdensity(dn_r8,'Bandwidth',0.000006);
[q,xix] = ksdensity(dn_r9,'Bandwidth',0.000006);
[w,xx] = ksdensity(dn_r10,'Bandwidth',0.000006);
[e,xxi] = ksdensity(dn_r11,'Bandwidth',0.000006);

figure
yline(1/(0.02 - 0.000069),'b','Linewidth',3.5);
xlim([0.000069 0.0005]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(dn_r7),0,'xg','markersize',20)
plot(mean(dn_r8),0,'xr','markersize',20)
plot(mean(dn_r9),0,'xc','markersize',20)
plot(mean(dn_r10),0,'xm','markersize',20)
plot(mean(dn_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n}$)','Post-round 6 density($d_{n}$)','Post-round 7 density($d_{n}$)','Post-round 8 density($d_{n}$)','Post-round 9 density($d_{n}$)', 'Post-round 10 density($d_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$d_{n}$ values')
ylabel('Probability density')

%% Density plot of gamma
gamma_init = table2array(paras_init_d3(:,3));
gamma_r2 = table2array(paras_r2_d3(:,3));
gamma_r3 = table2array(paras_r3_d3(:,3));
gamma_r4 = table2array(paras_r4_d3(:,3));
gamma_r5 = table2array(paras_r5_d3(:,3));
gamma_r6 = table2array(paras_r6_d3(:,3));
gamma_r7 = table2array(paras_r7_d3(:,3));
gamma_r8 = table2array(paras_r8_d3(:,3));
gamma_r9 = table2array(paras_r9_d3(:,3));
gamma_r10 = table2array(paras_r10_d3(:,3));
gamma_r11 = table2array(paras_r11_d3(:,3));

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
lgd = legend({'Initial density($\gamma$)','Post-round 1 density($\gamma$)','Post-round 2 density($\gamma$)','Post-round 3 density($\gamma$)','Post-round 4 density($\gamma$)', 'Post-round 5 density($\gamma$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$\gamma$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(gamma_r7,'Bandwidth',0.0002);
[l,xviii] = ksdensity(gamma_r8,'Bandwidth',0.0002);
[q,xix] = ksdensity(gamma_r9,'Bandwidth',0.0002);
[w,xx] = ksdensity(gamma_r10,'Bandwidth',0.0002);
[e,xxi] = ksdensity(gamma_r11,'Bandwidth',0.0002);

figure
yline(1/(0.26 - 0.005),'b','Linewidth',3.5);
xlim([0.005 0.008]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(gamma_init),0,'xb','markersize',20)
plot(mean(gamma_r7),0,'xg','markersize',20)
plot(mean(gamma_r8),0,'xr','markersize',20)
plot(mean(gamma_r9),0,'xc','markersize',20)
plot(mean(gamma_r10),0,'xm','markersize',20)
plot(mean(gamma_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\gamma$)','Post-round 6 density($\gamma$)','Post-round 7 density($\gamma$)','Post-round 8 density($\gamma$)','Post-round 9 density($\gamma$)', 'Post-round 10 density($\gamma$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$\gamma$ values')
ylabel('Probability density')

%% Density plot of rn
rn_init = table2array(paras_init_d3(:,4));
rn_r2 = table2array(paras_r2_d3(:,4));
rn_r3 = table2array(paras_r3_d3(:,4));
rn_r4 = table2array(paras_r4_d3(:,4));
rn_r5 = table2array(paras_r5_d3(:,4));
rn_r6 = table2array(paras_r6_d3(:,4));
rn_r7 = table2array(paras_r7_d3(:,4));
rn_r8 = table2array(paras_r8_d3(:,4));
rn_r9 = table2array(paras_r9_d3(:,4));
rn_r10 = table2array(paras_r10_d3(:,4));
rn_r11 = table2array(paras_r11_d3(:,4));

[g,xii] = ksdensity(rn_r2,'Bandwidth',0.003);
[h,xiii] = ksdensity(rn_r3,'Bandwidth',0.003);
[j,xiv] = ksdensity(rn_r4,'Bandwidth',0.003);
[z,xv] = ksdensity(rn_r5,'Bandwidth',0.003);
[m,xvi] = ksdensity(rn_r6,'Bandwidth',0.003);

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
lgd = legend({'Initial density($r_{n}$)','Post-round 1 density($r_{n}$)','Post-round 2 density($r_{n}$)','Post-round 3 density($r_{n}$)','Post-round 4 density($r_{n}$)', 'Post-round 5 density($r_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$r_{n}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(rn_r7,'Bandwidth',0.001);
[l,xviii] = ksdensity(rn_r8,'Bandwidth',0.001);
[q,xix] = ksdensity(rn_r9,'Bandwidth',0.001);
[w,xx] = ksdensity(rn_r10,'Bandwidth',0.001);
[e,xxi] = ksdensity(rn_r11,'Bandwidth',0.001);

figure
yline(1/(0.08 - 0.0008),'b','Linewidth',3.5);
xlim([0.0008 0.04]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(rn_init),0,'xb','markersize',20)
plot(mean(rn_r7),0,'xg','markersize',20)
plot(mean(rn_r8),0,'xr','markersize',20)
plot(mean(rn_r9),0,'xc','markersize',20)
plot(mean(rn_r10),0,'xm','markersize',20)
plot(mean(rn_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n}$)','Post-round 6 density($r_{n}$)','Post-round 7 density($r_{n}$)','Post-round 8 density($r_{n}$)','Post-round 9 density($r_{n}$)', 'Post-round 10 density($r_{n}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$r_{n}$ values')
ylabel('Probability density')

%% Density plots of eta
eta_init = table2array(paras_init_d3(:,5));
eta_r2 = table2array(paras_r2_d3(:,5));
eta_r3 = table2array(paras_r3_d3(:,5));
eta_r4 = table2array(paras_r4_d3(:,5));
eta_r5 = table2array(paras_r5_d3(:,5));
eta_r6 = table2array(paras_r6_d3(:,5));
eta_r7 = table2array(paras_r7_d3(:,5));
eta_r8 = table2array(paras_r8_d3(:,5));
eta_r9 = table2array(paras_r9_d3(:,5));
eta_r10 = table2array(paras_r10_d3(:,5));
eta_r11 = table2array(paras_r11_d3(:,5));


[g,xii] = ksdensity(eta_r2,'Bandwidth',0.12);
[h,xiii] = ksdensity(eta_r3,'Bandwidth',0.12);
[j,xiv] = ksdensity(eta_r4,'Bandwidth',0.12);
[z,xv] = ksdensity(eta_r5,'Bandwidth',0.12);
[m,xvi] = ksdensity(eta_r6,'Bandwidth',0.12);

figure
yline(1/(18 - 7),'b','Linewidth',3.5);
xlim([12 18]);
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

[k,xvii] = ksdensity(eta_r7,'Bandwidth',0.12);
[l,xviii] = ksdensity(eta_r8,'Bandwidth',0.12);
[q,xix] = ksdensity(eta_r9,'Bandwidth',0.12);
[w,xx] = ksdensity(eta_r10,'Bandwidth',0.12);
[e,xxi] = ksdensity(eta_r11,'Bandwidth',0.12);

figure
yline(1/(18 - 7),'b','Linewidth',3.5);
xlim([12 18]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(eta_init),0,'xb','markersize',20)
plot(mean(eta_r7),0,'xg','markersize',20)
plot(mean(eta_r8),0,'xr','markersize',20)
plot(mean(eta_r9),0,'xc','markersize',20)
plot(mean(eta_r10),0,'xm','markersize',20)
plot(mean(eta_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\eta$)','Post-round 6 density($\eta$)','Post-round 7 density($\eta$)','Post-round 8 density($\eta$)','Post-round 9 density($\eta$)', 'Post-round 10 density($\eta$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\eta$ values')
ylabel('Probability density')

%% Density plots of dm
dm_init = table2array(paras_init_d3(:,6));
dm_r2 = table2array(paras_r2_d3(:,6));
dm_r3 = table2array(paras_r3_d3(:,6));
dm_r4 = table2array(paras_r4_d3(:,6));
dm_r5 = table2array(paras_r5_d3(:,6));
dm_r6 = table2array(paras_r6_d3(:,6));
dm_r7 = table2array(paras_r7_d3(:,6));
dm_r8 = table2array(paras_r8_d3(:,6));
dm_r9 = table2array(paras_r9_d3(:,6));
dm_r10 = table2array(paras_r10_d3(:,6));
dm_r11 = table2array(paras_r11_d3(:,6));

[g,xii] = ksdensity(dm_r2,'Bandwidth',0.005);
[h,xiii] = ksdensity(dm_r3,'Bandwidth',0.005);
[j,xiv] = ksdensity(dm_r4,'Bandwidth',0.005);
[z,xv] = ksdensity(dm_r5,'Bandwidth',0.005);
[m,xvi] = ksdensity(dm_r6,'Bandwidth',0.005);

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
lgd = legend({'Initial density($d_{m}$)','Post-round 1 density($d_{m}$)','Post-round 2 density($d_{m}$)','Post-round 3 density($d_{m}$)','Post-round 4 density($d_{m}$)', 'Post-round 5 density($d_{m}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$d_{m}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(dm_r7,'Bandwidth',0.005);
[l,xviii] = ksdensity(dm_r8,'Bandwidth',0.005);
[q,xix] = ksdensity(dm_r9,'Bandwidth',0.005);
[w,xx] = ksdensity(dm_r10,'Bandwidth',0.005);
[e,xxi] = ksdensity(dm_r11,'Bandwidth',0.005);

figure
yline(1/(0.033 - 0.0001),'b','Linewidth',3.5);
xlim([0.0001 0.033]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(dm_init),0,'xb','markersize',20)
plot(mean(dm_r7),0,'xg','markersize',20)
plot(mean(dm_r8),0,'xr','markersize',20)
plot(mean(dm_r9),0,'xc','markersize',20)
plot(mean(dm_r10),0,'xm','markersize',20)
plot(mean(dm_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($d_{m}$)','Post-round 6 density($d_{m}$)','Post-round 7 density($d_{m}$)','Post-round 8 density($d_{m}$)','Post-round 9 density($d_{m}$)', 'Post-round 10 density($d_{m}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',8);
xlabel('$d_{m}$ values')
ylabel('Probability density')

%% Density plots of alpha
alpha_init = table2array(paras_init_d3(:,7));
alpha_r2 = table2array(paras_r2_d3(:,7));
alpha_r3 = table2array(paras_r3_d3(:,7));
alpha_r4 = table2array(paras_r4_d3(:,7));
alpha_r5 = table2array(paras_r5_d3(:,7));
alpha_r6 = table2array(paras_r6_d3(:,7));
alpha_r7 = table2array(paras_r7_d3(:,7));
alpha_r8 = table2array(paras_r8_d3(:,7));
alpha_r9 = table2array(paras_r9_d3(:,7));
alpha_r10 = table2array(paras_r10_d3(:,7));
alpha_r11 = table2array(paras_r11_d3(:,7));

[g,xii] = ksdensity(alpha_r2,'Bandwidth',0.005);
[h,xiii] = ksdensity(alpha_r3,'Bandwidth',0.005);
[j,xiv] = ksdensity(alpha_r4,'Bandwidth',0.005);
[z,xv] = ksdensity(alpha_r5,'Bandwidth',0.005);
[m,xvi] = ksdensity(alpha_r6,'Bandwidth',0.005);

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

[k,xvii] = ksdensity(alpha_r7,'Bandwidth',0.005);
[l,xviii] = ksdensity(alpha_r8,'Bandwidth',0.005);
[q,xix] = ksdensity(alpha_r9,'Bandwidth',0.005);
[w,xx] = ksdensity(alpha_r10,'Bandwidth',0.005);
[e,xxi] = ksdensity(alpha_r11,'Bandwidth',0.005);

figure
yline(1/(0.18 - 0.07),'b','Linewidth',3.5);
xlim([0.07 0.18]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(alpha_init),0,'xb','markersize',20)
plot(mean(alpha_r7),0,'xg','markersize',20)
plot(mean(alpha_r8),0,'xr','markersize',20)
plot(mean(alpha_r9),0,'xc','markersize',20)
plot(mean(alpha_r10),0,'xm','markersize',20)
plot(mean(alpha_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($\alpha$)','Post-round 6 density($\alpha$)','Post-round 7 density($\alpha$)','Post-round 8 density($\alpha$)','Post-round 9 density($\alpha$)', 'Post-round 10 density($\alpha$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$\alpha$ values')
ylabel('Probability density')

%% Density plots of r.init
r_init_init = table2array(paras_init_d3(:,8));
r_init_r2 = table2array(paras_r2_d3(:,8));
r_init_r3 = table2array(paras_r3_d3(:,8));
r_init_r4 = table2array(paras_r4_d3(:,8));
r_init_r5 = table2array(paras_r5_d3(:,8));
r_init_r6 = table2array(paras_r6_d3(:,8));
r_init_r7 = table2array(paras_r7_d3(:,8));
r_init_r8 = table2array(paras_r8_d3(:,8));
r_init_r9 = table2array(paras_r9_d3(:,8));
r_init_r10 = table2array(paras_r10_d3(:,8));
r_init_r11 = table2array(paras_r11_d3(:,8));

[g,xii] = ksdensity(r_init_r2,'Bandwidth',0.1);
[h,xiii] = ksdensity(r_init_r3,'Bandwidth',0.1);
[j,xiv] = ksdensity(r_init_r4,'Bandwidth',0.1);
[z,xv] = ksdensity(r_init_r5,'Bandwidth',0.1);
[m,xvi] = ksdensity(r_init_r6,'Bandwidth',0.1);

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

[k,xvii] = ksdensity(r_init_r7,'Bandwidth',0.1);
[l,xviii] = ksdensity(r_init_r8,'Bandwidth',0.1);
[q,xix] = ksdensity(r_init_r9,'Bandwidth',0.1);
[w,xx] = ksdensity(r_init_r10,'Bandwidth',0.1);
[e,xxi] = ksdensity(r_init_r11,'Bandwidth',0.1);

figure
yline(1/(5 - 1),'b','Linewidth',3.5);
xlim([1 5]);
%ylim([0 13]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(r_init_init),0,'xb','markersize',20)
plot(mean(r_init_r7),0,'xg','markersize',20)
plot(mean(r_init_r8),0,'xr','markersize',20)
plot(mean(r_init_r9),0,'xc','markersize',20)
plot(mean(r_init_r10),0,'xm','markersize',20)
plot(mean(r_init_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($R_{init.}$)','Post-round 6 density($R_{init.}$)','Post-round 7 density($R_{init.}$)','Post-round 8 density($R_{init.}$)','Post-round 9 density($R_{init.}$)', 'Post-round 10 density($R_{init.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$R_{init.}$ values')
ylabel('Probability density')

%% Density plots of p.ext
p_ext_init = table2array(paras_init_d3(:,9));
p_ext_r2 = table2array(paras_r2_d3(:,9));
p_ext_r3 = table2array(paras_r3_d3(:,9));
p_ext_r4 = table2array(paras_r4_d3(:,9));
p_ext_r5 = table2array(paras_r5_d3(:,9));
p_ext_r6 = table2array(paras_r6_d3(:,9));
p_ext_r7 = table2array(paras_r7_d3(:,9));
p_ext_r8 = table2array(paras_r8_d3(:,9));
p_ext_r9 = table2array(paras_r9_d3(:,9));
p_ext_r10 = table2array(paras_r10_d3(:,9));
p_ext_r11 = table2array(paras_r11_d3(:,9));

[g,xii] = ksdensity(p_ext_r2,'Bandwidth',0.005);
[h,xiii] = ksdensity(p_ext_r3,'Bandwidth',0.005);
[j,xiv] = ksdensity(p_ext_r4,'Bandwidth',0.005);
[z,xv] = ksdensity(p_ext_r5,'Bandwidth',0.005);
[m,xvi] = ksdensity(p_ext_r6,'Bandwidth',0.005);

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
lgd = legend({'Initial density($P_{ext.}$)','Post-round 1 density($P_{ext.}$)','Post-round 2 density($P_{ext.}$)','Post-round 3 density($P_{ext.}$)','Post-round 4 density($P_{ext.}$)', 'Post-round 5 density($P_{ext.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

[k,xvii] = ksdensity(p_ext_r7,'Bandwidth',0.005);
[l,xviii] = ksdensity(p_ext_r8,'Bandwidth',0.005);
[q,xix] = ksdensity(p_ext_r9,'Bandwidth',0.005);
[w,xx] = ksdensity(p_ext_r10,'Bandwidth',0.005);
[e,xxi] = ksdensity(p_ext_r11,'Bandwidth',0.005);

figure
yline(1/(0.1 - 0.01),'b','Linewidth',3.5);
xlim([0.01 0.1]);
%ylim([0 13]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(p_ext_init),0,'xb','markersize',20)
plot(mean(p_ext_r7),0,'xg','markersize',20)
plot(mean(p_ext_r8),0,'xr','markersize',20)
plot(mean(p_ext_r9),0,'xc','markersize',20)
plot(mean(p_ext_r10),0,'xm','markersize',20)
plot(mean(p_ext_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($P_{ext.}$)','Post-round 6 density($P_{ext.}$)','Post-round 7 density($P_{ext.}$)','Post-round 8 density($P_{ext.}$)','Post-round 9 density($P_{ext.}$)', 'Post-round 10 density($P_{ext.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

%% Density plots of p.mit
p_mit_init = table2array(paras_init_d3(:,10));
p_mit_r2 = table2array(paras_r2_d3(:,10));
p_mit_r3 = table2array(paras_r3_d3(:,10));
p_mit_r4 = table2array(paras_r4_d3(:,10));
p_mit_r5 = table2array(paras_r5_d3(:,10));
p_mit_r6 = table2array(paras_r6_d3(:,10));
p_mit_r7 = table2array(paras_r7_d3(:,10));
p_mit_r8 = table2array(paras_r8_d3(:,10));
p_mit_r9 = table2array(paras_r9_d3(:,10));
p_mit_r10 = table2array(paras_r10_d3(:,10));
p_mit_r11 = table2array(paras_r11_d3(:,10));


[g,xii] = ksdensity(p_mit_r2,'Bandwidth',0.03);
[h,xiii] = ksdensity(p_mit_r3,'Bandwidth',0.03);
[j,xiv] = ksdensity(p_mit_r4,'Bandwidth',0.03);
[z,xv] = ksdensity(p_mit_r5,'Bandwidth',0.03);
[m,xvi] = ksdensity(p_mit_r6,'Bandwidth',0.03);

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

[k,xvii] = ksdensity(p_mit_r7,'Bandwidth',0.03);
[l,xviii] = ksdensity(p_mit_r8,'Bandwidth',0.03);
[q,xix] = ksdensity(p_mit_r9,'Bandwidth',0.03);
[w,xx] = ksdensity(p_mit_r10,'Bandwidth',0.03);
[e,xxi] = ksdensity(p_mit_r11,'Bandwidth',0.03);

figure
yline(1/(1 - 0.2),'b','Linewidth',3.5);
xlim([0.2 1]);
%ylim([0 13]);
hold on;
plot(xvii,k,'g-',xviii,l,'r-',xix,q,'c-', xx, w, 'm-', xxi, e, 'k-')
plot(mean(p_mit_init),0,'xb','markersize',20)
plot(mean(p_mit_r7),0,'xg','markersize',20)
plot(mean(p_mit_r8),0,'xr','markersize',20)
plot(mean(p_mit_r9),0,'xc','markersize',20)
plot(mean(p_mit_r10),0,'xm','markersize',20)
plot(mean(p_mit_r11),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit.}$)','Post-round 6 density($P_{mit.}$)','Post-round 7 density($P_{mit.}$)','Post-round 8 density($P_{mit.}$)','Post-round 9 density($P_{mit.}$)', 'Post-round 10 density($P_{mit.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit.}$ values')
ylabel('Probability density')
