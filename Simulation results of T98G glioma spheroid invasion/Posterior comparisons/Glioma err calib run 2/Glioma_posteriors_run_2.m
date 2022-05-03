% Glioma_posteriors_run2.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of probability densities of the
% parameter estimations at the end of each round of the second run of 
% applying error-calibrated ABC scheme on the invasion pattern of the SCC 
% reference dataset.

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
paras_init = readtable("Round 1 initial parameters.txt");
paras_r2 = readtable("Round 2 parameters log transform.txt");
paras_r3 = readtable("Round 3 parameters log transform.txt");
paras_r4 = readtable("Round 4 parameters log transform.txt");
paras_r5 = readtable("Round 5 parameters log transform.txt");
paras_r6 = readtable("Round 6 parameters log transform.txt");


%% Density plot of dn
dn_init = table2array(paras_init(:,2));
dn_r2 = table2array(paras_r2(:,2));
dn_r3 = table2array(paras_r3(:,2));
dn_r4 = table2array(paras_r4(:,2));
dn_r5 = table2array(paras_r5(:,2));
dn_r6 = table2array(paras_r6(:,2));

[g,xii] = ksdensity(dn_r2,'Bandwidth',0.00001);
[h,xiii] = ksdensity(dn_r3,'Bandwidth',0.00001);
[j,xiv] = ksdensity(dn_r4,'Bandwidth',0.00001);
[z,xv] = ksdensity(dn_r5,'Bandwidth',0.00001);
[m,xvi] = ksdensity(dn_r6,'Bandwidth',0.00001);

figure
yline(1/(0.02 - 0.000069),'b','Linewidth',3.5);
xlim([0.000069 0.0005]);
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

figure
yline(1/(0.02 - 0.000069),'b','Linewidth',3.5);
xlim([0.000069 0.0005]);
hold on;
plot(xvi,m, 'k-')
plot(mean(dn_init),0,'xb','markersize',20)
plot(mean(dn_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($d_{n}$)', 'Post-round 5 density($d_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$d_{n}$ values')
ylabel('Probability density')

%% Density plot of rn
rn_init = table2array(paras_init(:,3));
rn_r2 = table2array(paras_r2(:,3));
rn_r3 = table2array(paras_r3(:,3));
rn_r4 = table2array(paras_r4(:,3));
rn_r5 = table2array(paras_r5(:,3));
rn_r6 = table2array(paras_r6(:,3));

[g,xii] = ksdensity(rn_r2,'Bandwidth',0.0006);
[h,xiii] = ksdensity(rn_r3,'Bandwidth',0.0006);
[j,xiv] = ksdensity(rn_r4,'Bandwidth',0.0006);
[z,xv] = ksdensity(rn_r5,'Bandwidth',0.0006);
[m,xvi] = ksdensity(rn_r6,'Bandwidth',0.0006);

figure
yline(1/(0.0035 - 0),'b','Linewidth',3.5);
xlim([0 0.0035]);
hold on;
plot(xii,g,'g-',xiii,h,'r-',xiv,j,'c-', xv, z, 'm-', xvi, m, 'k-')
plot(mean(rn_init),0,'xb','markersize',20)
plot(mean(rn_r2),0,'xg','markersize',20)
plot(mean(rn_r3),0,'xr','markersize',20)
plot(mean(rn_r4),0,'xc','markersize',20)
plot(mean(rn_r5),0,'xm','markersize',20)
plot(mean(rn_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n}$)','Post-round 1 density($r_{n}$)','Post-round 2 density($r_{n}$)','Post-round 3 density($r_{n}$)','Post-round 4 density($r_{n}$)', 'Post-round 5 density($r_{n}$)'},'Interpreter','latex','Location','southwest','Orientation','vertical','Fontsize',10);
xlabel('$r_{n}$ values')
ylabel('Probability density')

figure
yline(1/(0.0035 - 0),'b','Linewidth',3.5);
xlim([0 0.0035]);
hold on;
plot(xvi,m, 'k-')
plot(mean(rn_init),0,'xb','markersize',20)
plot(mean(rn_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($r_{n}$)','Post-round 5 density($r_{n}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$r_{n}$ values')
ylabel('Probability density')

%% Density plots of r.init
r_init_init = table2array(paras_init(:,4));
r_init_r2 = table2array(paras_r2(:,4));
r_init_r3 = table2array(paras_r3(:,4));
r_init_r4 = table2array(paras_r4(:,4));
r_init_r5 = table2array(paras_r5(:,4));
r_init_r6 = table2array(paras_r6(:,4));

[g,xii] = ksdensity(r_init_r2,'Bandwidth',0.008);
[h,xiii] = ksdensity(r_init_r3,'Bandwidth',0.008);
[j,xiv] = ksdensity(r_init_r4,'Bandwidth',0.008);
[z,xv] = ksdensity(r_init_r5,'Bandwidth',0.008);
[m,xvi] = ksdensity(r_init_r6,'Bandwidth',0.008);

figure
yline(1/(0.15 - 0.01),'b','Linewidth',3.5);
xlim([0.01 0.15]);
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
lgd = legend({'Initial density($R_{init.}$)','Post-round 1 density($R_{init.}$)','Post-round 2 density($R_{init.}$)','Post-round 3 density($R_{init.}$)','Post-round 4 density($R_{init.}$)', 'Post-round 5 density($R_{init.}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$R_{init.}$ values')
ylabel('Probability density')

figure
yline(1/(0.15 - 0.01),'b','Linewidth',3.5);
xlim([0.01 0.25]);
%ylim([0 13]);
hold on;
plot(xvi,m, 'k-')
plot(mean(r_init_init),0,'xb','markersize',20)
plot(mean(r_init_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($R_{init.}$)','Post-round 5 density($R_{init.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$R_{init.}$ values')
ylabel('Probability density')

%% Density plots of p.ext
p_ext_init = table2array(paras_init(:,5));
p_ext_r2 = table2array(paras_r2(:,5));
p_ext_r3 = table2array(paras_r3(:,5));
p_ext_r4 = table2array(paras_r4(:,5));
p_ext_r5 = table2array(paras_r5(:,5));
p_ext_r6 = table2array(paras_r6(:,5));

[g,xii] = ksdensity(p_ext_r2,'Bandwidth',0.013);
[h,xiii] = ksdensity(p_ext_r3,'Bandwidth',0.013);
[j,xiv] = ksdensity(p_ext_r4,'Bandwidth',0.013);
[z,xv] = ksdensity(p_ext_r5,'Bandwidth',0.013);
[m,xvi] = ksdensity(p_ext_r6,'Bandwidth',0.013);

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
lgd = legend({'Initial density($P_{ext.}$)','Post-round 1 density($P_{ext.}$)','Post-round 2 density($P_{ext.}$)','Post-round 3 density($P_{ext.}$)','Post-round 4 density($P_{ext.}$)', 'Post-round 5 density($P_{ext.}$)'},'Interpreter','latex','Location','southwest','Orientation','vertical','Fontsize',10);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

figure
yline(1/(0.1 - 0.01),'b','Linewidth',3.5);
xlim([0.01 0.1]);
ylim([0 20]);
hold on;
plot(xvi,m, 'k-')
plot(mean(p_ext_init),0,'xb','markersize',20)
plot(mean(p_ext_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($P_{ext.}$)','Post-round 5 density($P_{ext.}$)'},'Interpreter','latex','Location','northwest','Orientation','vertical','Fontsize',10);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

%% Density plots of p.mit
p_mit_init = table2array(paras_init(:,6));
p_mit_r2 = table2array(paras_r2(:,6));
p_mit_r3 = table2array(paras_r3(:,6));
p_mit_r4 = table2array(paras_r4(:,6));
p_mit_r5 = table2array(paras_r5(:,6));
p_mit_r6 = table2array(paras_r6(:,6));

[g,xii] = ksdensity(p_mit_r2,'Bandwidth',0.07);
[h,xiii] = ksdensity(p_mit_r3,'Bandwidth',0.07);
[j,xiv] = ksdensity(p_mit_r4,'Bandwidth',0.07);
[z,xv] = ksdensity(p_mit_r5,'Bandwidth',0.07);
[m,xvi] = ksdensity(p_mit_r6,'Bandwidth',0.07);

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

figure
yline(1/(1 - 0.2),'b','Linewidth',3.5);
xlim([0.2 1]);
%ylim([0 13]);
hold on;
plot(xvi,m, 'k-')
plot(mean(p_mit_init),0,'xb','markersize',20)
plot(mean(p_mit_r6),0,'xk','markersize',20)
hold off;
lgd = legend({'Initial density($P_{mit.}$)','Post-round 5 density($P_{mit.}$)'},'Interpreter','latex','Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit.}$ values')
ylabel('Probability density')
