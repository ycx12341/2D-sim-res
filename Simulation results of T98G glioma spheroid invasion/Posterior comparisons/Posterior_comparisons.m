% Posterior_comparisons.m
% Author: Yunchen Xiao
% This MATLAB file generates the plots of posterior probability densities 
% of the parameter estimates obtained using non-error and error-calibrated
% ABC on the T98G glioma invasion pattern dataset.

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
post_ec = readtable("Round 6 parameters log transform.txt");
post_non_ec = readtable("Round 7 parameters.txt");

dn_ec = table2array(post_ec(:,2));
rn_ec = table2array(post_ec(:,3));
r_init_ec = table2array(post_ec(:,4));
p_ext_ec = table2array(post_ec(:,5));
p_mit_ec = table2array(post_ec(:,6));

dn_non_ec = table2array(post_non_ec(:,2));
rn_non_ec = table2array(post_non_ec(:,3));
r_init_non_ec = table2array(post_non_ec(:,4));
p_ext_non_ec = table2array(post_non_ec(:,5));
p_mit_non_ec = table2array(post_non_ec(:,6));

%% dn posteriors
[g,xii] = ksdensity(dn_ec, 'bandwidth', 0.00001);
[h,xiii] = ksdensity(dn_non_ec, 'bandwidth', 0.00001);

figure
xlim([0.000069 0.0003]);
hold on;
plot(xii,g,'k-',xiii,h,'r-')
plot(mean(dn_ec),0,'xk','markersize',20)
plot(mean(dn_non_ec),0,'xr','markersize',20)
hold off;
lgd = legend({'Posterior density (EC)','Posterior density (non-EC)'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$d_{n}$ values')
ylabel('Probability density')

%% rn posteriors
[s,xii] = ksdensity(rn_ec, 'bandwidth', 0.0006);
[d,xiii] = ksdensity(rn_non_ec, 'bandwidth', 0.0006);

figure
xlim([0 0.0035]);
hold on;
plot(xii,s,'k-',xiii,d,'r-')
plot(mean(rn_ec),0,'xk','markersize',20)
plot(mean(rn_non_ec),0,'xr','markersize',20)
hold off;
lgd = legend({'Posterior density (EC)','Posterior density (non-EC)'},'Location','northeast','Orientation','vertical','Fontsize',12);
xlabel('$r_{n}$ values')
ylabel('Probability density')

%% r_init posteriors
[l,xii] = ksdensity(r_init_ec, 'bandwidth', 0.008);
[k,xiii] = ksdensity(r_init_non_ec, 'bandwidth', 0.008);

figure
xlim([0.01 0.15]);
hold on;
plot(xii,l,'k-',xiii,k,'r-')
plot(mean(r_init_ec),0,'xk','markersize',20)
plot(mean(r_init_non_ec),0,'xr','markersize',20)
hold off;
lgd = legend({'Posterior density (EC)','Posterior density (non-EC)'},'Location','northwest','Orientation','vertical','Fontsize',12);
xlabel('$R_{init.}$ values')
ylabel('Probability density')

%% P_ext posteriors
[o,xii] = ksdensity(p_ext_ec, 'bandwidth', 0.02);
[p,xiii] = ksdensity(p_ext_non_ec, 'bandwidth', 0.02);

figure
xlim([0.01 0.1]);
hold on;
plot(xii,o,'k-',xiii,p,'r-')
plot(mean(p_ext_ec),0,'xk','markersize',20)
plot(mean(p_ext_non_ec),0,'xr','markersize',20)
hold off;
lgd = legend({'Posterior density (EC)','Posterior density (non-EC)'},'Location','northeast','Orientation','vertical','Fontsize',9);
xlabel('$P_{ext.}$ values')
ylabel('Probability density')

%% P_mit posteriors
[n,xii] = ksdensity(p_mit_ec, 'bandwidth', 0.07);
[m,xiii] = ksdensity(p_mit_non_ec, 'bandwidth', 0.07);

figure
xlim([0.02 1]);
hold on;
plot(xii,n,'k-',xiii,m,'r-')
plot(mean(p_mit_ec),0,'xk','markersize',20)
plot(mean(p_mit_non_ec),0,'xr','markersize',20)
hold off;
lgd = legend({'Posterior density (EC)','Posterior density (non-EC)'},'Location','northeast','Orientation','vertical','Fontsize',10);
xlabel('$P_{mit.}$ values')
ylabel('Probability density')