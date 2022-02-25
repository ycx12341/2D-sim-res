% Parameter_estimations.m
% Author: Yunchen Xiao

% This MATLAB file generates the plots of actual effective sample size
% (ESS) and bandwidth factors for resampling weights calculation for each  
% round of the simulations in the case study on T98G glioma invasion 
% pattern. 
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
ess_bw_run1 = readtable('Run 1 ESS BW.txt');
ess_bw_run2 = readtable('Run 2 ESS BW.txt');
ess_bw_run3 = readtable('Run 3 ESS BW.txt');

%% ESS
ess_run1 = table2array(ess_bw_run1(:,2));
ess_run2 = table2array(ess_bw_run2(:,2));
ess_run3 = table2array(ess_bw_run3(:,2)); 

%% BW factors for resampling weights calculation
bw_run1 = table2array(ess_bw_run1(:,3));
bw_run2 = table2array(ess_bw_run2(:,3));
bw_run3 = table2array(ess_bw_run3(:,3)); 

x = [1 2 3 4 5 6];

%% Plot of actual ESS variations in three different runs.
figure
title('Variation of actual ESS (glioma)')
hold on;
plot(1, ess_run1(1), 'xk', 'markersize', 20)
plot(2, ess_run1(2), 'xk', 'markersize', 20)
plot(3, ess_run1(3), 'xk', 'markersize', 20)
plot(4, ess_run1(4), 'xk', 'markersize', 20)
plot(5, ess_run1(5), 'xk', 'markersize', 20)
plot(6, ess_run1(6), 'xk', 'markersize', 20)
plot(1, ess_run2(1), 'xb', 'markersize', 20)
plot(2, ess_run2(2), 'xb', 'markersize', 20)
plot(3, ess_run2(3), 'xb', 'markersize', 20)
plot(4, ess_run2(4), 'xb', 'markersize', 20)
plot(5, ess_run2(5), 'xb', 'markersize', 20)
plot(6, ess_run2(6), 'xb', 'markersize', 20)
plot(1, ess_run3(1), 'xr', 'markersize', 20)
plot(2, ess_run3(2), 'xr', 'markersize', 20)
plot(3, ess_run3(3), 'xr', 'markersize', 20)
plot(4, ess_run3(4), 'xr', 'markersize', 20)
plot(5, ess_run3(5), 'xr', 'markersize', 20)
plot(6, ess_run3(6), 'xr', 'markersize', 20)
yline(1500,'k--','Linewidth',3.5);
xlim([0.5 6.5])
xlabel('Rounds')
ylim([1400 2600])
ylabel('Actual Effective sample size')

%% Plot of BW factor variations in three different runs
figure
title('Variation of bandwidth factors (glioma)')
hold on;
plot(1, bw_run1(1), 'xk', 'markersize', 20)
plot(2, bw_run1(2), 'xk', 'markersize', 20)
plot(3, bw_run1(3), 'xk', 'markersize', 20)
plot(4, bw_run1(4), 'xk', 'markersize', 20)
plot(5, bw_run1(5), 'xk', 'markersize', 20)
plot(6, bw_run1(6), 'xk', 'markersize', 20)
plot(1, bw_run2(1), 'xb', 'markersize', 20)
plot(2, bw_run2(2), 'xb', 'markersize', 20)
plot(3, bw_run2(3), 'xb', 'markersize', 20)
plot(4, bw_run2(4), 'xb', 'markersize', 20)
plot(5, bw_run2(5), 'xb', 'markersize', 20)
plot(6, bw_run2(6), 'xb', 'markersize', 20)
plot(1, bw_run3(1), 'xr', 'markersize', 20)
plot(2, bw_run3(2), 'xr', 'markersize', 20)
plot(3, bw_run3(3), 'xr', 'markersize', 20)
plot(4, bw_run3(4), 'xr', 'markersize', 20)
plot(5, bw_run3(5), 'xr', 'markersize', 20)
plot(6, bw_run3(6), 'xr', 'markersize', 20)
xlim([0.5 6.5])
xlabel('Rounds')
ylim([2 11])
ylabel('Bandwidth factors')