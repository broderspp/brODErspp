close all
clear
clc

% Load data
dd = load('data.csv');

x = dd(:,end-2);
y = dd(:,end-1);
T = dd(:,1);

% Diffusion coefficient
alpha = 0.0749138/1142.4;

% PLOT

eta = linspace(0,2);

figure
hold on
plot(erf(eta), eta, 'r', 'linewidth', 2)
plot((T-T(1))/(T(end)-T(1)), y/sqrt(4*alpha*x(end)), 'xb', 'linewidth', 2)
xlabel('\theta')
ylabel('\eta')
