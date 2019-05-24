close all
clear
clc

pkg load optim;

page_screen_output(0)

% ==========  Freestream values  ===========

T_fs = 220;  % [K]

% Load data file storing initial condition
%dd = load('TOTALLY_INITIAL.csv');
dd = load('FINAL_vertExt.csv');

% Extract it
yy = dd(:,end-1);
TT = dd(:,1);
ni = dd(:,[2:12]);

n = sum(ni, 2);

Xi = ni./repmat(n, 1, size(ni,2));

% Crop solution
pos_y_stop = find(yy > 0.03);
pos_y_stop = pos_y_stop(1);

pos_y_start = find(yy > 0.022);
pos_y_start = pos_y_start(1);

% Plot solution
figure
plot(yy, TT)
hold on
plot(yy(pos_y_start:pos_y_stop), TT(pos_y_start:pos_y_stop), 'g')

% Fit temperature
T_funct   = @(p, x) T_fs + p(1)*exp(p(2)*(x-yy(pos_y_start))); % Temp model function

coefs_T = nonlin_curvefit(T_funct, [6000;-3], yy(pos_y_start:pos_y_stop), TT(pos_y_start:pos_y_stop));

y_extended = linspace(yy(pos_y_start), yy(end));

hold on
plot(y_extended, T_funct(coefs_T, y_extended), 'r', 'linewidth', 2)

% Create one only solution
T_ext = T_funct(coefs_T, yy(pos_y_stop:end));

T_new = [TT(1:pos_y_stop-1); T_ext];

figure
plot(yy, T_new, 'k', 'linewidth', 2)

[yy, T_new]
