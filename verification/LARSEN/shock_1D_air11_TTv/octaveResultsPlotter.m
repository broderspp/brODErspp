%
% This script reads an "external flowfield file" computed by Stagline
% and compares with results from LARSEN
%

close all
clear
clc

page_screen_output(0);
graphics_toolkit('gnuplot');

% Parameters
nSpecies = 11; % Number of species

% Load external stagline flowfield file
dd = load('output.dat');

% Load LARSEN result
ddL = load('outpL');

% Extract info..
xx  = dd(:,1);
TT  = dd(:,2);
Tv  = dd(:,3);
rho = dd(:,4);
uu  = dd(:,5);

yi = dd(:,6:end);


% PLOT
figure
hold on
% T
plot(xx, TT, '+-k', 'linewidth',1)
plot(ddL(:,1), ddL(:,end-1), 'b', 'linewidth',2)
% Tv
plot(xx, Tv, '+-k', 'linewidth',1)
plot(ddL(:,1), ddL(:,end), 'r', 'linewidth',2)
xlabel('x [m]', 'fontsize', 14)
ylabel('T [K]', 'fontsize', 14)
legend('Stagline', 'LARSEN')
title('RED: stagline, BLUE: LARSEN results', 'fontsize', 14)

print('Temp.png');

figure
plot(xx, yi, '-r', 'linewidth', 2)
hold on
plot(ddL(:,1), ddL(:,2:end-2), 'b', 'linewidth', 2)
xlabel('x [m]', 'fontsize', 14)
ylabel('Yi', 'fontsize', 14)
title('RED: stagline, BLUE: LARSEN results', 'fontsize', 14)

print('Yi.png')

% Block the execution to this point
fprintf('\n##########################\n\n=========> Press RETURN to quit <=========\n')
pause()
