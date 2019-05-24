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
nSpecies = 5; % Number of species

% Load external stagline flowfield file
dd = load('output.dat');

% Load LARSEN result
ddL = load('outpL');

% Extract info..
xx  = dd(:,1);
TT  = dd(:,2);
rho = dd(:,3);
uu  = dd(:,4);

yi = dd(:,5:5+nSpecies-1);
rhoi = repmat(rho, 1, size(yi,2)).*yi;

% PLOT
figure
plot(xx, TT, '+-r', 'linewidth',1)
hold on
plot(ddL(:,1), ddL(:,end), 'b', 'linewidth',2)
xlabel('x [m]', 'fontsize', 14)
ylabel('T [K]', 'fontsize', 14)
legend('Stagline', 'LARSEN')
title('RED: stagline, BLUE: LARSEN results', 'fontsize', 14)

print('Temp.png');

figure
plot(xx, yi, '-r', 'linewidth', 2)
hold on
plot(ddL(:,1), ddL(:,2:end-1), 'b', 'linewidth', 2)
xlabel('x [m]', 'fontsize', 14)
ylabel('Yi', 'fontsize', 14)
title('RED: stagline, BLUE: LARSEN results', 'fontsize', 14)

print('Yi.png')

% Block the execution to this point
fprintf('\n##########################\n\n=========> Press RETURN to quit <=========\n')
pause()
