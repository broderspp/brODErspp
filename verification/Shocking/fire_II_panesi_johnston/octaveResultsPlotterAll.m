close all
clear
clc

% Load the solution
xT_Panesi1  = load('data_reference/T_pan_shocking1.dat');
xTv_Panesi1 = load('data_reference/Tv_pan_shocking1.dat');

xT_Panesi2  = load('data_reference/T_pan_shocking2.dat');
xTv_Panesi2 = load('data_reference/Tv_pan_shocking2.dat');

xT_Johnston  = load('data_reference/T_johnston.dat');
xTv_Johnston = load('data_reference/Tv_johnston.dat');

% Load shocking result
dd_S = load('shocking_output.dat');

xx_S = dd_S(:,1);
Xi_S = dd_S(:,2:end-3);
uu_S = dd_S(:,end-2);
TT_S = dd_S(:,end-1);
Tv_S = dd_S(:,end);

% PLOT
figure
hold on
% translational temperatures:
plot(xT_Panesi1(:,1), xT_Panesi1(:,2), '-k', 'linewidth', 1)
plot(xT_Panesi2(:,1), xT_Panesi2(:,2), '-g', 'linewidth', 1)
plot(xT_Johnston(:,1), xT_Johnston(:,2), '-m', 'linewidth', 1)
plot(xx_S, TT_S, 'b', 'linewidth', 2)
% vibrational temperatures
plot(xTv_Panesi1(:,1), xTv_Panesi1(:,2), '-k', 'linewidth', 1)
plot(xTv_Panesi2(:,1), xTv_Panesi2(:,2), '-g', 'linewidth', 1)
plot(xTv_Johnston(:,1), xTv_Johnston(:,2), '-m', 'linewidth', 1)
plot(xx_S, Tv_S, 'b', 'linewidth', 2)

legend('Panesi Shocking 1', 'Panesi - Shocking 1', 'Johnston', 'CURRENT - SHOCKING++')

ylabel('Temperatures [K]')
xlabel('distance from shock [m]')

print('Temp.png')


fprintf('=====>  Press RETURN to quit  <=====\n')
pause()
