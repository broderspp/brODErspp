close all
clear
clc

% Load a streamline
str_num = 3; % PICK HERE THE STREAMLINE <<<----------------------------------------

filename_out = sprintf('output_streamlines/streamline_%05d', str_num);
fprintf('Reading streamline %s\n');

dd = load(filename_out);

% Extract data
xx = dd(:,1);
yy = dd(:,2);
TT = dd(:,3);
rho = dd(:,4);
uu  = dd(:,5);
vv  = dd(:,6);


% PLOT information
figure
subplot(5,1,1)
plot(xx, yy, 'k', 'linewidth', 2)
xlabel('x [m]')
ylabel('y [m]')

subplot(5,1,2)
plot(xx, TT, 'r', 'linewidth', 2)
ylabel('T [K]')

subplot(5,1,3)
plot(xx, rho, 'g', 'linewidth', 2)
ylabel('rho [kg/m3]')

subplot(5,1,4)
plot(xx, uu, 'b', 'linewidth', 2)
ylabel('u [m/s]')

subplot(5,1,5)
plot(xx, vv, 'm', 'linewidth', 2)
ylabel('v [m/s]')
xlabel('x [m]')

figure
plot(xx, dd(:,7:end), 'k')
xlabel('x [m]')
ylabel('Xi')

