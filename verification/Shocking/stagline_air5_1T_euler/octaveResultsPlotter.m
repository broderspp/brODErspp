close all
clear
clc

% Load the solution
dd = load('stagline_result.dat');

xx = dd(end:-1:1,1);
uu = -dd(end:-1:1,2); % U is inverted in stagline
TT = dd(end:-1:1,4);
Xi = dd(end:-1:1,9:13);

% Load shocking result
dd_S = load('shocking_output.dat');

xx_S = dd_S(:,1);
Xi_S = dd_S(:,2:end-2);
uu_S = dd_S(:,end-1);
TT_S = dd_S(:,end);

% Detect position where the shock starts
posStart = find(TT > 1.1*TT(1));
posStart = posStart(3);

% Rescale position
xx = xx(posStart)-xx;

% PLOT
figure
hold on
plot(xx, TT, 'om', 'linewidth', 2)
plot(xx_S, TT_S, 'b', 'linewidth', 2)
legend('Stagline', 'CURRENT - Shocking')
ylabel('Temperature [K]')
xlabel('x [m]')

print('Temp.png')

% ---
figure
hold on
semilogy(xx, Xi, 'm', 'linewidth', 2)
semilogy(xx_S, Xi_S, 'b', 'linewidth', 2)
title('GREEN: Stagline - BLUE: Shocking')
ylabel('Mole fractions')
xlabel('x [m]')
ylim([1e-10, 1])

print('Xi.png')

% ---
figure
hold on
plot(xx, uu, 'm', 'linewidth', 2)
plot(xx_S, uu_S, 'b', 'linewidth', 2)
legend('Stagline', 'CURRENT - Shocking')
ylabel('velocity')
xlabel('x [m]')

print('Velocity.png')

fprintf('=====>  Press RETURN to quit  <=====\n')
pause()
