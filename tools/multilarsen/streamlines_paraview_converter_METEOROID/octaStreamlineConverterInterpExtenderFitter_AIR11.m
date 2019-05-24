%%%%%
%%%%%
%%%%%   P R O P E R T I E S     O F    T H E    A I R 11     M I X T U R E
%%%%%
%%%%%   #  Name  e-   N   O          Mw    Charge       Phase     Hf(298)
%%%%%                            [g/mol]                          [kJ/mol]
%%%%%Gas Species (11):
%%%%%   1: e-     1   0   0     0.00055        -1         gas        0.00
%%%%%   2: N      0   1   0     14.0067         0         gas      472.44
%%%%%   3: N+    -1   1   0     14.0061         1         gas     1881.90
%%%%%   4: O      0   0   1     15.9994         0         gas      249.23
%%%%%   5: O+    -1   0   1     15.9989         1         gas     1568.84
%%%%%   6: NO     0   1   1     30.0061         0         gas       91.09
%%%%%   7: N2     0   2   0     28.0134         0         gas        0.00
%%%%%   8: N2+   -1   2   0     28.0129         1         gas     1509.51
%%%%%   9: O2     0   0   2     31.9988         0         gas        0.00
%%%%%  10: O2+   -1   0   2     31.9982         1         gas     1171.41
%%%%%  11: NO+   -1   1   1     30.0055         1         gas      990.65


close all
clear
clc

pkg load optim; % Load optimization package

page_screen_output(0);

% ==========   PARAMETERS! MODIFY HERE!   ===============
strNumStart = 1;
strNumEnd   = 18;
strStep     = 1; % Direct ordering

% ==========   FLAGS   ==================================

PLOT_FLAG = 0;

% ==========   FREE STREAM VALUES    ====================

U_fs = 20000;
T_fs = 220;
rho_fs = 5.65e-5;

% ==========   INTERPOLATION ALONG STREAMLINE    ========

x_start_interp = 0.0255;
x_end_interp   = 0.34;
num_x_interp   = 30;

% ==========   EXTRAPOLATING STREAMLINE    ==============

x_start_fitting = 0.05; 
x_end_fitting   = 10;
x_step_fitting  = 0.05;

% ==========   NUMBER OF VERTICAL ADDITIONS    ==========

nClonesTop = 8;

% ==========   START CYCLING OVER STREAMLINES   =========

streamlines_plot = figure;

strCount = 0; % Intialize counter
for(strNumber = strNumStart:strStep:strNumEnd)
  fprintf('Reading streamline %d of %d\n', strNumber, strNumEnd);

  filename_in  = sprintf('exported_streamlines/streamline_%05d.dat', strNumber);

  strCount = strCount + 1; % Increment the counter

  dd = load(filename_in);
  
  %x    = dd(:, end-2);
  %y    = dd(:, end-1);
  %nrho = dd(:, 5);
  %T    = dd(:, 6);
  %u    = dd(:, 7);
  %v    = dd(:, 8);
  
  x    = dd(:, end-2);
  y    = dd(:, end-1);
  nrho = dd(:, 12);
  T    = dd(:, 9);
  u    = dd(:, 6);
  v    = dd(:, 7);

  ni   = dd(:, [13:24]);  
  % Sparta order: 
  % 1   2   3  4  5   6  7   8   9    10   11  12
  % O2  N2  O  N  NO  e  Mg  Np  N2p  NOp  Op  O2p
  %
  % e  N  N+  O   O+   NO   N2   N2+   O2   O2+   NO+
  ni_MPPorder = [ni(:, 6), ni(:, 4), ni(:, 8), ni(:, 3), ni(:, 11), ni(:, 5), ni(:, 2), ni(:, 9), ni(:, 1), ni(:, 12), ni(:, 10)];
  
  % Some properties
  Navo = 6.02214E+23;
  % e-  N   N+  O   O+  NO  N2  N2+ O2  O2+ NO+ 
  Mw_vect_MPPorder = [0.00055, 14.0067, 14.0061, 15.9994, 15.9989, 30.0061, 28.0134, 28.0129, 31.9988, 31.9982, 30.0055]/1000; % [kg/mol]

  % Computing mole fractions
  nSpecies = size(ni_MPPorder,2);
  n  = sum(ni_MPPorder,2);
  Xi_MPPorder = ni_MPPorder./repmat(n,1,nSpecies);

  % Computing rho
  rho = 1/Navo*sum(repmat(Mw_vect_MPPorder,size(ni_MPPorder,1),1).*ni_MPPorder, 2);
  
  % ====================================================================
  % ======    INTERPOLATION PHASE!!    =================================
  % ====================================================================
  % This phase might be required so that Larsen have as input some files
  % that do not differ in number of points! Thus it's more likely that 
  % results are unaffected by the change of number of cells
  %
  % ALTRIMENTI POSSO FARE UNO STUDIO DI SENSITIVITA' PER OGNI TESTCASE:
  % PRIMA LANCIO LA SIMULAZIONE CON IL NUMERO COMPLETO DI PUNTI SULLA
  % STAGNATION LINE, POI LO RIFACCIO CON LA META' DEI PUNTI E VEDO SE
  % TORNA!!! Se non torna, allora probabilmente la densita' di punti 
  % viene 
  
  xx = linspace(x_start_interp,x_end_interp, num_x_interp); % THINK ABOUT MAKING IT LINEARLY GROW!!!!!!
  x_interp   = interp1(x, x, xx);
  y_interp   = interp1(x, y, xx);
  T_interp   = interp1(x, T, xx);
  rho_interp = interp1(x, rho, xx);
  u_interp   = interp1(x, u, xx);
  v_interp   = interp1(x, v, xx);

  Xi_interp = [];
  for(s_id = 1:nSpecies)
    Xi_interp(:,s_id) = interp1(x, Xi_MPPorder(:,s_id), xx);
  end 
  
  % -----  Now fit  -----

  x_plot = [x_start_fitting:x_step_fitting:x_end_fitting];

  T_funct   = @(p, x) T_fs + p(1)*exp(p(2)*(x-x_start_fitting)); % Temp model function
  U_funct   = @(p, x) U_fs + p(1)*exp(p(2)*(x-x_start_fitting)); % Vel model function
  %rho_funct = @(p, x) rho_fs + p(1)*exp(p(2)*(x-x_start_fitting) ...
  %                           + p(3)*exp(p(4)*(x-x_start_fitting))); % Dens model function
  rho_funct = @(p, x) rho_fs + p(1)*exp(p(2)*(x-x_start_fitting)); % Dens model function

  pos_x = find(x_interp > x_start_fitting);
  pos_x = pos_x(1); % first

  % First extimate the coefficients, then use them for the real computation
  coefs_T = nonlin_curvefit(T_funct, [6000;-3], x_interp(pos_x:end), T_interp(pos_x:end));
  coefs_T = nonlin_curvefit(T_funct, coefs_T, x_interp(pos_x:end), T_interp(pos_x:end));

  coefs_U = nonlin_curvefit(U_funct, [-1000;-2], x_interp(pos_x:end), u_interp(pos_x:end));
  coefs_U = nonlin_curvefit(U_funct, coefs_U, x_interp(pos_x:end), u_interp(pos_x:end));

  coefs_rho = nonlin_curvefit(rho_funct, [0.3;20;-10;1], x_interp(pos_x:end), rho_interp(pos_x:end));
  %coefs_rho = nonlin_curvefit(rho_funct, [0.3;20], x_interp(pos_x:end), rho_interp(pos_x:end));
  coefs_rho = nonlin_curvefit(rho_funct, coefs_rho, x_interp(pos_x:end), rho_interp(pos_x:end));


  %Plot it
  if PLOT_FLAG == 1
    figExtend = figure;
    subplot(3,1,1)
    plot(x_interp, u_interp, 'k', 'linewidth', 2)
    hold on
    plot(x_plot, U_funct(coefs_U,x_plot), 'r', 'linewidth', 2)
    title(['y = ', num2str(y_interp(1))]);
    xlabel('x [m]')
    ylabel('u [m/s]')
  
    subplot(3,1,2)
    plot(x_interp, rho_interp, 'k', 'linewidth', 2)
    hold on
    plot(x_plot, rho_funct(coefs_rho,x_plot), 'r', 'linewidth', 2)
    xlabel('x [m]')
    ylabel('rho [kg/m3]')
  
    subplot(3,1,3)
    plot(x_interp, T_interp, 'k', 'linewidth', 2)
    hold on
    plot(x_plot, T_funct(coefs_T,x_plot), 'r', 'linewidth', 2)
    xlabel('x [m]')
    ylabel('T [K]')
  
    fprintf('Press enter to continue\n')
    pause()
    close(figExtend);
  end

  % Extending the streamline
  x_extrap   = [x_interp(end)+x_step_fitting : x_step_fitting : x_end_fitting];
  y_extrap   = y_interp(end)*ones(size(x_extrap));
  u_extrap   = U_funct(coefs_U, x_extrap);
  v_extrap   = 0.0*u_extrap;
  T_extrap   = T_funct(coefs_T, x_extrap);
  rho_extrap = rho_funct(coefs_rho, x_extrap);

  Xi_extrap = [];
  for(s_id = 1:nSpecies)
    Xi_extrap(:,s_id) = Xi_interp(end,s_id)*ones(size(x_extrap));
  end

  % And merge the old and the new part

  x_extended   = [x_interp, x_extrap];
  y_extended   = [y_interp, y_extrap];
  u_extended   = [u_interp, u_extrap];
  v_extended   = [v_interp, v_extrap];
  T_extended   = [T_interp, T_extrap];
  rho_extended = [rho_interp, rho_extrap];

  Xi_extended = [];
  for(s_id = 1:nSpecies)
    Xi_extended(:,s_id) = [Xi_interp(:,s_id); Xi_extrap(:,s_id)];
  end

  x_all(:,strNumber)   = x_extended';
  y_all(:,strNumber)   = y_extended';
  u_all(:,strNumber)   = u_extended';
  v_all(:,strNumber)   = v_extended';
  rho_all(:,strNumber) = rho_extended';
  T_all(:,strNumber)   = T_extended';

  for(s_id = 1:nSpecies)
    Xi_all(:,s_id, strNumber) = Xi_extended(:,s_id);
  end

  figure(streamlines_plot);
  hold on
  plot(x_interp, y_interp, 'k', 'linewidth', 2);
  plot(x_extrap, y_extrap, 'r', 'linewidth', 2);
  axis equal
  grid on
  pause(0.01);
 
end 

% ==============================================================================
% ==============================================================================
% ==============================================================================
% ==============================================================================
% ==============================================================================
% 
% Now the script creates a number nClonesTop of streamlines and puts them to the top.
% Such streamlines have the same shape as the last one and copy everything except
% the composition, that is set to freestream N2, O2.
%

% FREESTREAM VALUES
X_fs_N2 = 0.79;
X_fs_O2 = 0.21;
X_fs_vect = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, X_fs_N2, 0.0, X_fs_O2, 0.0, 0.0];

% Parabolic spacing
deltay_0 = 0.002;
a = 0.0003;
b = 0.00005;

deltay = @(id) a*id.^2 + b*id + deltay_0;

for(ii = 1:nClonesTop)
  x_all(:, strNumEnd + ii)   = x_all(:,end);
  y_all(:, strNumEnd + ii)   = y_all(:,end) + deltay(ii-1);
  u_all(:, strNumEnd + ii)   = u_all(:,end); %% ???????????????????????????????????
  v_all(:, strNumEnd + ii)   = v_all(:,end);
  T_all(:, strNumEnd + ii)   = T_all(:,end); %%%%%%%%%%%%% HERE PICK FROM DSMC DATA
  rho_all(:, strNumEnd + ii) = rho_all(:,end); %%%%%%%%%%% HERE PICK FROM DSMC DATA

  for(s_id = 1:nSpecies)
    Xi_all(:, s_id, strNumEnd+ii) = X_fs_vect(s_id)*ones(size(x_all,1),1);
  end

 
  plot(x_all(:, strNumEnd + ii), y_all(:, strNumEnd + ii), 'g', 'linewidth', 2)

end

% #########################   Write streamlines   ##############################

strCount = 0;
for str_id = 1:size(x_all,2)


  filename_out = sprintf('output_streamlines/streamline_%05d', strCount);

  strCount = strCount + 1; % Increment counter

  fid = fopen(filename_out, 'w');

  for ii = 1: size(x_all,1) % for each point along the streamline

    fprintf(fid, '%e %e %e %e %e %e ', x_all(ii,str_id), y_all(ii, str_id), T_all(ii, str_id), rho_all(ii, str_id), u_all(ii,str_id), v_all(ii,str_id));

    % Air 11:  e-  N   N+  O   O+  NO  N2  N2+ O2  O2+ NO+ 
    for(s_id = 1:nSpecies) % Write species
      fprintf(fid, '%e ', Xi_all(ii, s_id, str_id));
    end
    fprintf(fid, '\n');

  end
  
  fclose(fid);

  fprintf('Streamline %d of %d written.\n', strCount, size(x_all,2));
end

