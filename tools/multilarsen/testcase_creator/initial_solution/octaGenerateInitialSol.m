close all
clear
clc

page_screen_output(0)

% ==============  FREESTREAM CONDITION  ==================
T_fs   = 220;      % [K]
U_fs   = 20000;    % [m/s]
V_fs   = 0;        % [m/s]
rho_fs = 5.65e-5;  % [kg/m3]

% ==============  DOMAIN PARAMETERS  =====================
Nelements    = 315;
Nstreamlines = 55;

% Parabolic spacing
deltax_0 = 0.08;
a_x = 0.1;
b_x = 0.005;

deltay_0 = 0.003;
a_y = 0.001;
b_y = 0.0001;

deltax = @(id) a_x*id.^2 + b_x*id + deltax_0;
deltay = @(id) a_y*id.^2 + b_y*id + deltay_0;

x_vect = [0];
y_vect = [0];

for(ele_id = 1:Nelements)
  x_vect(ele_id + 1) = x_vect(ele_id) + deltax(ele_id-1);
end

for(str_id = 1:Nstreamlines)
  y_vect(str_id + 1) = y_vect(str_id) + deltay(str_id-1);
end


% PLOT the generated streamlines
[XX,YY] = meshgrid(x_vect, y_vect);

figure;
hold on
for(str_id = 1:Nstreamlines)
  plot(x_vect, y_vect(str_id)*ones(size(x_vect)), '-ok', 'linewidth', 2)
end

% =======  Load data file storing initial condition  ======
dd = load('FINAL_vertExt.csv');

% Extract it
yy_sim = dd(:,end-1);
TT_sim = dd(:,1);
ni_sim = dd(:,[2:12]);

% CORRECT: NO O2+ and its electrons
ions = [3, 5, 8, 11];
ni_sim(:,1) = sum(ni_sim(:, ions),2);
ni_sim(:,10) = zeros(size(ni_sim(:,11)));

n_sim = sum(ni_sim, 2);

Xi_sim = ni_sim./repmat(n_sim, 1, size(ni_sim,2));

yy_sim = [0; yy_sim];
TT_sim = [TT_sim(1); TT_sim];
Xi_sim = [Xi_sim(1,:); Xi_sim];

% Plot solution
initTemp = figure;
plot(yy_sim, TT_sim)

figure
plot(yy_sim, Xi_sim, 'k')

%%%% % Extended solution
%%%% yy_ext = [yy; 100];
%%%% TT_ext = [TT; T_fs];
%%%% 
%%%% for s_id = 1:size(Xi,2)
%%%%   Xi_ext(:,s_id) = []
%%%% end

% =======  Create and write streamlines   =========

for(str_id = 1:Nstreamlines)

  % Create a streamline file
  filename = sprintf('./streamlines/str_%05d', str_id-1);
  filename 
  fid = fopen(filename,'w');

  fprintf('Writing streamline %d\n', str_id);

  % March along streamline
  for(point_id = 1:Nelements)
    x = x_vect(point_id);
    y = y_vect(str_id);

    Xi = zeros(1,11);

    if(y < yy_sim(end))  % if I'm inside the previously simulated region
      T_now = interp1(yy_sim, TT_sim, y);
  
      % Mole fractions
      for s_id = 1:11
        Xi(s_id) = interp1(yy_sim, Xi_sim(:,s_id), y);
      end

    else % Freestream
      T_now = T_fs;

      Xi(7) = 0.79;
      Xi(9) = 0.21;
    end

    % x, y, T, rho, U, V, Xi
    fprintf(fid, '%e %e %e %e %e %e ', x, y, T_now, rho_fs, U_fs, V_fs);
    
    for s_id = 1:11
      fprintf(fid, '%e ', Xi(s_id))
    end
    fprintf(fid, '\n')

  end

  fclose(fid);

  T_vect(str_id) = T_now;
end


% Superimpose the numerical initial temperature profile
figure(initTemp)
hold on
plot(y_vect(1:end-1), T_vect, 'xr');
xlim([yy_sim(1), yy_sim(end)])
legend('Initial condition', 'Simulated points')
