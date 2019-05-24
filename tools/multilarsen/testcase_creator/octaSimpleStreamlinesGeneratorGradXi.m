close all
clear
clc

% Parameters 
xMin = 0.0; % [m]
xMax = 100.0; % [m]

yMin = 0.0; % [m]
yMax = 1.0; % [m]

U   = 1000;  % mixture horz velocity [m/s]
V   = 0;     % mixture vert velocity [m/s]
rho = 0.001; % mixture density [kg/m3]
T   = 1000;  % mixture temperature [K]

Nelements    = 30;
Nstreamlines = 10;

x_vect = linspace(xMin, xMax, Nelements);
y_vect = linspace(yMin, yMax, Nstreamlines);

% Create and write streamlines
for(str_id = 1:Nstreamlines)

  % Create a streamline file
  filename = sprintf('./streamlines/str_%05d', str_id-1);
  filename 
  fid = fopen(filename,'w');

  % March along streamline
  for(point_id = 1:Nelements)
    x = x_vect(point_id);
    y = y_vect(str_id);

    if str_id == 5
      X_N = 0.9;
      X_O = 0.1;
    else
      X_N = 1.0;
      X_O = 0.0;
    end

    % x, y, T, rho, U, V, Xi
    fprintf(fid, '%e %e %e %e %e %e %e %e %e %e %e\n',x,y,T,rho,U,V,X_N,X_O,0,0,0);
  end

  % Close streamline file
  fclose(fid);
end
