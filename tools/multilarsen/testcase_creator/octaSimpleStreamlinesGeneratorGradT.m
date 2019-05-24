close all
clear
clc

page_screen_output(0);

% Parameters 
xMin = 0.0; % [m]
xMax = 100.0; % [m]

yMin = 0.0; % [m]
yMax = 1.0; % [m]

U   = 1000;  % mixture horz velocity [m/s]
V   = 0;     % mixture vert velocity [m/s]
rho = 0.001; % mixture density [kg/m3]
T   = 5000;  % mixture temperature [K]

Nelements    = 30;
Nstreamlines = 30;

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

    if (str_id == 4) || (str_id == 5) || (str_id == 6)
      T_now = T + 3000;
    else
      T_now = T;
    end

    % x, y, T, rho, U, V, Xi
    fprintf(fid, '%e %e %e %e %e %e %e %e %e %e %e\n',x,y,T_now,rho,U,V,0.,0.,0.,0.79,0.21);
  end

  % Close streamline file
  fclose(fid);
end
