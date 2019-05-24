close all
clear
clc

page_screen_output(0);

% ==========   PARAMETERS! MODIFY HERE!   ===============
strNumStart = 20;
strNumEnd   = 4;
strStep     = -1; % Reverse order: from streamline 20 to 4

% ==========   START CYCLING OVER STREAMLINES   =========
strCount = 0; % Intialize counter
for(strNumber = strNumStart:strStep:strNumEnd)
  fprintf('Processing streamline %d\n', strNumber);
  fprintf('Printing streamline %d\n', strCount);

  filename_in  = sprintf('exported_streamlines/streamline_%05d.dat', strNumber);
  filename_out = sprintf('output_streamlines/streamline_%05d', strCount);

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
  nrho = dd(:, 2);
  T    = dd(:, 3);
  u    = dd(:, 6);
  v    = dd(:, 7);
  
  % Some properties
  Navo = 6.02214E23;
  Mmol = 39.948;
  
  % Computing rho
  rho  = nrho.*Mmol/Navo;
  
  % Computing Yi (it's trivial... rhoi / rho = 1 in this case!!)
  Yi   = ones(size(rho));
  
  % Computing module of velocity
  modU = sqrt(u.^2 + v.^2); % Module of velocity
  
  % Curvilinear abscissa
  ds = sqrt(diff(x).^2 + diff(y).^2); % Length of small streamline element
  s = [0; cumsum(ds)];
  
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
  
  
  % ====================================================================
  % =======   SAVING RESULTS    ========================================
  % ====================================================================
  % Saving results in file output_SPARTA.dat, in the format that Larsen
  % wishes to eat.
  
  %fid = fopen('output_SPARTA.dat','w');
  fid = fopen(filename_out, 'w');
  
  xx = linspace(0.25, 0.8, 30);
  x_interp   = interp1(x, x, xx);
  y_interp   = interp1(x, y, xx);
  T_interp   = interp1(x, T, xx);
  rho_interp = interp1(x, rho, xx);
  u_interp   = interp1(x, u, xx);
  v_interp   = interp1(x, v, xx);
  Yi_interp  = interp1(x, Yi, xx);
  
  for(ii = 1: numel(x_interp))
    % Air 5: N, O, NO, N2, O2
    fprintf(fid, '%e %e %e %e %e %e %e %e %e %e %e\n', x_interp(ii), y_interp(ii), T_interp(ii), rho_interp(ii), u_interp(ii), v_interp(ii), 1.0, 0.0, 0.0, 0.0, 0.0);

    %% Argon, 3 species: el, Ar+, Ar
    %fprintf(fid, '%e %e %e %e %e %e %e %e %e\n', x_interp(ii), y_interp(ii), T_interp(ii), rho_interp(ii), u_interp(ii), v_interp(ii), 0.0, 0.0, Yi_interp(ii));

%%%%%    %fprintf(fid, '%e %e %e %e %e %e %e %e %e %e %e\n', x(ii), y(ii), T(ii), rho(ii), u(ii), v(ii), Yi(ii), 0.0, 0.0, 0.0, 0.0);
%%%%%  %  fprintf(fid, '%e %e %e %e %e %e\n', s(ii), T(ii), 0.0, rho(ii), modU(ii), Yi(ii));
  end
  
  fclose(fid);
end 
