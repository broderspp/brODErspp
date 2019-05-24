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

page_screen_output(0);

% ==========   PARAMETERS! MODIFY HERE!   ===============
strNumStart = 1;
strNumEnd   = 18;
%strStep     = -1; % Reverse order: from streamline 20 to 4
strStep     = 1; % Direct ordering

% ==========   START CYCLING OVER STREAMLINES   =========
strCount = 0; % Intialize counter
for(strNumber = strNumStart:strStep:strNumEnd)
  fprintf('Processing streamline %d\n', strNumber);
  fprintf('Printing streamline %d\n\n', strCount);

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
  Navo = 6.02214E23;
  % e-  N   N+  O   O+  NO  N2  N2+ O2  O2+ NO+ 
  Mw_vect_MPPorder = [0.00055, 14.0067, 14.0061, 15.9994, 15.9989, 30.0061, 28.0134, 28.0129, 31.9988, 31.9982, 30.0055];

  % Computing mole fractions
  nSpecies = size(ni_MPPorder,2);
  n  = sum(ni_MPPorder,2);
  Xi_MPPorder = ni_MPPorder./repmat(n,1,nSpecies);

  % Computing rho
  rho = 1/Navo*sum(repmat(Mw_vect_MPPorder,size(ni_MPPorder,1),1).*ni_MPPorder, 2);
  
%% NOT NEEDED %%  % Computing module of velocity
%% NOT NEEDED %%  modU = sqrt(u.^2 + v.^2); % Module of velocity
  
%% NOT NEEDED %%  % Curvilinear abscissa
%% NOT NEEDED %%  ds = sqrt(diff(x).^2 + diff(y).^2); % Length of small streamline element
%% NOT NEEDED %%  s = [0; cumsum(ds)];
  
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
  
  xx = linspace(0.025, 0.34, 30); % THINK ABOUT MAKING IT LINEARLY GROW!!!!!!
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
  
  for(ii = 1: numel(x_interp))

    % Air 11:  e-  N   N+  O   O+  NO  N2  N2+ O2  O2+ NO+ 
    fprintf(fid, '%e %e %e %e %e %e ', x_interp(ii), y_interp(ii), T_interp(ii), rho_interp(ii), u_interp(ii), v_interp(ii));
    for(s_id = 1:nSpecies) % Write species
      fprintf(fid, '%e ', Xi_interp(ii, s_id));
    end
    fprintf(fid, '\n');

    %% Air 5: N, O, NO, N2, O2
    %fprintf(fid, '%e %e %e %e %e %e %e %e %e %e %e\n', x_interp(ii), y_interp(ii), T_interp(ii), rho_interp(ii), u_interp(ii), v_interp(ii), 1.0, 0.0, 0.0, 0.0, 0.0);

    %% Argon, 3 species: el, Ar+, Ar
    %fprintf(fid, '%e %e %e %e %e %e %e %e %e\n', x_interp(ii), y_interp(ii), T_interp(ii), rho_interp(ii), u_interp(ii), v_interp(ii), 0.0, 0.0, Yi_interp(ii));

%%%%%    %fprintf(fid, '%e %e %e %e %e %e %e %e %e %e %e\n', x(ii), y(ii), T(ii), rho(ii), u(ii), v(ii), Yi(ii), 0.0, 0.0, 0.0, 0.0);
%%%%%  %  fprintf(fid, '%e %e %e %e %e %e\n', s(ii), T(ii), 0.0, rho(ii), modU(ii), Yi(ii));
  end
  
  fclose(fid);
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

nClonesTop = 10; % Number of streamlines to be cloned above

dy = MUORIIIII

for(ii = 1:nClonesTop)
  % Create new 
end
