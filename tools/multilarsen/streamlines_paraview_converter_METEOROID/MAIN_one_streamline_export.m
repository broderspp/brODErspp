close all
clear
clc

% ===============================  LOAD PARAVIEW OUTPUT  ================================
% Load streamlines file
strData = load('venti_streamlines.csv');

xAll = strData(:,end-2);
yAll = strData(:,end-1);
zAll = strData(:,end);

% ====================================  PROCESSING  =====================================
% Domain limits (in the HP that the low goes from left to right)
xMin = min(xAll);
xMax = max(xAll);

LxDomain = xMax - xMin;

DiffLim  = LxDomain / 6; % In the hope that any streamline is longer than 3 points..

% Find when streamlines split, that is when the difference of two points is almost
% big as the x domain
diffBoolVect = abs(diff(xAll)) > DiffLim;

numStreamlines = nnz(diffBoolVect) + 1; % The number of streamlines is the number of
                                        % jumps plus 1.
startStreamlines = find(diffBoolVect == 1) + 1; % Starting point of each streamline
streamlinesExtremes = [1; startStreamlines; numel(xAll)]; % Add the last point.

%%% SOME BULLSHIT HERE
%%%  % Now, it's possible that paraview gave me one more jump after the last streamline.. if it's
%%%  % the case then remove it! And also remove the last idStreamlines element
%%%  if idStreamlines(end) == numel(xAll)
%%%    numStreamlines = numStreamlines - 1;
%%%    idStreamlines  = idStreamlines(1:end-1);
%%%  end

% Save streamlines into a struct of cell arrays
streamlines.x = {};
streamlines.y = {};
streamlines.z = {};

streamlinesplot = figure;
hold on

for strID = 1:numStreamlines

  streamlines.x{strID} = xAll(streamlinesExtremes(strID) : streamlinesExtremes(strID+1)-1);
  streamlines.y{strID} = yAll(streamlinesExtremes(strID) : streamlinesExtremes(strID+1)-1);
  streamlines.z{strID} = zAll(streamlinesExtremes(strID) : streamlinesExtremes(strID+1)-1);

  plot(streamlines.x{strID}, streamlines.y{strID},'-k', 'linewidth', 2)
end

% SUPERIMPOSING CYLINDER
R_cyl = 0.152;
Nth = 100;
thvec = linspace(0,2*pi, Nth + 1);
thvec = thvec(1:end-1);
xvec = R_cyl * cos(thvec);
yvec = R_cyl * sin(thvec);

figure(streamlinesplot)
fill(xvec, yvec, 'b')
%ylim([0, 0.5])
%xlim([-0.5, 0])
ylim([0, 0.95])
xlim([0, 1])

axis equal

% ===============================  EXPORT STREAMLINES  ====================================

% Which streamline do you wish to export?
fprintf('I have found %d streamlines.\n',numStreamlines);
user_choice = input('Exporting a streamline? If yes type the number, otherwise type "N" ','s');

if user_choice ~= 'N' % Then it's hopefully a number
  num_str = str2num(user_choice)
  if num_str > numStreamlines
    fprintf('\nMORON! You exceeded the available streamlines!\n\n');
    error('ABORTING.')
  end

  % Over-plotting the streamline
  figure(streamlinesplot)
  plot( streamlines.x{num_str}, streamlines.y{num_str}, 'g', 'linewidth',2 )

else
  fprintf('Grazie di aver volato con noi.\n');
  return
end


% Exporting that streamline
fid = fopen('user_choosed_streamline.dat','w');

for ii = streamlinesExtremes(num_str): streamlinesExtremes(num_str+1)-1
  DataNow = strData(ii, :); % Extract all the data on this streamline point
  for jj = 1:numel(DataNow)
    fprintf(fid, '%e ', DataNow(jj));
  end
  fprintf(fid, '\n');
end

fclose(fid);

