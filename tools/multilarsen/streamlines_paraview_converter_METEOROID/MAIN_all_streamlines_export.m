close all
clear
clc

% ===============================  LOAD PARAVIEW OUTPUT  ================================
% Load streamlines file
%strData = load('/home/starlight/tmp/met_20kms/some_streamlines_20kms.csv');
strData = load('/mydir/VKI/RM/PROJ/meteor_trail_final_simulations/meteor_20kms_1cm/initial_data/trail_20kms_18_streamlines.csv');

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

figure(streamlinesplot)
axis equal

% ===============================  EXPORT STREAMLINES  ====================================

% Which streamline do you wish to export?
fprintf('I have found %d streamlines.\n',numStreamlines);
user_choice = input('Exporting streamlines? [Y/N] ','s');

if user_choice == 'Y' % Then it's hopefully a number

  % Exporting that streamline
  for str_id = 1:numel(streamlines.x)
  %for str_id = 4:numel(streamlines.x)

    str_name = sprintf('exported_streamlines/streamline_%05d.dat', str_id);
    fid = fopen(str_name, 'w');
    
    for ii = streamlinesExtremes(str_id): streamlinesExtremes(str_id+1)-1
      DataNow = strData(ii, :); % Extract all the data on this streamline point
      for jj = 1:numel(DataNow)
        if strData(ii, end-2) < 0.85
          fprintf(fid, '%e ', DataNow(jj));
        end
      end
      fprintf(fid, '\n');
    end

    fclose(fid);

    % Over-plotting the streamline
    figure(streamlinesplot)
    hold on
    plot( streamlines.x{str_id}, streamlines.y{str_id}, 'g', 'linewidth',2 )

  end

else
  fprintf('Grazie di aver volato con noi.\n');
  return
end



