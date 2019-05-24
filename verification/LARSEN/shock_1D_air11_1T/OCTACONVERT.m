close all
clear
clc

% Load the data
dd = load('output.dat');

% extract
dd_new = dd(:,[1,2,4:end]);

% Print
fid = fopen('outputNEW.dat', 'w');

for ii = 1:size(dd_new,1)
  for jj = 1:size(dd_new,2)
    fprintf(fid, '%.10e   ', dd_new(ii,jj))
  end
  fprintf(fid, '\n')
end

% Close file
fclose(fid);
