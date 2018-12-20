function [X,map]=rastread(filename);
%RASTREAD Read a BMP SUN raster file
%  [X,MAP]=RASTREAD('filename') reads the file 
%  'filename' and returns the indexed image X and
%  associated colormap MAP. If no extension is given 
%  for the filename, the extension '.rast' is assumed.
%
%  See also: BMPREAD, BMPWRITE, GIFREAD, HDFREAD, 
%            PCXREAD, TIFFREAD, XWDREAD.
%  Jeffery J. Faneuff  March 29, 1994
%  Copyright (c) 1994 by The MathWorks, Inc.
if nargin~=1 | isstr(filename)~=1
  error('Requires a string filename as an argument.');
end;
if (isempty(findstr(filename,'.'))==1)
  filename=[filename,'.ras'];
end;
fid=fopen(filename,'rb','b');		% Wuzz a bug: it read 'l'
if (fid==-1)
  disp('The file should have an extension of .ras ');
  error(['Error opening ',filename,'.']);
end;
magicnum=fread(fid,1,'uint32');
if (magicnum~= hex2dec('59a66a95'))
  fclose(fid);
  error('Not a SUN raster file.');
end;
width = fread(fid,1,'uint32');
height = fread(fid,1,'uint32');
depth = fread(fid,1,'uint32');
length = fread(fid,1,'uint32');
type = fread(fid,1,'uint32');
maptype = fread(fid,1,'uint32');
maplength = fread(fid,1,'uint32');
if (type~=1)
  fclose(fid);
  error('Not a standard encoded SUN raster file. ');
end;
if (maptype==1)
  rawmap = fread(fid,maplength,'uchar');
  map = reshape(rawmap,maplength/3,3);    % make a Nx3
                                          % colormap
  map = map./255;                         % normalize 
                                          % to 255
end
X = fread(fid,[width, height],'uchar')';
X = X(height:-1:1,:);
fclose(fid);


