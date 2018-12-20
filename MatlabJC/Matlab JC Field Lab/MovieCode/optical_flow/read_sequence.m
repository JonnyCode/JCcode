% Function to read an image sequence
% Supported: Fleet's binary format and SUN rasterfiles
%
% Usage: II = read_sequence (stem, middle_frame, st, [sx, sy]);
%	stem	common stem of filenames, e.g. "yosemite/yosemite256."
%  middle_frame number of middle frame [10]
%	st	temporal span [5]
%	sx,sy	spatial dimensions, only provided for the binary format

function II = read_sequence (stem1, mid, st, sx, sy)

if (nargin<1)
	stem1 = 'yosemite/yosemite256.';
end
if (nargin<2)
	mid = 10;
end
if (nargin<3)
	st = 5;
end
if (nargin<4)
	bin_flag = 0;
else
	bin_flag = 1;
end

aux1 = floor(st/2);

if (rem(st,2)==0)
	IND = mid-aux1-1:mid+aux1;
else
	IND = mid-aux1:mid+aux1;
end

if (bin_flag)
	II = zeros (sy,sx,st);
	for t=1:st
		fid = fopen ([stem1 num2str(IND(t))], 'r');

		% THERE IS NO HEADER IN TREE DATA
		X = fread(fid,[sx,sy],'uchar')';
		fclose(fid);

		II(:,:,t) = X;
	end
else
	% Load First Frame
	name_in = [stem1 num2str(IND(1))];
	X = rastread (name_in);
	[sy sx] = size(X);
	II = zeros (sy,sx,st);
	II(:,:,1) = X;

	% Load Rest
	for t=1:st
		name_in = [stem1 num2str(IND(t))];
		II(:,:,t) = rastread (name_in);
	end
end


