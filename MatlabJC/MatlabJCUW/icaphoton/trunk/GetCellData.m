function mymatrix = GetCellData(mycells,CellInfo)
%  mymatrix = GetCellData(mycells, allcells)
%  
%   Given an array of the cells you want to pull out, it pulls that data out of
%   CellInfo, and converts it to a matrix.  All cells must be same length!!!
%   If you feed this routine a cell array with at least one cell containing
%   anything but a numeric scalar or an empty array, this routine will barf.
%
%   Oh and please don't feed it cell arrays with more than two dimensions.

% RCS info: $Id$
%
% created by:
% jpg 5/20/98 - out of *extreme* frustration and irritation
% adapted to cellinfo structure by:
% mkm 4/01

data = CellInfo.EpochData.Data(mycells);
xsize = length(data{1});
ysize = length(mycells);
padding = zeros(xsize,1);
mymatrix = zeros(xsize,ysize);  
for y = 1:ysize
	if isempty(data{y}), mymatrix(:,y)=padding;
	else
		mymatrix(:,y)=data{y}';
	end
end

