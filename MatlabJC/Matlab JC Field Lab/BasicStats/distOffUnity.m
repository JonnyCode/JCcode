function DoU = distOffUnity(x,y) 

% this function will find the distance of a point off of unity
% x, y are cooridinates

% JC 12/29/2015

DoU = sqrt(2*(x-y)^2)/2 ; % JC math

%DoU = abs(y*x-y^2)/sqrt(2*y^2) ; % this seemed equivilent, found on
%wikipedia - distance of point to line