function A = acuteAngle(A1,A2) 

% this function will find the acute angle between two values in degrees
% JC 11/19/2015

A = abs(atan2d(sind(A1 - A2), cosd(A1 - A2))) ;