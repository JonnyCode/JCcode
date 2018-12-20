function angle_degrees = convertRad2Deg(angle_rad) 

% converts radians to degrees and fixes negatives
% JC 10/28/16

angle_degrees = angle_rad*180/pi ;

angle_degrees(angle_degrees<0) = angle_degrees(angle_degrees<0)+360 ;