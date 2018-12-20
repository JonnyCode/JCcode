% looking at trig functions


% two points in space
pA = [0,0] ; 
pB = [-5,5] ;

% some direction in space
nd = 300 ; 

% angle of two points
cA = atan2d(pB(1),pB(2)) ;

% acute angle between points vector and nd vector
dA = atan2d(sind(cA-nd),cosd(cA-nd)) ;

