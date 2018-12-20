function PV = PolarVectorAddition(Vectors) 

%this function will take a set of polar vectors in a matrix, Vectors, and add
%them. Vectors has a row for every vector, first column is angle in degrees,
%second column is magnitude.  The output, PV is [degrees, magnitude]. 360
%is not an output, only 0 degrees.

% JC 9/10/15

for v=1:size(Vectors,1) ; % for each vector
    x(v) = cosd(Vectors(v,1))*Vectors(v,2) ;
    y(v) = sind(Vectors(v,1))*Vectors(v,2) ;
end

x_sum = sum(x) ;
y_sum = sum(y) ;

A = atan2d(y_sum, x_sum) ;
if A>=0 ;
    PV(1) = A ;
else
    PV(1) = 360 + A ;
end

PV(2) = sqrt(y_sum^2 + x_sum^2) ;
