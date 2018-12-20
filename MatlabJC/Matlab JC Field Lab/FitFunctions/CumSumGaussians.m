function [fit] = CumSumGaussians(coef, x) 

% the cumulative sum of the difference between two gaussians (center and surround)
% JC 12/2/2015

% coef(1) = Center amp
% coef(2) = Center std
% coef(3) = Surround amp
% coef(4) = Surround std

longPoint = 10000 ;

centerTemp = (normcdf(x,0,coef(2))-normcdf(0,0,coef(2))) ; % cumulative gaussian
center = coef(1)*centerTemp/(normcdf(longPoint,0,coef(2))-normcdf(0,0,coef(2))) ; % scaled

surroundTemp = (normcdf(x,0,coef(4))-normcdf(0,0,coef(4))) ;
surround = -coef(3)*surroundTemp/(normcdf(longPoint,0,coef(4))-normcdf(0,0,coef(4))) ;

fit = center+surround ;

