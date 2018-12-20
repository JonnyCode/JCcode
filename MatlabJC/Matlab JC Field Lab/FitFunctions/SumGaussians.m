function  [fit] = SumGaussians(coef, x) 

% sum of two gaussians (see Demb et al 2001, Enroth-Cugell and Roabson 1966)
% JC 2/19/2016

% coef(1) = Center amp
% coef(2) = Center std
% coef(3) = Surround amp
% coef(4) = Surround std

fit = (coef(1)*(1-exp(-x/coef(2)).^2))-(coef(3)*(1-exp(-x/coef(4)).^2)) ;