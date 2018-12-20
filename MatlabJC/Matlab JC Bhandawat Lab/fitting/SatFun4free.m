function fit = SatFun4free(coef,x) 

% adapted from "SatFun"
% added two more free variables

% coef1 - max point
% coef2 - exponent
% coef3 - mid point
% coef4 - offset

fit = (coef(1)-coef(4))*((x.^coef(2))./(x.^coef(2) +coef(3).^coef(2)))+coef(4) ;


