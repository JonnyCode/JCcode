function fit = ThresholdingNL(coef,x)

% thresholding nonlinearity on log-log axis 
% weber like funtion when ceof(2)=1 ; rose-devries when s=2 see Angueyra and Rieke 2013)
% coef(1) sets threshold at which curve begins to become linear
% coef(2) sets slope of linear portion
% coef(3) sets offset of cure

fit = ((1+x/coef(1)).^(1/coef(2)))*coef(3) ;


