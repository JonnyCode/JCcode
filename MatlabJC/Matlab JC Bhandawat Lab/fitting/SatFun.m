function fit = SatFun(coef,x) ;

% a saturating powerlaw function 
expon = 1.5 ;
fit = coef(1)*(x.^expon)./(x.^expon +coef(2)) ;