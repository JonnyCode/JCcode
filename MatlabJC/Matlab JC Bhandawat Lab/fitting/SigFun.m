function fit = SigFun(coef,Input) 

% a sigmoid from cumulative gaussian
% fiting only with commented code does not allow acurate interpolation later (not clear exactly how to fix this properly)
%  coef1=mean, coef2=var
x = Input.x ; % x for data points to fit
Amp = Input.amp ; % amplitude of cumulative gauss
Xinterp = Input.Xinterp ; % the interpolated x you want to fit data to (must intersect all x x points)

[v,Xpoints] = intersect(Xinterp,x) ;

fitInterp = Amp*cumsum(exp(-((Xinterp-coef(1)).^2)/(2*coef(2))))/sum(exp(-((Xinterp-coef(1)).^2)/(2*coef(2)))) ;
fit = fitInterp(Xpoints) ;

%fit = coef(1)*cumsum(exp(-((x-coef(2)).^2)/(2*coef(3))))/sum(exp(-((x-coef(2)).^2)/(2*coef(3)))) ;
