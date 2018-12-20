function fit = SigFunInterp2free2fixed(coef,Input) 

% a sigmoid from cumulative gaussian
% fiting only with commented code does not allow acurate interpolation later (not clear exactly how to fix this properly)
%  coef1=mean, coef2=var
x = Input.x ; % x for data points to fit
MaxSat = Input.MaxSat ; % amplitude of cumulative gauss
MinSat = Input.MinSat ; % offset of cumulative gauss 
Xinterp = Input.Xinterp ; % the interpolated x you want to fit data to (must intersect all x x points)

[v,Xpoints] = intersect(Xinterp,x) ;

fitInterp = (MaxSat-MinSat)*cumsum(exp(-((Xinterp-coef(1)).^2)/(2*coef(2))))/sum(exp(-((Xinterp-coef(1)).^2)/(2*coef(2))))+MinSat ;
fit = fitInterp(Xpoints) ;

