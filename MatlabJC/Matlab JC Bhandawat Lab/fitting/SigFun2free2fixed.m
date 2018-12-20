function fit = SigFun2free2fixed(coef,Input) 

% a sigmoid from cumulative gaussian
% Input should be structure with fields
%  coef1=mean, coef2=var

x = Input.x ; % x for data points to fit
MaxSat = Input.MaxSat ; % amplitude of cumulative gauss
MinSat = Input.MinSat ; % min amplitude of cumulative gauss

fit = (MaxSat-MinSat)*(cumsum(exp(-((x-coef(1)).^2)/(2*coef(2))))/sum(exp(-((x-coef(1)).^2)/(2*coef(2)))))+MinSat ;