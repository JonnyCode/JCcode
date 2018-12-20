function fit = SigFun2free(coef,input) 

% a sigmoid from cumulative gaussian
% input.x = x values
% input.amp = amplitude of sigmoid
% coef1=(amp-offset), coef2=mean,
x = input.x ;

%fit = input.amp*cumsum(exp(-((x-coef(1)).^2)/(2*coef(2))))/sum(exp(-((x-coef(1)).^2)/(2*coef(2)))) ;
fit = input.amp*normcdf(x,[coef(1),coef(2)]) ;