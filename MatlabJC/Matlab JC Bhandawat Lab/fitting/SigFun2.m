function fit = SigFun2(coef,x) 

% sigmoid fit from cumulative gaussian
% JC 8/8/2012
% coef1 = amp, coef2=mean, coef3=variance, coef4 = offset

fit = coef(1)*cdf('norm',x,coef(2),coef(3)) + coef(4) ;