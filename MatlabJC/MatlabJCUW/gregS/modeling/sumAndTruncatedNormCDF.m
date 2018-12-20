function y = sumAndTruncatedNormCDF(x,mu,sigma,A,C,T)
%T is transform (not used)
x = sum(x);
y = (A*normcdf(x,mu,sigma) + C) / (A*normcdf(1,mu,sigma) + C) ;
y(y<0) = 0;
%y = y./max(y);
%y(y>1) = 1;
