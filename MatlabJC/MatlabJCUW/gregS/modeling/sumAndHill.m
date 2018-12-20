function y = sumAndHill(x,Kd,n,T)
%T is transform (not used)
x = sum(x);
%pause;
y = hillFunc(x,Kd,n);