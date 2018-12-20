function y = sumAndSigmoid(x,A,K,B,v,Q,M,T)
%T is transform (not used)
x = sum(x);
%pause;
y = generalizedSigmoid(x,A,K,B,v,Q,M);