function y = addConeNoise(x,sampleRate,T)
%T is transform
x = T.revert(x); %revert to time domain
L = length(x);
y = x + ConeNoise(L,sampleRate);
y = T.invert(y);

