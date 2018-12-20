function T = makeTuningCurves(N)
%N is how many sine wave tuning curves to make
%we use 1 degree resolution for now

T = zeros(N,360);
x = 0:1:359;

phaseOffset = 360./N;
for i=1:N
   T(i,:) = (sin(((i-1).*phaseOffset+x)*pi/180) + 1)/2; %so minimum is zero max is 1       
end