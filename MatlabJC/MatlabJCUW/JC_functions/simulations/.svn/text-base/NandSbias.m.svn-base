x=[0:.0001:5] ;
numTrials = 3 ;

var = .0004 ;
m = 2.5 ;
g = exp(-(x-m).^2/(2*var))/(sqrt(var)*2*pi) ;


s= normrnd(4,1,1,length(x)) ;
sg= ifft(fft(s).*fft(g)) ;
signal = sg ;%/max(sg) ; 


for a=1:numTrials ;
    N = normrnd(0,300,1,length(x)) ;
    Noise(a,:) = lowPassFilter(N,10000,60) ;
    Trial(a,:) = signal + Noise(a,:) ;
end

for a=1:numTrials ;
    Res(a,:) = Trial(a,:)- mean(Trial) ;
end

[xP,meanP] = PowerSpectrumFinder(mean(Trial),10000) ;
[xP,sigP] = PowerSpectrumFinder(signal,10000) ;

[xP,ResP] = PowerSpectrumFinder(Res,10000) ;
[xP,NoiseP] = PowerSpectrumFinder(Noise,10000) ;

[xP,TrialP] = PowerSpectrumFinder(Trial,10000) ;

figure
plot(x,signal,'b')
hold on
plot(x,Noise(1,:),'r')

figure
plot(xP,meanP)
hold on
plot(xP,sigP,'r')
a=gca
set(a,'xscale','log')
set(a,'yscale','log')

figure
plot(xP,ResP)
hold on
plot(xP,NoiseP,'r')
a=gca
set(a,'xscale','log')
set(a,'yscale','log')


figure
plot(xP,TrialP)
a=gca
set(a,'xscale','log')
set(a,'yscale','log')

figure
plot(xP,NoiseP./ResP)
a=gca
set(a,'xscale','log')
set(a,'yscale','log')

