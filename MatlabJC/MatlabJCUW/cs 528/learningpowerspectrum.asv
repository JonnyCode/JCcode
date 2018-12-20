%this script to help me understand fourier transforms, filtering, and power
%spectrums
%JC 6/11/07


clear all
close all

samplepoints=10000 %Hz
t = 0:.1:1000; %ms time points evaluated
x = sin(2*pi*.01*t)+sin(2*pi*1*t)+sin(2*pi*.1*t); %sinwaves
figure,plot(t,x);
maxfreq = floor(.5*length(x))

F = fft(x);
powerF = F.*conj(F)/length(x);
f= samplepoints*(0:maxfreq)/length(x);
figure, plot(f,powerF(1:5001));



%filtering
lowpass = 500
highpass = 1001
filterF=F;
filterF(lowpass+1:highpass)=0;        % fft has two signals, with the highest frequencies in the middle of F
filterF(10001-highpass+1:10001-lowpass)=0;       % therefore you must filter out both or it won't work


filterpowerF = filterF.*conj(filterF)/length(x) ;
f= samplepoints*(0:maxfreq)/length(x);
figure, plot(f,filterpowerF(1:5001));

filterx = ifft(filterF)
figure, plot(t,filterx)


