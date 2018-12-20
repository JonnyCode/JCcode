% simple simulation to help understand the frequency domain, make sure to
% put lowPassFilter on your path if running the noiseSignal
% JC 1/20/12

% select signal type
sineSignal = true ; % if true than create a sum of sine waves
noiseSignal = false ; % if true than create gaussian noise
gaussianSignal = false ; % if true than create a gaussian

% time parameters
sampleRate = 10000 ; %hz
maxTime = 50 ; %sec

%sumed sine wave parameters
SineFrequency = [1,3,5,7,9] ; %hz
SinePhase = [90,0,0,0,0] ; % degrees
SineAmplitude = [1,1/3,1/5,1/7,1/9] ; 

% noise parameters
noiseVar = 500 ; % variance
noiseMaxFreq = 4 ; % maximum high freq

% gaussian parameters
gaussMean = 10 ;
gaussVar = .5 ;

% time
time = [0:1/sampleRate:maxTime] ;
numPnts = length(time) ;

% construct signal
if sineSignal
    for a=1:length(SineFrequency) ;
        sineWave(a,:) = SineAmplitude(a)*sin(SineFrequency(a)*2*pi*time-SinePhase(a)*pi/180) ;
    end
    signal = sum(sineWave,1) ; 
    
elseif noiseSignal 
    signalTemp = normrnd(0,sqrt(noiseVar),1,length(time)) ;
    signal = lowPassFilter(signalTemp,sampleRate,noiseMaxFreq) ;
    
elseif gaussianSignal
    signal = exp(-((time-gaussMean).^2)/(2*gaussVar)) ;
    
end

% analyzable frequencies
niquist = sampleRate/2 ;
fstep = 1/maxTime ;
powerSpecX = [0:fstep:niquist] ;

% dft of signal
signal_dft_temp = fft(signal) ;
signal_dft = signal_dft_temp(1:length(powerSpecX)) ; % only the positive frequency components

% analyze signal variance
signal_dft_amplitude = sqrt(real(signal_dft).^2+imag(signal_dft).^2) ;
% signal_dft_amplitude = sqrt(signal_dft*conj(signal_dft)) ; % equivilent to above
scaleFactor = 2/(numPnts^2) ;
signal_dft_power = signal_dft_amplitude.^2*scaleFactor ;

signal_dft_power_sum = sum(signal_dft_power) 
signal_var = var(signal) 

% analyze signal phase
signal_dft_phase = atand(real(signal_dft)./imag(signal_dft)) ;
signal_dft_phase = signal_dft_phase(1:length(powerSpecX)) ;

% figures
figure 
subplot(3,2,1) % time domain
plot(time,signal) 
xlabel('time')
ylabel('amplitude')
title('time domain')

if sineSignal
    hold on
    plot(time,sineWave,'r')
end

subplot(3,2,2) % frequency domain
plot(real(signal_dft),imag(signal_dft),'*')
xlabel('dft real parts')
ylabel('dft imag parts')
title('frequency domain')

subplot(3,2,3:4) % power spectrum
loglog(powerSpecX,signal_dft_power,'*')
xlabel('frequency')
ylabel('amplitude^2')
title('power spectrum')

subplot(3,2,5:6) % phase spectrum
stem(powerSpecX,signal_dft_phase,'filled')
h = gca ;
set(h,'xscale','log')
xlabel('frequency')
ylabel('phase (degrees)')
title('phase spectrum')

