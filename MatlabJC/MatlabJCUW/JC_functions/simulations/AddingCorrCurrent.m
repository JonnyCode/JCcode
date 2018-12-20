% simulation of power spectrum changes and as aperiodic signals are
% shifted and added together.

% parameters
freqcutoff = 50 ; % hz
rect = -10 ; % anything below this is removed to rectify x and y
anticorr = 1 ; % use negative 1 in order to make x and y anticorrelated

% set time (ms)  
tmax = 1000 ;           % max of time
h= .1 ;                 % time step 
time = h:h:tmax ;       % time points assessed
tpoints=length(time) ;  % number of time points

% low pass filter original wave
orig=normrnd(0,5,1,tpoints) ; % mean,std,m,n
freqcutoff = freqcutoff/((1000/h)/tpoints) ; % cutoff adjusted for length of signal 
orig_fft = fft(orig) ; % take fft of exc signal
orig_fft(:,1+freqcutoff:length(orig_fft)-freqcutoff) = 0 ; % cut out high frequencies in first and second half of fft
origfiltered = real(ifft(orig_fft)) ; % inverse fft

% make x wave the original wave filtered and rectified
x = origfiltered  ; 
x(x<rect) = rect ;

% power spectrum of x wave
[powerx_xvalues, mean_powerx] = PowerSpectrumFinder(x,10000) ; 

% for each round make y wave from orig wave add lag, sign, and rectify 
round = 0;
for shift=0:10:1000 ;
    round=round+1 ;
    y(round,:)=circshift(origfiltered,[0,shift])*anticorr ; % lag*anticorr
    y(round,y(round,:)<rect) = rect ; % rectify
    z(round,:)=x+y(round,:) ; 
    varz(round)=var(z(round,:)) ;
    [powery_xvalues(round,:), mean_powery(round,:)] = PowerSpectrumFinder(y(round,:),10000) ;
    [powerz_xvalues(round,:), mean_powerz(round,:)] = PowerSpectrumFinder(z(round,:),10000) ;

    cc(round,:) = xcorr(x,y(round,:)) ; % cross correlation of x and y

    varzPred(round) = ((abs(cc(round,tpoints))*2)/tpoints)+var(x)+var(y(round,:)) ;
end

tl = length(powerx_xvalues) ;
varx= var(x);
varx2 = 2*varx ;
varx4 = 4*varx ;

xcc = [1:length(cc)] - (length(cc)+1)/2 ; % x values for the cross corelation

figure
plot(time,origfiltered,'k')
hold on
plot(time,x)
plot(time,y(1,:),'r--')
plot(time,y(end,:),'r')

figure
plot(xcc,cc(1,:),'r--') 
hold on
plot(xcc,cc(end,:),'r')

figure
plot(0:10:1000,varz)
hold on
plot(0:10:1000,varzPred,'k')
plot(0:10:1000,varx2,'r-')
plot(0:10:1000,varx4,'r-')


round =10
figure
plot(powery_xvalues(round,:), mean_powery(round,1:tl),'r')
hold on
plot(powerz_xvalues(round,:), mean_powerz(round,1:tl),'k')
plot(powerx_xvalues, mean_powerx(1,1:tl),'g--')
H = gca ; %get current axis handle
set(H,'Xscale','log') ;
set(H,'Yscale','log') ;



% A simpler case of simulation of amplitude changes as a periodic wave is 
% shifted and added to original wave.

round = 0 ;
x=[1:.01:100];
T=3 ;

for shift=0:.1:T/2 ;
     round = round+1 ;
     y1=cos(2*pi*(x+shift)/T) ; y2 = cos(2*pi*x/T)  ;
     A(round)=var(y1+y2) ;
end

figure
plot([0:.1:T/2],A,'*-')

shift = .4
T=3, y1=cos(2*pi*(x+T/shift)/T), y2 = cos(2*pi*x/T)  ;
figure
plot(x,y1)
hold on
plot(x,y2,'r')
