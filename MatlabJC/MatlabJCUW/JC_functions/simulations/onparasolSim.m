% simulate on and off bipolars onto on parsol

% parmeters linear filters
Exc_LF_tpeak = .05 ;
Exc_LF_peakRise = .008 ;
Exc_LF_peakAmp = .12 ;
Exc_LF_ttrough = .1 ;
Exc_LF_troughDecay = .05 ;
Exc_LF_troughAmp = .02 ;

Inh_LF_tpeak = .05 ;
Inh_LF_peakRise = .008 ;
Inh_LF_peakAmp = -.12 ;
Inh_LF_ttrough = .1 ;
Inh_LF_troughDecay = .05 ;
Inh_LF_troughAmp = -.02 ;

% parameters non-linear filters
Exc_NLF_mean = 20;  % mean nonlinearity
Exc_NLF_std = 10;   % std nonlinearity
Exc_NLF_amp = 5;     % nonlinearity amplitude

Inh_NLF_mean = 20;  % mean nonlinearity
Inh_NLF_std = 5;   % std nonlinearity
Inh_NLF_amp = 8;     % nonlinearity amplitude

% parameters time (sec)
si=.0001 ;  % sample interval
maxtime = 5 ;    % sec of time 
freq = 4 ;  % sine wave frequency (hz)

% LIF model params
params.Einh = -80 ;
params.Eexc = 0 ;
params.Eleak = -60 ;
params.cap = .06 ;%nF
params.Vthresh = -45 ;
params.AbsRef  = .002 ; %sec
params.RelRefTau = .008 ;%sec
params.RelRefAmp = 4 ; %mV
params.Gleak = 0 ;
params.Vrest = -60 ;

% run simulation

% time
time = [0:si:maxtime] ; % time vector
lt= length(time) ;

% stim
Exc_light = sin((2*pi*freq)*time) ; % temporal light
Inh_light = sin((2*pi*freq)*time+pi) ; % 
    
% linear filters
Exc_LF = simFilter(time,Exc_LF_tpeak,Exc_LF_peakRise,Exc_LF_peakAmp,Exc_LF_ttrough,Exc_LF_troughDecay,Exc_LF_troughAmp) ;
Inh_LF = simFilter(time,Inh_LF_tpeak,Inh_LF_peakRise,Inh_LF_peakAmp,Inh_LF_ttrough,Inh_LF_troughDecay,Inh_LF_troughAmp) ;

% response conductances after linear filter
Exc_G_preNLF = ifft(fft(Exc_light).*fft(Exc_LF)) ;
Inh_G_preNLF = ifft(fft(Inh_light).*fft(Inh_LF)) ;

% response g after nonlinear filter
Exc_G = Exc_NLF_amp*normcdf(Exc_G_preNLF,Exc_NLF_mean,Exc_NLF_std) ;
Inh_G = Inh_NLF_amp*normcdf(Inh_G_preNLF,Inh_NLF_mean,Inh_NLF_std) ;

% spike generation
V_trace = LIFmodelG(Exc_G,Inh_G,params.Gleak,1/si,params) ;



%figures

figure
plot(time,Exc_G)
hold on
plot(time,Inh_G,'r')
plot(time,Exc_light,'b--')
plot(time,Inh_light,'r--')
plot(time,V_trace)

figure
plot(sort(Exc_G_preNLF),Exc_NLF_amp*normcdf(sort(Exc_G_preNLF),Exc_NLF_mean,Exc_NLF_std)) 
hold on
plot(sort(Inh_G_preNLF),Inh_NLF_amp*normcdf(sort(Inh_G_preNLF),Inh_NLF_mean,Inh_NLF_std),'r') 
;