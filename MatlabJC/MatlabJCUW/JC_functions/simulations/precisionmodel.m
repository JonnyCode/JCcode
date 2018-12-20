%PARAMETERS
% time parameters (ms)  
tmax = 10000 ;             % max of time
h = .01 ;                 % time step 
time = 1:h:tmax ;         % time points assessed
tpoints = length(time) ;  % number of time points

% conductance parameters (ns) 
mu_exc = 300 ;
mu_inh = 300 ;
sigma_exc = 100 ;
sigma_inh = 100 ;

% conductance filter/kinetics parameters
tau_g = -2; %ms       % the decay constant of exponetial to be convolved w/ gaussian noise
filtert = [0:h:10]; %ms    % the time of the exponetial filter
lfilter = length(filtert); % the number of noise points
filter_g = exp(tau_g*filtert) ; % exponetial decay to be convolved w/ conductance noise (divided by norm in order to maintain stdev of noise signal) 

%F = exp(-((x-20).^2)/(20*3^2))/(3*sqrt(2*pi)); % POSSIBLE gaussian filter


% conductance noise parameters (ns)
noise_max = 2 ;     % the point at which mean noise stops increasing as function of g amp
noise_min = .5 ;    % the baseline noise    
noise_k = .1 ;      % the rate at which the noise increases
sigma_noise = .5 ;  % the stdev of the noise

% conductance noise filter/kinetic parameters
tau_noise = -2 ; %ms       % the decay constant of exponetial to be convolved w/ gaussian noise
filter_noise = exp(tau_noise*filtert)/max(exp(tau_noise*filtert)) ; % exponetial decay to be convolved w/ conductance noise (divided by norm in order to maintain stdev of noise signal) 

% intrinsic cell parameters 
E_inh = -70 ;%mV       % reversal potential for inhibitory current
E_exc = 10 ;           % reversal potential for excitatory current
E_leak = -60 ;         % reversal potential of leak conductance (this can make cell spike spontaneously)

C = 1 ; %nF            % Set capacitance of cell
R = .1 ; %Giga Ohms    % Set input resistance of membrane at rest (leak resistance)    

V = -60 ;            % initial potential of cell (resting potential mV) 
V_thr = -40 ;        % spike threshold
V_trace = [] ;       % voltage points    
V_trace = [V_trace, V] ;
abs_ref = 2/h ;      % absolute refractory period (2ms)
ref = 0 ;

%CONDUCTANCES
% create gaussian conductances (ns)
G_exc_gs = normrnd(0,sigma_exc,1,tpoints+lfilter) ;  
G_inh_gs = normrnd(0,sigma_inh,1,tpoints+lfilter) ;

G_exc_gs2 = conv(G_exc_gs,filter_g) ; % convolve w/ filter to create appropriate kinetics
G_inh_gs2 = conv(G_inh_gs,filter_g) ;

G_exc_stable = ((G_exc_gs2/std(G_exc_gs2))*mu_exc)+mu_exc ;



for trial = 1:100 ; % for each new set of conductance noise ... 

% create conductance noise
mu_exc_noise = noise_max-((noise_max-noise_min)./((noise_c*G_exc_stable)+1)) ;%ns  %noise is based on relation ship between g-amp and noise amp 
mu_inh_noise = noise_max-((noise_max-noise_min)./((noise_c*G_inh_stable)+1)) ;%ns  %noise is based on relation ship between g-amp and noise amp 

noise_exc = conv(filter_noise, normrnd(mu_exc_noise,sigma_noise)) ; % create gaussian noise and convolve w/ filter
noise_inh = conv(filter_noise, normrnd(mu_inh_noise,sigma_noise)) ;

% add noise to stable G
G_exc = G_exc_stable + noise_exc ;
G_inh = G_inh_stable + noise_inh ;

for t = 1:tpoints % for each point in time

% LIF SPIKE GENERATION 
% calculate synaptic current
I_syn = G_inh(t)*(V-E_inh) + G_exc(t)*(V-E_exc) ;%pA

% implement integrate and fire model c*dV/dt=Sum(g(V-E))
% dV/dt = - V/RC + I/C
if ~ref                    % if not within refrectory period... 
    V = V - (h*(-(E_leak/R*C) +(V/(R*C)) + (I_syn/C))) ;    %calculate the voltage
else
    V = V_trace(1) ;       % reset voltage to rest
    ref = ref-1 ;
end

% Spike
if (V > V_thr)      % if voltage is higher than threshold...
    V = 50 ;        % spike and make the next point at rest
    ref = abs_ref ;
end

% track voltage
V_trace = [V_trace, V] ;

end % this ends t loop 

% Note times of spikes
Spiketimes{trial}= find(V_trace == 50)*h ;

end % this ends trial loop

% SPIKE PRECISION ANALYSIS
% create PSTH of spiketimes


% find peak times in PSTH
PSTH_thrsh =

% find spike-triggered conductances

%

