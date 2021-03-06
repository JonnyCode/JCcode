function stats = precision_model(f_inh,f_exc,amp_inh,amp_exc,phase) ;


% This function will implement a basic integrate and fire model for a set of
% noisy conductances.  Hopefully, this model will help answer the
% question: What synaptic inputs will permit high spike precision regardless
% of intrinsic mechanisms?

%NOTE- sections between lines are helpful when running function with only a
%single set of conducatance parameters.

% J. Cafaro 5/13/07

% clear all
% close all

% set time (ms)  
tmax = 1000 ;           % max of time
h= .01 ;                 % time step 
time = 1:h:tmax ;       % time points assessed
tpoints=length(time) ;  % number of time points
       
% %__________________________________________________________________________
% %set conductance parameters
% f_inh = .03 ;          % frequency of inhibitory conductance 
% f_exc = .03 ;          % frequency of excitatory conductance
% 
% amp_inh = 2 ;%ns       % amplitude of inhibitory synaptic coductance 
% amp_exc = 10 ;%ns       % amplitude of excitatory synaptic conductance 
% 
% phase = pi ;           % relative phase of sinusiods of conductances
% %__________________________________________________________________________

% set conductance noise parameters 
noise_max = 2 ;%ns   % the point at which noise stops increasing as function of g amp
noise_min = .5 ;%ns      
noise_start = 0 ;%ns   % the amp of the conductance at which the noise begins growing 
noise_stop = 5 ;%ns    % the amp of the conductance at which the noise maxes out (sigmoid points)  

tau = -2; %ms           % the decay constant of exponetial to be convolved w/ gaussian noise
noiset = [0:h:10]; %ms  % the time of the exponetial filter
lnoiset = length(noiset); % the number of noise points

mu_inh = 0 ;%nS       % mean of noise for inhibitory conductance
mu_exc = 0 ;          % mean of noise for excitatory conductance

sigma_inh = (smf(amp_inh,[noise_start,noise_stop])*(noise_max-noise_min))+noise_min ;  % standard deviation of noise for inhibitory conductance
sigma_exc = (smf(amp_exc,[noise_start,noise_stop])*(noise_max-noise_min))+noise_min ;  % these values derived from sigmoidal relationship between g-amp and noise amp

% set intrinsic cell properties
E_inh = -70 ;%mV      % reversal potential for inhibitory current
E_exc = 10 ;          % reversal potential for excitatory current
E_leak = -60 ;         % reversal potential of leak conductance (this can make cell spike spontaneously)

C = 1 ; %nF            % Set capacitance of cell
R = .1 ; %Giga Ohms    % Set input resistance of membrane at rest (leak resistance)    

V = -60 ;            % initial potential of cell (resting potential mV) 
V_thr = -40 ;        % spike threshold
V_trace = [] ;       % voltage points    
V_trace = [V_trace, V] ;
abs_ref = 2/h ;        % absolute refractory period (2ms)
ref = 0 ;

% set excititoatory and inhibitory sinusiodal conductances
G_sin_inh = amp_inh*sin(time*f_inh)+amp_inh ;

G_sin_exc = amp_exc*sin((time*f_exc)+phase)+amp_exc ;

% set the number of trials
max_rounds = 2 ;

for round = 1:max_rounds 
% set noise for conductances
G_noise_inh = normrnd(mu_inh,sigma_inh,1,tpoints+lnoiset) ;  % normrnd = draws from normal distribution with mean(mu) & stdev(sigma) 
G_noise_exc = normrnd(mu_exc,sigma_exc,1,tpoints+lnoiset) ;  

noise_smoother = exp(tau*noiset)/norm(exp(tau*noiset)) ; % exponetial decay to be convolved w/ conductance noise (divided by norm in order to maintain stdev of noise signal)  
G_noise_inh_smooth = conv(G_noise_inh,noise_smoother);
G_noise_exc_smooth = conv(G_noise_exc,noise_smoother);

G_inh = max((G_sin_inh + G_noise_inh_smooth(lnoiset+1:(tpoints+lnoiset))),0) ; % max(x,0)= rectifies conductance 
G_exc = max((G_sin_exc + G_noise_exc_smooth(lnoiset+1:(tpoints+lnoiset))),0) ; % this selects only middle of noise to avoid convolution effects

for t = 1:tpoints

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
Spiketimes{round}= find(V_trace == 50)*h ;

% %__________________________________________________________________________
% % plot voltage
% figure(1)
% plot([0,time],V_trace + (ones(1,length(V_trace))*150*round))
% hold on
% title('voltage')
% %__________________________________________________________________________

% reset variables to go for next round
V = -60 ;  
V_trace = [] ; 
V_trace = [V_trace, V] ;
ref = 0 ;

end % this ends the round loop

% assess spike precision
cost = .05 ;                                % set cost of spike distance metric (shift if <40ms)
count = 1 ;                                 
for a=1:size(Spiketimes,2)-1 ;               % for each spike train trial...                                       
    for b=a+1:size(Spiketimes,2) ;           % compare it with another spike train only once
d(count)=spkd(Spiketimes{a},Spiketimes{b},cost) ;    % runs spike distance function at a given cost
d_norm(count)=d(count)/max(length(Spiketimes{a}),length(Spiketimes{b})) ;    % normalize by the number of spikes    
count = count+1 ;
    end
end
average_spkd = mean(d);
stdev_spkd = std(d) ;
average_spkd_norm = mean(d_norm) ;
stdev_spkd_norm = std(d_norm) ;

% assemble important numbers 

stats=[f_inh;f_exc;amp_inh;amp_exc;phase;average_spkd_norm;stdev_spkd_norm];

% %__________________________________________________________________________
% % plot conductances
% figure
% plot(time,G_sin_inh)
% hold on
% plot(time,G_sin_exc,'r')
% title('base conductances')
% 
% figure
% plot(G_noise_inh_smooth)
% hold on
% plot(G_noise_exc_smooth,'r')
% title('conductance noise')
% 
% figure
% plot(time,G_inh)
% hold on
% plot(time,G_exc,'r')
% title('conductances with noise')
% %__________________________________________________________________________

end % this ends the function

