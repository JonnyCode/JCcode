function stats = precision_model(f_inh,f_exc,amp_inh,amp_exc,phase)

% This function will implement a basic integrate and fire model for a set of
% noisy conductances.  Hopefully, this model will help answer the
% question: What synaptic inputs will permit high spike precision regardless
% of intrinsic mechanisms?

%NOTE- sections between lines are helpful when running function with only a
%single set of conducatance parameters.

% J. Cafaro 5/13/07

clear all

%set time  
tmax = 200 ;            % max of time
time = 1:.1:tmax ;      % time points assessed
tpoints=length(time) ;  % number of time points

%__________________________________________________________________________
% %set conductance parameters
% f_inh = .1 ;          % frequency of inhibitory conductance 
% f_exc = .1 ;          % frequency of excitatory conductance
% 
% amp_inh = .01 ;        % amplitude of inhibitory synaptic coductance
% amp_exc = .05 ;        % amplitude of excitatory synaptic conductance
% 
% phase = pi ;          % relative phase of sinusiods of conductances
%__________________________________________________________________________

E_inh = -70 ;         % reversal potential for inhibitory current
E_exc = 10 ;          % reversal potential for excitatory current

mu_inh = .002 ;          % mean of noise for inhibitory conductance
mu_exc = .002 ;          % mean of noise for excitatory conductance

sigma_inh = .002 ;        % standard deviation of noise for inhibitory conductance
sigma_exc = .002 ;        % standard deviation of noise for excitatory conductance

% set intrinsic cell properties
E_leak = -60 ;         % reversal potential of leak conductance (this can make cell spike spontaneously)
C = 1 ; %nF           Set capacitance of cell
R = 50 ; %Mega Ohms   Set input resistance of membrane at rest (leak resistance)    

V = -60 ;            % initial potential of cell (resting potential) 
V_thr = -40 ;        % spike threshold
V_trace = [] ;       % voltage points    
V_trace = [V_trace, V] ;
abs_ref = 2 ;        % absolute refractory period
ref = 0 ;

% set excititoatory and inhibitory sinusiodal conductances
G_sin_inh = amp_inh*sin(time*f_inh)+amp_inh ;

G_sin_exc = amp_exc*sin((time*f_exc)+phase)+amp_exc ;

% set the number of trials
max_rounds = 10 ;

for round = 1:max_rounds 
% set noise for conductances
G_noise_inh = abs(normrnd(mu_inh,sigma_inh,1,tpoints)) ; 
G_noise_exc = abs(normrnd(mu_exc,sigma_exc,1,tpoints)) ;

G_inh = G_sin_inh + G_noise_inh ; 
G_exc = G_sin_exc + G_noise_exc ;

for t = 1:length(time)

% calculate synaptic current
I_syn = G_inh(t)*(V-E_inh) + G_exc(t)*(V-E_exc) ;

% implement integrate and fire model c*dV/dt=Sum(g(V-E))
% dV/dt = - V/RC + I/C
if ~ref                    % if not within refrectory period... 
    V = V - (-(E_leak/R*C) +(V/(R*C)) + (I_syn/C)) ;    %calculate the voltage
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
Spiketimes{round}= find(V_trace == 50) ;

%__________________________________________________________________________
% % plot voltage
% figure(1)
% plot(V_trace + (ones(1,length(V_trace))*80*round))
% hold on
% title('voltage')
%__________________________________________________________________________

% reset variables to go for next round
V = -60 ;  
V_trace = [] ; 
V_trace = [V_trace, V] ;
ref = 0 ;

end % this ends the round loop

% assess spike precision
cost = 3 ;                                   % set cost of spike distance metric
count = 1 ;                                 
for a=1:size(Spiketimes,2)-1 ;               % for each spike train trial...                                       
    for b=a+1:size(Spiketimes,2) ;           % compare it with another spike train only once
d(count)=spkd(Spiketimes{a},Spiketimes{b},cost) ;    % runs spike distance function at a given cost
    count = count+1 ;
    end
end
average_spkd = mean(d);
stdev_spkd = std(d) ;

% assemble important numbers 

stats=[f_inh;f_exc;amp_inh;amp_exc;phase;average_spkd;stdev_spkd];

%__________________________________________________________________________
% % plot conductances
% figure
% plot(time,G_sin_inh)
% hold on
% plot(time,G_sin_exc,'r')
% hold on
% plot(time,G_noise_inh)
% hold on
% plot(time,G_noise_exc,'r')
% title('sin conductances and noise')
% 
% figure
% plot(time,G_inh)
% hold on
% plot(time,G_exc,'r')
% title('conductances with noise')
%__________________________________________________________________________

end % this ends the function

