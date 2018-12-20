%% original script
% Fire a neuron via alpha function synapse and random input spike train
% R Rao 2007

clear
rand('state',0)
% I & F implementation dV/dt = - V/RC + I/C
h = 1; % step size, Euler method, = dt ms
t_max= 200; % ms, simulation time period
tstop = t_max/h; % number of time steps
ref = 0; % refractory period counter

% Generate random input spikes
% Note: This is not entirely realistic - no refractory period
% Also: if you change step size h, input spike train changes too...
spike_train = rand(tstop,1);
thr = 0.9; % threshold for random spikes
spike_train(find(spike_train > thr)) = ones(size((find(spike_train > thr))));
spike_train(find(spike_train < thr)) = zeros(size((find(spike_train < thr))));

% alpha func synaptic conductance
t_a = 100; % Max duration of syn conductance
t_peak = 1; % ms
g_peak = 0.05; % nS (peak synaptic conductance)
const = g_peak/(t_peak*exp(-1)); 
t_vec = 0:h:t_a;
alpha_func = const*t_vec.*(exp(-t_vec/t_peak));
clf
plot(t_vec(1:80),alpha_func(1:80))
xlabel('t (in ms)')
title('Alpha Function (Synaptic Conductance for Spike at t=0)')
pause(2) 

% capacitance and leak resistance
C = 0.5 % nF
R = 40 % M ohms

% conductance and associated parameters to simulate spike rate adaptation
g_ad = 0; 
G_inc = 1/h;
tau_ad = 2;

% Initialize basic parameters
E_leak = -60; % mV, equilibrium potential
E_syn = 0; % Excitatory synapse (why is this excitatory?)
g_syn = 0; % Current syn conductance
V_th = -40; % spike threshold mV
V_spike = 50; % spike value mV
ref_max = 4/h; % Starting value of ref period counter
t_list = [];
V = E_leak;
V_trace = [V];
t_trace = [0];

clf
subplot(2,1,1)
plot(0:h:t_max,[0; spike_train])
title('Input spike train')

for t = 1:tstop

   % Compute input
   if (spike_train(t) > 0) % check for input spike
    t_list = [t_list; 1];
   end   
   % Calculate synaptic current due to current and past input spikes
   g_syn = sum(alpha_func(t_list));
   I_syn = g_syn*(E_syn - V); 

   % Update spike times
   if t_list
     t_list = t_list + 1;
     if (t_list(1) == t_a)  % Reached max duration of syn conductance
       t_list = t_list(2:max(size(t_list)));
     end
   end

   % Compute membrane voltage
   % Euler method: V(t+h) = V(t) + h*dV/dt
   if ~ref
     V = V + h*(- ((V-E_leak)*(1+R*g_ad)/(R*C)) + (I_syn/C));
     g_ad = g_ad + h*(- g_ad/tau_ad); % spike rate adaptation
   else
     ref = ref - 1;
     V = V_th-10; % reset voltage after spike
     g_ad = 0;
   end
   
   % Generate spike
   if ((V > V_th) & ~ref)
     V = V_spike;
     ref = ref_max;
     g_ad = g_ad + G_inc;
   end

   V_trace = [V_trace V];
   t_trace = [t_trace t*h];   
end

subplot(2,1,2)
plot(t_trace,V_trace)
drawnow
title('Output spike train')







%% QUESTION 3b

% Fire a neuron via alpha function synapse and random input spike train
% R Rao 2007

clear
count=1
round=1
rand('state',0)


%for t_peak=.5:.5:10     %for question 3a

for round = 1:5
t_peak=.5;
for thr=.7:.05:.95      %for question 3b
    
% I & F implementation dV/dt = - V/RC + I/C
h = 1; % step size, Euler method, = dt ms
t_max= 200; % ms, simulation time period
tstop = t_max/h; % number of time steps
ref = 0; % refractory period counter

% Generate random input spikes
% Note: This is not entirely realistic - no refractory period
% Also: if you change step size h, input spike train changes too...
spike_train = rand(tstop,1);
%thr = 0.9; % threshold for random spikes
spike_train(find(spike_train > thr)) = ones(size((find(spike_train > thr))));
spike_train(find(spike_train < thr)) = zeros(size((find(spike_train < thr))));

% alpha func synaptic conductance
t_a = 100; % Max duration of syn conductance
%t_peak = 1; % ms
g_peak = 0.05; % nS (peak synaptic conductance)
const = g_peak/(t_peak*exp(-1)); 
t_vec = 0:h:t_a;
alpha_func = const*t_vec.*(exp(-t_vec/t_peak));
clf
plot(t_vec(1:80),alpha_func(1:80))
xlabel('t (in ms)')
title('Alpha Function (Synaptic Conductance for Spike at t=0)')
pause(2) 

% capacitance and leak resistance
C = 0.5; % nF
R = 40; % M ohms

% conductance and associated parameters to simulate spike rate adaptation
g_ad = 0; 
G_inc = 1/h;
tau_ad = 2;

% Initialize basic parameters
E_leak = -60; % mV, equilibrium potential
E_syn = 0; % Excitatory synapse (why is this excitatory?)
g_syn = 0; % Current syn conductance
V_th = -40; % spike threshold mV
V_spike = 50; % spike value mV
ref_max = 4/h; % Starting value of ref period counter
t_list = [];
V = E_leak;
V_trace = [V];
t_trace = [0];

clf
subplot(2,1,1)
plot(0:h:t_max,[0; spike_train])
title('Input spike train')

for t = 1:tstop

   % Compute input
   if (spike_train(t) > 0) % check for input spike
    t_list = [t_list; 1];
   end   
   % Calculate synaptic current due to current and past input spikes
   g_syn = sum(alpha_func(t_list));
   I_syn = g_syn*(E_syn - V); 

   % Update spike times
   if t_list
     t_list = t_list + 1;
     if (t_list(1) == t_a)  % Reached max duration of syn conductance
       t_list = t_list(2:max(size(t_list)));
     end
   end

   % Compute membrane voltage
   % Euler method: V(t+h) = V(t) + h*dV/dt
   if ~ref
     V = V + h*(- ((V-E_leak)*(1+R*g_ad)/(R*C)) + (I_syn/C));
     g_ad = g_ad + h*(- g_ad/tau_ad); % spike rate adaptation
   else
     ref = ref - 1;
     V = V_th-10; % reset voltage after spike
     g_ad = 0;
   end
   
   % Generate spike
   if ((V > V_th) & ~ref)
     V = V_spike;
     ref = ref_max;
     g_ad = g_ad + G_inc;
   end

   V_trace = [V_trace V];
   t_trace = [t_trace t*h];   
end

num_spikes = length(find(V_trace==50));
Firing_rate(round,count) = num_spikes*(1000/tstop);

num_input_spikes = sum(spike_train);
Firing_rate_input(round,count) = num_input_spikes*(1000/tstop);
count=count+1;
end
count=1
round=round+1;
end

mean_Firing_rate_input = mean(Firing_rate_input);
mean_Firing_rate = mean(Firing_rate);
figure, plot(mean_Firing_rate_input,mean_Firing_rate)
xlabel('mean input firing rate')
ylabel('mean output firing rate')



figure, plot([.5:.5:10],Firing_rate)
xlabel('t peak')
ylabel('Firing rate (Hz)')


subplot(2,1,2)
plot(t_trace,V_trace)
drawnow
title('Output spike train')









%% Question 3c
% Fire a neuron via alpha function synapse and random input spike train
% R Rao 2007

clear
count=1
rand('state',0)

t_peak_range=1:2:15
for t_peak=t_peak_range 

    % I & F implementation dV/dt = - V/RC + I/C
h = 1; % step size, Euler method, = dt ms
t_max= 200; % ms, simulation time period
tstop = t_max/h; % number of time steps
ref = 0; % refractory period counter

% Generate random input spikes
% Note: This is not entirely realistic - no refractory period
% Also: if you change step size h, input spike train changes too...
spike_train = rand(tstop,1);
thr = 0.9; % threshold for random spikes
spike_train(find(spike_train > thr)) = ones(size((find(spike_train > thr))));
spike_train(find(spike_train < thr)) = zeros(size((find(spike_train < thr))));

% alpha func synaptic conductance
t_a = 100; % Max duration of syn conductance
%t_peak = 1; % ms
g_peak = 0.05; % nS (peak synaptic conductance)
const = g_peak/(t_peak*exp(-1)); 
t_vec = 0:h:t_a;
alpha_func = const*t_vec.*(exp(-t_vec/t_peak));
clf
plot(t_vec(1:80),alpha_func(1:80))
xlabel('t (in ms)')
title('Alpha Function (Synaptic Conductance for Spike at t=0)')
pause(2) 

% capacitance and leak resistance
C = 0.5; % nF
R = 40; % M ohms

% conductance and associated parameters to simulate spike rate adaptation
g_ad = 0; 
G_inc = 1/h;
tau_ad = 2;

% Initialize basic parameters
E_leak = -60; % mV, equilibrium potential
E_syn = -70; % Excitatory synapse (why is this excitatory?)
g_syn = 0; % Current syn conductance
V_th = -40; % spike threshold mV
V_spike = 50; % spike value mV
ref_max = 4/h; % Starting value of ref period counter
t_list = [];
V = E_leak;
V_trace = [V];
t_trace = [0];

clf
subplot(2,1,1)
plot(0:h:t_max,[0; spike_train])
title('Input spike train')

for t = 1:tstop

   % Compute input
   if (spike_train(t) > 0) % check for input spike
    t_list = [t_list; 1];
   end   
   % Calculate synaptic current due to current and past input spikes
   g_syn = sum(alpha_func(t_list));
   constant_current_inj = 2.7 ; % this gets it spiking around 150Hz
   I_syn = g_syn*(E_syn - V) + constant_current_inj;    % JC added to I_syn to simulate constant current input 
   %I_syn = constant_current_inj
   
   % Update spike times
   if t_list
     t_list = t_list + 1;
     if (t_list(1) == t_a)  % Reached max duration of syn conductance
       t_list = t_list(2:max(size(t_list)));
     end
   end

   % Compute membrane voltage
   % Euler method: V(t+h) = V(t) + h*dV/dt
   if ~ref
     V = V + h*(- ((V-E_leak)*(1+R*g_ad)/(R*C)) + (I_syn/C));
     g_ad = g_ad + h*(- g_ad/tau_ad); % spike rate adaptation
   else
     ref = ref - 1;
     V = V_th-10; % reset voltage after spike
     g_ad = 0;
   end
   
   % Generate spike
   if ((V > V_th) & ~ref)
     V = V_spike;
     ref = ref_max;
     g_ad = g_ad + G_inc;
   end

   V_trace = [V_trace V];
   t_trace = [t_trace t*h];   
end

num_spikes = length(find(V_trace==50));
Firing_rate(count) = num_spikes*(1000/tstop);

num_input_spikes = sum(spike_train);
Firing_rate_input(count) = num_input_spikes*(1000/tstop);
count=count+1;
end

figure, plot(t_peak_range,Firing_rate)
xlabel('t peak')
ylabel('Output Firing rate')


%%
