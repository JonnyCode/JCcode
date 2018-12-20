% Basic integrate-and-fire neuron 
% R Rao 2007

clear
tstop=200;
count=1
% input current
f_range=1:.1:100;
for f= f_range 
I = sin((1:tstop)*f); %for question 2d
%figure, plot(I)

% I_range = .25:.1:10 ;
%for I = I_range   % for question 2c

%I = 1000 % nA   % for question 2

% capacitance and leak resistance
C = 1 % nF
R = 40 % M ohms

% I & F implementation dV/dt = - V/RC + I/C
% Using h = 1 ms step size, Euler method

V = 0;
tstop = 200;
abs_ref = 5; % absolute refractory period 
ref = 0; % absolute refractory period counter
V_trace = []; % voltage trace for plotting
V_th = 10; % spike threshold

for t = 1:tstop
  
   if ~ref
     V = V - (V/(R*C)) + (I(t)/C);
   else
     ref = ref - 1;
     V = 0.2*V_th; % reset voltage
   end
   
   if (V > V_th)
     V = 50;  % emit spike
     ref = abs_ref; % set refractory counter
   end

   V_trace = [V_trace V];

end

  %figure, plot(V_trace)

  
 % Jon's additions:

 % calculate FI curve
num_spikes = length(find(V_trace==50));
Firing_rate(count) = num_spikes*(1000/tstop);
count=count+1
end

figure, plot(f_range,Firing_rate)
xlabel('input frequency')
ylabel('firing rate (Hz)')


% plot FI curve
figure, plot(I_range,Firing_rate)
xlabel('current')
ylabel('Firing rate')

figure, plot(f_range,Firing_rate)
xlabel('frequency of current input')
ylabel('Firing rate')







