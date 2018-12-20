% Basic integrate-and-fire neuron
% B Lundstrom, Summer 2003
%
% Given an input current, membrane potential (V) is determined 
% as follows:
% 1) V = 50 mV during a spike, i.e. when V > V_th = 10 mV
% 2) V = V_th - 15 mV during refractory period
% 3) V = V - (V-E)/(R*C) + I/C, otherwise
%
% (3) is Euler apporx derived from Iin = Icap + Ires,
% which yields CdV/dt = -(V-E)/R + Iin;
% further, this ODE can be solved via variation of parameters 
% leading to V(t) = E + IR + (Vi - E - IR)*exp(-t/RC)

clear

% capacitance (nF) and resistance (M ohms)
C = 1; %nF, assuming the standard 10 nF/mm2, we have a .1 mm2 neuron
R = 40; %Mohm
El = -65;
time = 100; %ms length of trial
ref = 0; % counter for refractory period
V_th = -50;

I = 2; % nA

% create 'look-up' table according to (3) above
% (more rows can be added for multiple currents)
V = zeros(size(I,1), time);
    
V(1) = El + I(1)/C;
for c = 2:time
        
    if ~ref
        V(c) = V(c-1) + (El - V(c-1))/(R*C) + I/C;
    else
        ref = ref - 1;
        V(c) = -75;
    end
    if V(c) > V_th
        V(c) = 30;
       ref = 5;
    end
end
  
x = 1:time;
plot(x,V)
title(['V (mV) vs. time (ms), I = ' int2str(I) ' nA/.1mm2'])
