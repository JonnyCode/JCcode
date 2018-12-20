function V_trace = LIFmodelGplusI_LinSynInt(Gexc,Ginh,Gleak,Iadd,Vint,samplerate,params) 

%implements a leaky intgrate and fire model, input Gexc and Ginh should be in rows,
%Gleak and Iadd should be a constant, samplerate should be in Hz

% adapted from LIFmodelG: added constant current injection (Iadd)

% adapted from LIFmodelGplusI: kept constant driving force on exc and inh conductances at Vint (integration voltage)

% JC 6/10/11

% set intrinsic cell parameters (uncomment to run as script)

% E_inh = -80 ;%mV      % reversal potential for inhibitory current
% E_exc = 0 ;          % reversal potential for excitatory current
% E_leak = -70 ;         % reversal potential of leak conductance (this can make cell spike spontaneously)
% 
% C = 0.06 ; %nF                          % Set capacitance of cell
% V_rest = -70 ;                          % initial potential of cell (resting potential mV) 
% V_threshold = -57 ;                     % spike threshold
% abs_ref = .002*samplerate ; %pnts       % 2 ms absolute refractory period (points)
% V_relref_tau = .008 ; %sec              % decay time constant of relative refractory period (these numbers taken from Trong and Rieke, 2008)
% V_relref_amp = 4 ; %mV                  % amplitude of relative refractory period threshold change 
% Iadd (pA)

E_inh = params.Einh ;
E_exc = params.Eexc ;
E_leak = params.Eleak ;
C = params.cap ;
V_rest = params.Vrest ;
V_threshold = params.Vthresh ;
abs_ref = round(params.AbsRef * samplerate) ; 
V_relref_tau = params.RelRefTau ;
V_relref_amp = params.RelRefAmp ;

% set time (ms)  
time = 1/samplerate:1/samplerate:length(Gexc)/samplerate ;       % time points assessed
tpoints=length(time) ;  % number of time points

for trial = 1:size(Gexc,1) ; % for each exc trace

% reset variables to go for trial
V_trace(trial,:) = time ; 
V_trace(trial,1) = V_rest ;
I_syn = time ;
I_syn(1) = 0;
ref = 0 ;   
V_thr(trial,:) = ones(1,length(time))*V_threshold ;

for t = 2:tpoints %for each time point

% current from conductance
I_syn(t) = Gexc(trial,t)*(Vint-E_exc) + Ginh(trial,t)*(Vint-E_inh) + Gleak*(V_trace(trial,t-1)-E_leak) ; % note the implementation of an unchanging voltage for each time point in exc and inh currents but not leak

if ref == 0 ; % if not within refractory period
V_trace(trial, t) = V_trace(trial, t-1) + (1/samplerate)*((-I_syn(t)+Iadd)/C) ; %h needs to be in units of seconds from ms

else
V_trace(trial, t) = V_rest ;
ref = ref - 1 ;

end

if V_trace(trial,t)>V_thr(trial,t)
V_trace(trial,t) = 50 ;
V_thr(trial,[t:end]) = V_thr(trial,[t:end]) + exp([0:length(V_thr)-t]./(samplerate*-V_relref_tau))*V_relref_amp ; % add on reletive refractory period
ref = abs_ref ;
end

end % end t loop

end % end trial loop




