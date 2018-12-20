function V_trace = EIFmodelGplusGleakAndIadd(Gexc,Ginh,Gleak,Iadd,V,fV,samplerate,params) 

% modified from EIFmodelG.m: added  

%implements a exponential intgrate and fire model, (see Badel  et. al.
%2008), input Gexc and Ginh should be in rows, samplerate should be in Hz,
%fV models the intrinsic conductances as a fucntion of V (voltage).

% JC 5/16/11

% set intrinsic cell parameters (uncomment to run as script)

% E_inh = -80 ;%mV      % reversal potential for inhibitory current
% E_exc = 0 ;          % reversal potential for excitatory current
% 
% C = 0.06 ; %nF                          % Set capacitance of cell
% V_rest = -70 ;                          % initial potential of cell (resting potential mV) 
% V_threshold = -57 ;                     % spike threshold
% abs_ref = .002*samplerate ; %pnts       % 2 ms absolute refractory period (points)
% V_relref_tau = .008 ; %sec              % decay time constant of relative refractory period (these numbers taken from Trong and Rieke, 2008)
% V_relref_amp = 4 ; %mV                  % amplitude of relative refractory period threshold change 

E_inh = params.Einh ;
E_exc = params.Eexc ;
E_leak = params.Eleak ;
C = params.cap ;
V_rest = params.Vrest ;
V_threshold = params.Vthresh ;
abs_ref = params.AbsRef * samplerate ; 
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
        I_syn(t) = Gexc(trial,t)*(V_trace(trial,t-1)-E_exc) + Ginh(trial,t)*(V_trace(trial,t-1)-E_inh) + Gleak*(V_trace(trial,t-1)-E_leak) ;

        if ref == 0 ; % if not within refractory period
            V_int = interp1(V,fV,V_trace(trial, t-1),'linear','extrap') ; % calculate the intrinsic V from fV (mV/sec)
            V_trace(trial, t) = V_trace(trial, t-1) + (1/samplerate)*(((-I_syn(t)+Iadd)/C) + V_int) ; %

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
