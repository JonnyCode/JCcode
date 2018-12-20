%% ===== Initializes cameras and data acquisition ============================

clear all; 
daqreset; 
global fromEphys

% parameters (TURN THIS INTO A STRUCT PASSED TO A FUNCTION SELECTOR)
experiment_number = '2' ;
numsecs =5 ; % for each trial
numRepeats = 5 ; % number of repeated trials
TimeBetweenRepeats = 10 ; % (sec)
sampleRate = 10000 ; % hz
InjPulseInitTime = .5 ; % (sec) start time of Rinput test 
InjPulseDur = .5 ; % (sec) duration of Rinput test
InjPulseAmp = -0.05 ; % (volts) from daq to external clamp signal 
OdorPulseInitTime = 3 ; % (sec) start time of odor pulse
OdorPulseDur = .5 ; % (sec) duration of odor pulse
OdorPulseAmp = 10 ; % (volts) from daq to air valve
updateFig = false ; % true if you want to update figures every second (false otherwise)

% set daq output
fromEphys.lt = numsecs*sampleRate;
outputData = zeros(fromEphys.lt,2) ; % ao0

% odor
outputData(OdorPulseInitTime*sampleRate:(OdorPulseInitTime+OdorPulseDur)*sampleRate,1)= OdorPulseAmp ;

% Vext
outputData(InjPulseInitTime*sampleRate:(InjPulseInitTime+InjPulseDur)*sampleRate,2)= InjPulseAmp ; %step
% outputData(:,2) = normrnd(0,InjPulseAmp,1,length(t)) ; % noise injection
% outputData(:,2) = lowPassFilter(outputData(:,2)',sampleRate,60) ;
% outputData(end-10:end,2) = 0 ;

% prep data file  
rootdir = ['C:\Cafaro Data\' datestr(date, 'yymmdd'),'Data'];
exptag = [datestr(now,'yymmdd') '_' experiment_number];

mkdir(rootdir); % make root directory
cd(rootdir);    

mkdir(rootdir,exptag); % make experiment directory
cd([rootdir,'\',exptag]);


for a=[1:numRepeats] ;

    D = dir(['voltage*.mat']);
    if isempty(D)
        n = 0 ;
    else 
        n = length(D) ;
    end   

    fromEphys.n = n ;
    fromEphys.voltagefile = ['voltage_' exptag '_' int2str(n)];
    fromEphys.currentfile = ['current_' exptag '_' int2str(n)];
    fromEphys.Ao0file = ['Ao0_' exptag '_' int2str(n)];
    fromEphys.Ao1file = ['Ao1_' exptag '_' int2str(n)];
    fromEphys.timefile = ['TrigTime_' exptag '_' int2str(n)];
    
    % daq object for session
    s = daq.createSession('ni');
    s.Rate= sampleRate;

    % set daq channels
    s.addAnalogOutputChannel('Dev1','ao0','Voltage');
    s.addAnalogOutputChannel('Dev1','ao1','Voltage');

    s.addAnalogInputChannel( 'Dev1','ai0','voltage');
    s.addAnalogInputChannel( 'Dev1','ai1','voltage');
    s.addAnalogInputChannel( 'Dev1','ai2','voltage');
    s.addAnalogInputChannel( 'Dev1','ai3','voltage');
    s.addAnalogInputChannel( 'Dev1','ai4','voltage');
    s.addAnalogInputChannel( 'Dev1','ai5','voltage');
    s.addAnalogInputChannel( 'Dev1','ai6','voltage');

    s.queueOutputData(outputData) ;
    
    if updateFig ;
        s.NotifyWhenDataAvailableExceeds = sampleRate ; 
    else
        s.NotifyWhenDataAvailableExceeds = fromEphys.lt ;
    end
        
    lh1 = s.addlistener('DataAvailable', @savedataJCv2);

    s.startBackground() ; 

    s.wait() ;
    
    clear lhl s
    pause(TimeBetweenRepeats) ;
end

