%% ===== Initializes cameras and data acquisition ============================

clear all; 
daqreset; 

% parameters
experiment_number = '1' ;
numsecs =5 ;
sampleRate = 10000 ; % hz
InjPulseInitTime = .5 ; % (sec) start time of Rinput test 
InjPulseDur = 0.5 ; % (sec) duration of Rinput test
InjPulseAmp = -0.05 ; % (volts) from daq to external clamp signal 
OdorPulseInitTime = 3 ; % (sec) start time of odor pulse
OdorPulseDur = .5 ; % (sec) duration of odor pulse
OdorPulseAmp = 7 ; % (volts) from daq to air valve

% set daq output
t = [0:(1/sampleRate):numsecs] ;
outputData(:,1)= zeros(1,length(t)) ; % ao0
outputData(OdorPulseInitTime*sampleRate:(OdorPulseInitTime+OdorPulseDur)*sampleRate,1)= OdorPulseAmp ;

outputData(:,2)= zeros(1,length(t)) ; % ao1
outputData(InjPulseInitTime*sampleRate:(InjPulseInitTime+InjPulseDur)*sampleRate,2)= InjPulseAmp ; %step
% outputData(:,2) = normrnd(0,InjPulseAmp,1,length(t)) ; % noise injection
% outputData(:,2) = lowPassFilter(outputData(:,2)',sampleRate,30) ;
% outputData(end-10:end,2) = 0 ;

% prep data file  
rootdir = ['C:\Cafaro Data\' datestr(date, 'yymmdd'),'Data'];
exptag = [datestr(now,'yymmdd') '_' experiment_number];

mkdir(rootdir); % make root directory
cd(rootdir);    

mkdir(rootdir,exptag); % make experiment directory
cd([rootdir,'\',exptag]);

D = dir(['voltage*.mat']);
if isempty(D)
    n = 0 ;
else 
    n = length(D) ;
end   
data.n=n;

data.voltagefile = ['voltage_' exptag '_' int2str(data.n)];
data.currentfile = ['current_' exptag '_' int2str(data.n)];
data.Ao0file = ['Ao0_' exptag '_' int2str(data.n)];
data.Ao1file = ['Ao1_' exptag '_' int2str(data.n)];
data.timefile = ['time_', exptag '_' int2str(data.n)];

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
s.addAnalogInputChannel( 'Dev1','ai7','voltage');
s.addAnalogInputChannel( 'Dev1','ai8','voltage');

s.queueOutputData(outputData) ;
s.NotifyWhenDataAvailableExceeds = s.Rate*numsecs ; 

lh1 = s.addlistener('DataAvailable', @(s,event) savedataJC(data,event.TriggerTime,event.TimeStamps, event.Data));

s.IsContinuous = true ;
s.startBackground() ; 




