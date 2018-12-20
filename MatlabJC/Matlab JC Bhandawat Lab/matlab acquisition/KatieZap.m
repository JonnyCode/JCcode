%params 
freq = 10 ; % hz
amp = 5 ; % volts
offset = 1 ; % volts
duration = 10 ; % sec
sampleRate = 10000 ; % Hz

% analog output
t=[0:1/sampleRate:duration] ;
OutputWave = amp*square(sin(t*2*pi*freq))+offset ;

% daq
s = daq.createSession('ni');
s.Rate= sampleRate;
s.addAnalogInputChannel( 'Dev1','ai0','voltage');
s.addAnalogOutputChannel('Dev1','ao0','Voltage');
s.queueOutputData(OutputWave');

data = s.startForeground() ; 