clear all; 
daqreset; 

s = daq.createSession('ni');
s.Rate= 10000;

% set daq channels
s.addAnalogOutputChannel('Dev1','ao0','Voltage');
s.addAnalogInputChannel( 'Dev1','ai0','voltage');

outputData = zeros(s.Rate,1) ;
outputData(500:1000) = 1 ;

s.queueOutputData(outputData) ;
s.NotifyWhenDataAvailableExceeds = s.Rate ; 

lh1 = s.addlistener('DataAvailable', @Checker);
s.startBackground() ; 
