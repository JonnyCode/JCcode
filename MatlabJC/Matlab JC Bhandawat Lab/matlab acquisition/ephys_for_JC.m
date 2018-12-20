%% ===== Initializes cameras and data acquisition ============================
global data
close all; clear all; daqreset;imaqreset;
rootdir = ['C:\Desktop\' datestr(date, 'yymmdd')];
vidprefix = [datestr(date,'yymmdd') '_triggertest2'];
datfilename = [vidprefix '_input.mat'];

numsecs =5;
framerate = 95;
numframes = numsecs*framerate; %floor(100*20); %75*120;
experiment_number =num2str(1);

exptag = [datestr(now,'yymmdd') '_' experiment_number];
data.exptag = exptag;
data.duration = numsecs;
data.framerate=framerate;

try,
    cd(rootdir);
catch,
    mkdir(rootdir);
    cd(rootdir);    
end;
try,
    delete(datfilename);
catch,
end;

D = dir(['voltage*.mat']);
if isempty(D)
    n = 0
else 
    n = length(D)
end   
data.n=n;

s1 = 'video1';
s2=int2str(data.n);
s3='.avi';
data.avifile1 = [data.exptag s1 '_' s2 s3];
s4 = 'video2';
data.avifile2 = [data.exptag s4 '_' s2 s3];
data.voltagefile = ['voltage_' data.exptag '_' int2str(data.n)];
data.currentfile = ['current_' data.exptag '_' int2str(data.n)];
data.voltage=0;
data.current=0;
data.odorpulse=0;
s = daq.createSession('ni');
s.Rate=10000;

s.addAnalogOutputChannel('Dev1','ao0','Voltage');
s.addAnalogOutputChannel('Dev1','ao1','Voltage');

s.addAnalogInputChannel( 'Dev1','ai0','voltage');
s.addAnalogInputChannel( 'Dev1','ai1','voltage');
s.addAnalogInputChannel( 'Dev1','ai2','voltage');
s.addAnalogInputChannel( 'Dev1','ai3','voltage');
s.addAnalogInputChannel( 'Dev1','ai4','voltage');
s.addAnalogInputChannel( 'Dev1','ai5','voltage');
s.addAnalogInputChannel( 'Dev1','ai6','voltage');

t = 0:(1/s.Rate):(numsecs+1);
tempOut = ones(1,length(t)) ;
outputData(:,1)= ones(1,length(t));
outputData(:,2)= onses(;
s.queueOutputData(outputData);
s.NotifyWhenDataAvailableExceeds =s.Rate; %*10000;


%%
lh1 = s.addlistener('DataAvailable', @(s,event) savedata(data,event.TriggerTime,event.TimeStamps, event.Data));
%sh = s.addlistener('DataAvailable',@(s,event) savedata());
%lh1 = s.addlistener('DataAvailable', @(s,event) saveData_v2(datfilename,event.TriggerTime,event.TimeStamps, event.Data));

s.IsContinuous = true;
s.startBackground(); 
display('before waittime');
tic;
pause(10);
%wait([vid2, vid1],Inf);
toc
[frames_v1, timeStamp_v1, metdata1] = getdata(vid1);
[frames_v2, timeStamp_v2, metdata2] = getdata(vid2);
%end


% 



delete(lh1);
clear all;

