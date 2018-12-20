%% ===== Initializes cameras and data acquisition ============================

clear all; 
daqreset; 
global fromEphys

% parameters
experiment_number = '1' ;
params.VextStimType = 'singleVextStep' ;
params.OdorValveStimType = 'singleValveStep' ;

params.numsecs = 5 ; % for each trial
params.numRepeats = 10 ; % number of repeated trials
params.TimeBetweenRepeats = 60 ; % (sec)
params.sampleRate = 10000 ; % hz
params.InjPulseInitTime = .5 ; % (sec) start time of Rinput test 
params.InjPulseDur = .5 ; % (sec) duration of Rinput test
params.InjPulseAmp = -0.05 ; % (volts) from daq to external clamp signal 
params.OdorPulseInitTime = 3 ; % (sec) start time of odor pulse
params.OdorPulseDur = .5 ; % (sec) duration of odor pulse
params.OdorPulseAmp = 10 ; % (volts) from daq to air valve
params.updateFig = false ; % true if you want to update figures every second (false otherwise)

params.BinaryTime = .2 ; % (sec)

% set daq output
fromEphys.lt = params.numsecs*params.sampleRate;
outputData = zeros(fromEphys.lt,2) ; % ao0

% odor
if strcmp(params.OdorValveStimType,'singleValveStep') ;
    outputData(params.OdorPulseInitTime*params.sampleRate:(params.OdorPulseInitTime+params.OdorPulseDur)*params.sampleRate,1)= params.OdorPulseAmp ;
end
    
if strcmp(params.OdorValveStimType,'ValveNoise') ;
    tempV = zeros(1,9) ;
    Mmat = [] ;
    
    for a = 1:9 ;
        tempV(1:a)=1 ;
        tempVperms = unique(perms(tempV),'rows') ;
        Mmat = [Mmat; tempVperms] ;
    end
    PermVec = randperm(2^9-1) ;
    MmatRand = Mmat(PermVec,:) ;
    Mvec = reshape(MmatRand',1,numel(Mmat)) ;
    
    BPnts = params.BinaryTime*params.sampleRate ;
    
    outputDataTemp = zeros(length(Mvec),1) ;
    for a=1:length(Mvec) ;
        outputDataTemp((a-1)*BPnts+1:a*BPnts,1)= params.OdorPulseAmp*Mvec(a) ;
    end
    outputData(:,1) = outputDataTemp(1:fromEphys.lt,1) ;
    outputData(end-2:end,1) = 0 ;
end

% Vext
if strcmp(params.VextStimType,'singleVextStep') ;
    outputData(params.InjPulseInitTime*params.sampleRate:(params.InjPulseInitTime+params.InjPulseDur)*params.sampleRate,2)= params.InjPulseAmp ; %step
elseif strcmp(params.StimType,'VextNoise') ;
    outputData(:,2) = normrnd(0,params.InjPulseAmp,1,fromEphys.lt) ; % noise injection
    outputData(:,2) = lowPassFilter(outputData(:,2)',params.sampleRate,60) ;
    outputData(end-10:end,2) = 0 ;
end

% prep data file  
rootdir = ['C:\Cafaro Data\' datestr(date, 'yymmdd'),'Data'];
exptag = [datestr(now,'yymmdd') '_' experiment_number];
notestag = [datestr(now,'yymmdd'), '_Notes.txt'] ;

mkdir(rootdir); % make root directory 
mkdir(rootdir,exptag); % make experiment directory
cd([rootdir,'\',exptag]);

% note data block in notes
fromEphys.NotesFile = [rootdir,'\',notestag] ;

NoteTime = datestr(now) ;
fieldnamesParams = fieldnames(params) ;

fid = fopen(fromEphys.NotesFile,'a');
for a=1:length(fieldnamesParams) ;
    if a==1 ;
        Note = [datestr(now),' data block'] ;
        fprintf(fid,'%s \r\n',Note);
    end
    
    if ~ischar(params.(fieldnamesParams{a})) ;
        Note = [fieldnamesParams{a},'=',num2str(params.(fieldnamesParams{a}))] ;
    else
        Note = [fieldnamesParams{a},'=',params.(fieldnamesParams{a})] ;
    end
    fprintf(fid,'%s \r\n',Note);
    
end
fclose(fid);

% run data block
for a=1:params.numRepeats ;

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
    s.Rate= params.sampleRate;

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
    
    if params.updateFig ;
        s.NotifyWhenDataAvailableExceeds = params.sampleRate ; 
    else
        s.NotifyWhenDataAvailableExceeds = fromEphys.lt ;
    end
        
    lh1 = s.addlistener('DataAvailable', @savedataJCv3);

    s.startBackground() ; 

    s.wait() ;
    
    clear lhl s
    if a~=params.numRepeats ;
        pause(params.TimeBetweenRepeats) ;
    end
end

