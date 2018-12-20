%% ===== Initializes cameras and data acquisition ============================

clear all; 
daqreset; 
global fromEphys

% parameters
experiment_number = '1' ;
params.VextStimType = 'BackgroundOdorAlternate' ;
params.OdorValveStimType = 'singleValveStep' ;

params.numsecs = 30 ; % for each trial
params.numRepeats = 1 ; % number of repeated trials
params.TimeBetweenRepeats = 60 ; % (sec)
params.sampleRate = 10000 ; % hz
params.InjPulseInitTime = 5 ; % (sec) start time of Rinput test 
params.InjPulseDur = 15 ; % (sec) duration of Rinput test
params.InjPulseAmp = 10 ; % -0.05default (volts) from daq to external clamp signal 
params.InhPulseStep = .05 ; % (volts) size of step taken on each trial from "InjPulseAmp"
params.OdorPulseInitTime = 15 ; % (sec) start time of odor pulse
params.OdorPulseDur = .5 ; % (sec) duration of odor pulse
params.OdorPulseAmp = 10 ; % (volts) from daq to air valve
params.updateFig = false ; % true if you want to update figures every second (false otherwise)

params.BinaryTime = .2 ; % (sec)

% from params
fromEphys.lt = params.numsecs*params.sampleRate ;

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

% prep output data
for a=1:params.numRepeats ;
    
    % set daq output
    outputData{a} = zeros(fromEphys.lt,2) ; % ao0

    % odor
    if strcmp(params.OdorValveStimType,'singleValveStep') ;
        outputData{a}(params.OdorPulseInitTime*params.sampleRate:(params.OdorPulseInitTime+params.OdorPulseDur)*params.sampleRate,1)= params.OdorPulseAmp ;
    end

    if strcmp(params.OdorValveStimType,'ValveNoise') ;
        tempV = zeros(1,9) ;
        Mmat = [] ;

        for b = 1:9 ;
            tempV(1:b)=1 ;
            tempVperms = unique(perms(tempV),'rows') ;
            Mmat = [Mmat; tempVperms] ;
        end
        PermVec = randperm(2^9-1) ;
        MmatRand = Mmat(PermVec,:) ;
        Mvec = reshape(MmatRand',1,numel(Mmat)) ;

        BPnts = params.BinaryTime*params.sampleRate ;

        outputDataTemp = zeros(length(Mvec),1) ;
        for b=1:length(Mvec) ;
            outputDataTemp((b-1)*BPnts+1:b*BPnts,1)= params.OdorPulseAmp*Mvec(b) ;
        end
        outputData{a}(:,1) = outputDataTemp(1:fromEphys.lt,1) ;
        outputData{a}(end-2:end,1) = 0 ;
    end

    % Vext
    if strcmp(params.VextStimType,'singleVextStep') ;
        outputData{a}(params.InjPulseInitTime*params.sampleRate:(params.InjPulseInitTime+params.InjPulseDur)*params.sampleRate,2)= params.InjPulseAmp ; %step
    end
    
    if strcmp(params.VextStimType,'VextNoise') ;
        outputData{a}(:,2) = normrnd(0,params.InjPulseAmp,1,fromEphys.lt) ; % noise injection
        outputData{a}(:,2) = lowPassFilter(outputData(:,2)',params.sampleRate,60) ;
        outputData{a}(end-10:end,2) = 0 ;
    end
    
    if strcmp(params.VextStimType,'VextStepFam') ;
        outputData{a}(params.InjPulseInitTime*params.sampleRate:(params.InjPulseInitTime+params.InjPulseDur)*params.sampleRate,2)= params.InjPulseAmp+(a-1)*params.InhPulseStep ; %step family
    end
    
    if strcmp(params.VextStimType,'BackgroundOdorAlternate') ;
         outputData{a}(params.InjPulseInitTime*params.sampleRate:(params.InjPulseInitTime+params.InjPulseDur)*params.sampleRate,2)= params.InjPulseAmp*rem(a+1,2) ; %step off for odd, on for even
    end
    
end
    

% run data block
for a=1:params.numRepeats ;    
    
    % trial number
    D = dir(['voltage*.mat']);
    if isempty(D)
        n = 0 ;
    else 
        n = length(D) ;
    end   

    % global variables 
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

    s.queueOutputData(outputData{a}) ;
    
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

