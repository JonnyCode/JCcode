classdef ephysJCoo < handle
    properties
        inputParamsDir = '/Users/Shared/SpatialStimParams/'        
        outputParamsDir = '~/acquisition/SpatialStimParams/'
        inputParamsFile = 'TempParams.h5'
        %savedParamsDir = '/Users/Greg/acquisition/'
        %savedParamsDir = '/Users/Shared/SpatialStimDataFiles/'
        MatlabCodeDir = '~/matlab/common/stim/'
        testParamsDir = ''
        tempFile = 'TempParams.h5'
        repositoryVersion
        frameRate = 60; %frames / sec
        flipInterval; %actual acurate flip interval (s)
        windPtr
        imageArray = [];
        movieFrames
        gammaTable
        lastEpoch = -1;
    end
    
    methods
        
        function self = ephysJCoo(varargin)
            
             if strcmp(varargin{1}, 'WriteDefaults')
                AnalogOutputClass = varargin{2};
                props = properties(AnalogOutputClass); 
                params = struct;
                stringParams = {};
                for i=1:length(props)
                    if strmatch(props{i}, stringParams)
                        params.(props{i}) = ' ';
                    else
                        params.(props{i}) = 0; %maybe read deafult from file
                    end
                end
                params.AnalogOutputClass = AnalogOutputClass ;
             end
             
            if strcmp(varargin{1}, 'Run') % = start button
                ProtocolClass = varargin{2}; % protocol
                ProtocolClass.PrepProtocol ; % prep protocol output
                ProtocolClass.RunProtocol ; % run protocol
            end
                
             end   
        end
    end
end
            

    

% clear all ; 
% daqreset; 
% global fromEphys
% 
% % parameters
% experiment_number = '2' ;
% params.VextStimType = 'singleVextStep' ;
% params.OdorValveStimType = 'singleValveStep' ;
% fromEphys.TrialTag{1} = [] ;
% 
% 
% 
% % from params
% fromEphys.lt = params.numsecs*params.sampleRate ;
% 
% % prep data file  
% fromEphys.rootdir = ['C:\Cafaro Data\' datestr(date, 'yymmdd'),'Data'];
% exptag = [datestr(now,'yymmdd') '_' experiment_number];
% notestag = [datestr(now,'yymmdd'), '_Notes.txt'] ;
% 
% mkdir(fromEphys.rootdir); % make root directory 
% mkdir(fromEphys.rootdir,exptag); % make experiment directory
% cd([fromEphys.rootdir,'\',exptag]);
% 
% % note data block in notes
% fromEphys.NotesFile = [fromEphys.rootdir,'\',notestag] ;
% 
% NoteTime = datestr(now) ;
% fieldnamesParams = fieldnames(params) ;
% 
% fid = fopen(fromEphys.NotesFile,'a');
% for a=1:length(fieldnamesParams) ;
%     if a==1 ;
%         Note = [datestr(now),' data block'] ;
%         fprintf(fid,'%s \r\n',Note);
%     end
%     
%     if ~ischar(params.(fieldnamesParams{a})) ;
%         Note = [fieldnamesParams{a},'=',num2str(params.(fieldnamesParams{a}))] ;
%     else
%         Note = [fieldnamesParams{a},'=',params.(fieldnamesParams{a})] ;
%     end
%     fprintf(fid,'%s \r\n',Note);
%     
% end
% fclose(fid);
% 
% % prep output data
% for a=1:params.numRepeats ;
%     
%   
% end
%     
% 
% % run data block
% for a=1:params.numRepeats ;    
%     
%     % trial number
%     D = dir('*.mat') ;
%     if isempty(D)
%         n = 0 ;
%     else 
%         n = length(D) ;
%     end   
% 
%     % global variables 
%     fromEphys.n = n ;
%     fromEphys.TrialId = ['data_' exptag '_' int2str(n)];
%     
%     % daq object for session
%     s = daq.createSession('ni');
%     s.Rate= params.sampleRate;
% 
%     % set daq channels
%     s.addAnalogOutputChannel('Dev1','ao0','Voltage');
%     s.addAnalogOutputChannel('Dev1','ao1','Voltage');
% 
%     s.addAnalogInputChannel( 'Dev1','ai0','voltage');
%     s.addAnalogInputChannel( 'Dev1','ai1','voltage');
%     s.addAnalogInputChannel( 'Dev1','ai2','voltage');
%     s.addAnalogInputChannel( 'Dev1','ai3','voltage');
%     s.addAnalogInputChannel( 'Dev1','ai4','voltage');
%     s.addAnalogInputChannel( 'Dev1','ai5','voltage');
%     s.addAnalogInputChannel( 'Dev1','ai6','voltage');
% 
%     s.queueOutputData(outputData{a}) ;
%     
%     if params.updateFig ;
%         s.NotifyWhenDataAvailableExceeds = params.sampleRate ; 
%     else
%         s.NotifyWhenDataAvailableExceeds = fromEphys.lt ;
%     end
%         
%     lh1 = s.addlistener('DataAvailable', @savedataJCv4);
% 
%     s.startBackground() ; 
% 
%     s.wait() ;
%     
%     clear lhl s
%     if a~=params.numRepeats ;
%         pause(params.TimeBetweenRepeats) ;
%     end
% end
