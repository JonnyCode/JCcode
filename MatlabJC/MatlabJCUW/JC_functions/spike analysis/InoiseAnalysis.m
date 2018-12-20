% This script will get the spike triggered current from current noise
% injection

% Get the following outputs (1= get output, 0= don't) 
% Prespike stimulus (STI)
perform.INa = 1 ;         
perform.INb = 1 ;    
perform.INc = 1 ;    
perform.INd = 1 ;    
perform.INe = 0 ;    
perform.INf = 0 ;    
perform.INg = 0 ;    

% prep matrix if it doesn't exist already
if ~exist('ForIgor')
    ForIgor=struct() ;
end

% using the following cells...
for A = [2] ;
%% Cells for analysis

Input(1).cellname = '092408Bc1b' ; % good patch from naked retina with no light response (previously in blockers but in wash during these recordings)
Input(1).INa = '[71:115]' ; % mean = 200, sd=400          (using letters to denote different current injection parameters)
Input(1).INb = '[46:70]' ; % mean=200, sd=200
Input(1).INc = '[6:10,31:45]' ; % mean=200, sd=100
Input(1).INd = '[11:30]' ; % mean=200, sd=50

Input(2).cellname = '092408Bc3b' ; % good patch from naked retina with no light response (previously in blockers but in wash during these recordings)
Input(2).INa = '[46:65]' ; % mean = 200, sd=400
Input(2).INb = '[26:45]' ; % mean=200, sd=200
Input(2).INc = '[0:24]' ; % mean=200, sd=100
Input(2).INd = '[66:105]' ; % mean=200, sd=50

Input(3).cellname = '092408Ac1' ; % Fred's cell from naked retina with no light response (in blockers during recoding)
Input(3).INe = '[5:19]' ; % mean=0 sd=10
Input(3).INf = '[20:39]' ; %mean=10 sd=5
Input(3).INg = '[43:98]' ; %mean=13 sd=8

% Parameters assuming sample rate at 10 kHz
Parameters.PrePnts = 10000 ;    % number of points before stimuli starts
Parameters.StmPnts = 60000 ;    % number of points during stimuli
Parameters.PostPnts = 1000 ;    % number of points after stimuli ends
Parameters.STAPnts = 3000 ;     % number of points used for prespike wave forms
Parameters.DecimatePnts = 10 ;  % number of points averaged together in order to downsample prespike waveforms
Parameters.SmoothPnts = 100 ;   % number of points used to smooth spike train to make PSTH
Parameters.QuietPnts = 200 ;    % number of points used to identify quiet time before burst
Parameters.MinSpikeNumber = 2 ; % minimum number of spikes in a burst

%%

if perform.INa == 1 ; 
    id = 'INa' ;
    Temp = SpkTrigCurr(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.INb == 1 ; 
    id = 'INb' ;
    Temp = SpkTrigCurr(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.INc == 1 ; 
    id = 'INc' ;
    Temp = SpkTrigCurr(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.INd == 1 ; 
    id = 'INd' ;
    Temp = SpkTrigCurr(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.INe == 1 ; 
    id = 'INe' ;
    Temp = SpkTrigCurr(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.INf == 1 ; 
    id = 'INf' ;
    Temp = SpkTrigCurr(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

if perform.INg == 1 ; 
    id = 'INg' ;
    Temp = SpkTrigCurr(Input,Parameters,id,A) ;
    ForIgor = mergeStruct(ForIgor,Temp) ;
end

end


