
% this script will plot out spike times detected and a give PSTH

% parameters
CellInfo_str = '111607Bc1_a' ;
epochsUwant_str ={'[82:89]','[92:101]','[102:111]'} ; % choose the set of epochs you want
epochCond_num = [8,9,10] ; % epochcondtion numbers for each epoch set
pre= 0  ;    % don't look for spike times before sample point 
post = 0 ;  % don't look for spikes after sample point
threshold =10 ; % spike threshold
color_str = {'b.','r.','g.','k.'} ; % color order of plotting
colorLine_str = {'b','r','g','k'} ; % color order of plotting

% Get appropriate cell structure format
cd ~/data_analysis/Index ;
load (CellInfo_str) ;  % parathesies allow it to load that varaible name

CellInfo = LoadSCIData(CellInfo,1) ;  % for two amps used 

for a = 1:length(epochsUwant_str) ; % for each set of epochs you want

epochsUwant = str2num(epochsUwant_str{a}) ;    % change str in numbers    
    
% this function excludes all the epochs but the ones you want
% Input = CellInfo, EpochConditionNumber, vector of the epochs you want
% Output = CellInfo, b = index of epochs in cellinfo
CellInfo = EpochExcluder(CellInfo,epochCond_num(a),epochsUwant) ;         

% this must come after EpochExcluder
EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ;

% This is Gabe's code to detect spike times from Cell Attached data
% Input = CellInfo file name, EpochCondtion(number), pre and post sample
% points to be ingored, spike threshold for detection
% Output = Cell of vector of spike times for each epoch
[SpikeTimeIndex{a}]=GetSpikeTimes_CellAttachedjc(CellInfo,EpochCondition(epochCond_num(a)), pre, post, threshold) ;

% Gabe's code to get PSTH and spiketrains
[SpikeTrain{a},PSTH{a}] = ConstructPSTH_JC(SpikeTimeIndex{a},length(CellInfo.EpochData.Data{epochsUwant(a)+1})) ;


% set up a rastor plot
for b=1:length(epochsUwant) ; % for each epoch you want to plot 
rastor{a}(b,:) = SpikeTrain{a}(b,:)*epochsUwant(b) ; % make the spike train rastor reflect the epoch number
end
c = find(rastor{a} == 0 ) ;    % make no spikes NaN instead of zero
rastor{a}(c) = NaN ;
clear b c


figure(1)
plot([1:length(rastor{a})],rastor{a}', color_str{a})
hold on

figure(2)
plot(smooth(PSTH{a},50),colorLine_str{a}) ;
hold on


% plot detected spikes ontop of real spikes
for b= 3:length(epochsUwant)+2 
% set(H,'Position',[-1590,549,2918,509]) ; % for 2 screens
figure
plot(CellInfo.EpochData.Data{epochsUwant(b-2)+1},'k') ;
hold on, plot([zeros(1,pre),SpikeTrain{a}(b-2,:)]*20,'r') ; % zeros is for graphing only, spike times will still read starting from pre time
A = gca ;
set(A,'xlim',([pre, (length(CellInfo.EpochData.Data{epochsUwant(b-2)+1})-post)])) ; 
end


end