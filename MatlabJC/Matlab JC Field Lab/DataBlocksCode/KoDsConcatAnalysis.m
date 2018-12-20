function ForIgor = KoDsConcatAnalysis(DataBlock, DB, Params)

% this function will compare DS cell DG responses from a concatinated datablock
% JC 2017-12-04

RunAsScript = true ;
if RunAsScript ;
    DB = 45 ;
    [DataBlock,Params] = DataBlocks_KO ;
end

DsPathNums = [1,2] ; % path numbers to compare
ConcatPathNum = 1 ;
numElectrodeLayers = 2 ;

% parameters
Color_list = {'k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y'} ; % order of colors for each 

% save path
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/KoCompare/DsIdsDb',num2str(DB),'ConcatPathNum',num2str(ConcatPathNum)] ;
 
% load data
dataRun = load_data(DataBlock(DB).DsPathConcat{ConcatPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION
dataRun = load_ei(dataRun, 'all') ;

% get ei center of mass 
for cells = 1:length(dataRun.spikes) ; % for each cell
    EiCtr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end

for DsSet = 1:length(DsPathNums) ; % for each DsSet 
    dataRunDg{DsSet} = load_data(DataBlock(DB).DsPath{DsPathNums(DsSet)}) ;
    dataRunDg{DsSet} = load_neurons(dataRunDg{DsSet}) ;
end

% load stim
TrialTrigInterval = 10 ;
for DsSet = 1:length(DsPathNums) ; % for each DsSet 
    slashi = strfind(DataBlock(DB).DsPath{DsPathNums(DsSet)},'/') ; % find the /
    dataRunDg{DsSet}.names.stimulus_path = [DataBlock(DB).DsPath{DsPathNums(DsSet)}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).DsPath{DsPathNums(DsSet)}(end-1:end),'.txt'] ;
    dataRunDg{DsSet} = load_stim(dataRunDg{DsSet},'user_defined_trigger_interval', TrialTrigInterval) ;
end

% stim trials
for DsSet = 1:length(DsPathNums) ; % for each DsSet 
    for t=1:length(dataRunDg{DsSet}.stimulus.trials) ; % for each DG trial
        DgParams{DsSet}(t,1) = dataRunDg{DsSet}.stimulus.trials(t).DIRECTION ; % direction
        DgParams{DsSet}(t,2) = dataRunDg{DsSet}.stimulus.trials(t).TEMPORAL_PERIOD; % temporal period
    end
end

% stim parameters (assume params are all the same as DsSet 1)
NumDirDg = length(dataRunDg{1}.stimulus.params.DIRECTION) ; % number of directions
NumTempPeriods = length(dataRunDg{1}.stimulus.params.TEMPORAL_PERIOD) ; % number of directions

% triggers
trigs = [] ;
for DsSet = 1:length(DsPathNums) ; % for each DsSet
    trigs = [trigs, dataRunDg{1}.duration*(DsSet-1) + dataRunDg{DsSet}.stimulus.triggers] ;
end
    
% spike times organized by stim direction and speed

NumCells = length(dataRun.spikes) ;
for DsSet = 1:length(DsPathNums) ; % for each DsSet 
    for cells=1:NumCells ; % for each cell
        SpikeNum_Mean{DsSet}{cells} = nan(NumTempPeriods,NumDirDg) ;
        SpikeNum_Std{DsSet}{cells} = nan(NumTempPeriods,NumDirDg) ;
        for st = 1:NumDirDg ; % for each direction 
            for tmp = 1:NumTempPeriods ; % for each temp period
                prm = [dataRunDg{DsSet}.stimulus.params.DIRECTION(st),...
                    dataRunDg{DsSet}.stimulus.params.TEMPORAL_PERIOD(tmp)] ; % params set
                ti = find((DgParams{DsSet}(:,1)==prm(1)).*(DgParams{DsSet}(:,2)==prm(2))==1) ; % index of triggers for that stim
                ti = ti+length(dataRunDg{1}.stimulus.triggers)*(DsSet-1) ; % adjust for previous trigs
                
                for t=1:length(ti) ; % for each trial
                    spk = dataRun.spikes{cells}-trigs(ti(t)) ;
                    spk = spk(spk>=0 & spk<=TrialTrigInterval) ;
                    SpikeTimes{DsSet}{cells}{st}{tmp}{t} = spk ; 
                    SpikeNum{DsSet}{cells}{tmp}(st,t) = length(spk) ;
                end
                SpikeNum_Mean{DsSet}{cells}(tmp,st) = mean(SpikeNum{DsSet}{cells}{tmp}(st,:)) ; % average across trials
                SpikeNum_Std{DsSet}{cells}(tmp,st) = std(SpikeNum{DsSet}{cells}{tmp}(st,:)) ; % std across trials
            end
        end
    end
end

% ei centers

%% parameters on putative DS

PutativeDsi = [5,11,16:18,24,25,40,42,44,47,48,56,61,70,73,74,87,105,126,130,131,133:135,140,144,...
        148,149,152,157,160:162,165,166,177,179,181,187,191,194,195,203,208,214,215,228,230,237,239,...
        242,251,263,269,272,277,282,287,291,294,299,312,322,324,327,339,340,353,355:357,359,367,373,...
        375,379,383,388]; % DS cell selected qualititatively
  
ArrayImagePath = '/Volumes/lab/Experiments/Array/Images/2017-11-30-0/10x td-tom 8.jpg' ;         
Tform = map_Ei_to_camera(ArrayImagePath, DataBlock(DB).DsPathConcat{ConcatPathNum}) ;


%% figures

figure ; % tuning curves - all cells
for cells=PutativeDsi ;% 1:NumCells ; % for each cell
    clf
    for tmp = 1:NumTempPeriods ; % for each temp period
        subplot(NumTempPeriods,1,tmp) 
        for DsSet = 1:length(DsPathNums) ; % for each DsSet 
            errorbar(dataRunDg{1}.stimulus.params.DIRECTION,SpikeNum_Mean{DsSet}{cells}(tmp,:),SpikeNum_Std{DsSet}{cells}(tmp,:))
            hold on
        end
        xlabel('direction (deg)')
        ylabel('spike number')
        title(num2str(cells))
    end
    pause
end

figure ; % min and max response compare
for cells=PutativeDsi ;% 1:NumCells ; % for each cell

    for tmp = 1:NumTempPeriods ; % for each temp period
        [minCntrl,mini] = min(SpikeNum_Mean{1}{cells}(tmp,:)) ; % min of control
        minDrug = SpikeNum_Mean{2}{cells}(tmp,mini) ; % min of +Drug
        
        [maxCntrl,maxi] = max(SpikeNum_Mean{1}{cells}(tmp,:)) ; % max of control
        maxDrug = SpikeNum_Mean{2}{cells}(tmp,maxi) ; % min of +Drug
        
        subplot(NumTempPeriods,2,(tmp-1)*2+1)
        plot(minCntrl/maxCntrl,minDrug/maxCntrl,'k*')
        hold on
        xlabel('cntrl min / cntrl max')
        ylabel('drug at cntrl min / cntrl max')
        
        subplot(NumTempPeriods,2,(tmp-1)*2+2)
        plot(maxCntrl/maxCntrl,maxDrug/maxCntrl,'k*')
        hold on
        xlabel('cntrl max / cntrl max')
        ylabel('drug at cntrl max / cntrl max')
    end
end

subplot(NumTempPeriods,2,1)
%plot([.1,80],[.1,80],'r') 
plot([.001,10],[.001,10],'r') 
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(NumTempPeriods,2,2)
%plot([10,200],[10,200],'r')
plot([.001,10],[.001,10],'r') 
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(NumTempPeriods,2,3)
%plot([.1,80],[.1,80],'r')
plot([.001,10],[.001,10],'r') 
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(NumTempPeriods,2,4)
%plot([10,200],[10,200],'r')
plot([.001,10],[.001,10],'r') 
set(gca,'xscale','log')
set(gca,'yscale','log')


h = figure ; % rasters
NumTrials = length(SpikeTimes{1}{1}{1}{1}) ; 
for cells= [5,11,16:18,24,25,40,42,44,47,48,56,61,70,73,74,87,105,126,130,131,133:135,140,144,...
        148,149,152,157,160:162,165,166,177,179,181,187,191,194,195,203,208,214,215,228,230,237,239,...
        242,251,263,269,272,277,282,287,291,294,299,312,322,324,327,339,340,353,355:357,359,367,373,...
        375,379,383,388]; % DS cell i
    clf
    Rnd = 1 ;
    for st = 1:NumDirDg ; % for each direction 
        for tmp = 1:NumTempPeriods ; % for each temp period
            subplot(NumDirDg,NumTempPeriods,Rnd)
            for DsSet = 1:length(DsPathNums) ; % for each DsSet
                for t=1:NumTrials ; % for each trial
                    SpkY = ones(1,length(SpikeTimes{DsSet}{cells}{st}{tmp}{t}))*(t+(DsSet-1)*NumTrials) ;
                    for spk = 1:length(SpkY) ; % for each spike
                        plot([1,1]*SpikeTimes{DsSet}{cells}{st}{tmp}{t}(spk),[SpkY(spk),SpkY(spk)+1],Color_list{DsSet})
                        hold on
                    end
                end 
            end
            axis([0,10,1,NumTrials*length(DsPathNums)+1])
            Rnd = Rnd+1 ;
        end
    end
    title(num2str(cells))
    %hgsave(['celli',num2str(cells),'.mat'])
    %print(h, '-dtiff',['celli',num2str(cells),'.dtiff'])
    pause
end
            

figure ; % 
imshow(ArrayImagePath)
hold on
for cells= PutativeDsi ;
    ctrMon = tformfwd(Tform,EiCtr(cells,1),EiCtr(cells,2)) ;
    plot(ctrMon(1),ctrMon(2),'ro')
end


