function ForIgor = ImSliderConcatAnalysis(DataBlock, DB, Params)


% this function is adapted from 'ImSliderAnalysis' to deal with
% concatinated data sets, instead of mapped data sets.

% JC 5/15/2017 

LoadOldMatFileFlag = false; % true skip calculations and just load
PlotFigsFlag = true ;
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/'] ;
numElectrodeLayers = 2 ; % number of electrode layers surrounding the max that the center is calculated over
eiDistBin = 100 ; % (um) bins 
eiDistMax = 1000 ; % (um) max distance ei can be separated

% load vaiables if saved
% if LoadOldMatFileFlag ;
%     load([saveFigPath,'entire'])
%     PlotFigsFlag = true ;
% else

ConcatPathNum = 3 ; %TEMP -SHOULD BE INPUT
ImPathNum = 3 ;
DsPathNum = 2 ;


% parameters
Color_list = {'k','r','b','g','c','y','k','r','b','g','c','y'} ; % order of colors for each 

saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/DsIdsDb',num2str(DB),'ConcatPathNum',num2str(ConcatPathNum)] ;

    
% load stimulus 
slashi = strfind(DataBlock(DB).ImSlide{ImPathNum},'/') ; % find the /
StimPath = [DataBlock(DB).ImSlide{ImPathNum}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).ImSlide{ImPathNum}(end-1:end)] ;
load(StimPath) ;

% load data
dataRun = load_data(DataBlock(DB).ImSlideConcat{ConcatPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;

% ImSlide stimulus time
dataRunIm = load_data(DataBlock(DB).ImSlide{ConcatPathNum}) ;
dataRunIm = load_neurons(dataRunIm) ;
ImStimTime = dataRunIm.duration ;
DgStimTime = dataRun.duration-ImStimTime ;
clear dataRunIm ;

% get ei center of mass 
for cells = 1:length(dataRun.spikes) ; % for each cell
    EiCtr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end

% get triggers from 
triggs = dataRun.triggers(dataRun.triggers>DgStimTime) ; % image slide

% stim params
NumStimShort = size(stimulus.directionsShown,1)*size(stimulus.directionsShown,2)-length(triggs) ;
if NumStimShort~=0 ; % if not every stim was shown
    NumCompleteTrials = floor(length(triggs)/size(stimulus.directionsShown,2)) ;
    stimulus.directionsShown = stimulus.directionsShown(1:NumCompleteTrials,:) ; % take out the entire unfinished trials
    triggs = triggs(1:length(stimulus.directionsShown(:))) ;
end 
StimDuration = mean(diff(triggs)) ; % average stimulus duration
StimFrameRate = StimDuration/stimulus.num_frames ; % average frame rate sec/frame
dirShownTranspose = stimulus.directionsShown' ;
NumTrialsPerDir = size(stimulus.directionsShown,1) ; % number of trial for each direction

% id ds cells in DsPath
try load(saveDsIdsPath) ; % if they are ds ids already saved
catch % if not find them

    Params.TimeBounds = [0,DgStimTime] ;
    Params.DsPathNum = DsPathNum ;
    DataBlock(DB).DsConcatPath = DataBlock(DB).ImSlideConcat{ConcatPathNum} ;
    ForIgor = DsCellFinder(DataBlock, DB, Params) ;

    save(saveDsIdsPath, 'ForIgor')
end

% On-Off cell from DsPath
cell_id = ForIgor.ds_id{2} ; 
DsTypeName = ForIgor.dsName{2} ; 
dsi = ForIgor.dsi{2} ;

for DsType=1:length(DsTypeName) ;
    DsCelli{DsType} = get_cell_indices(dataRun,cell_id{DsType}) ;
end

% spike count and psth
for cells=1:length(dataRun.spikes) ; % for each cell
    for t=1:length(triggs) ;
        spk = dataRun.spikes{cells}-triggs(t) ;
        spk = spk(spk>=0 & spk<=StimDuration) ;
        spkCount(cells,t) = length(spk) ;
    end
end

% spike count for each image direction
for d = 1:length(stimulus.directions) ; % for each diection
   di = find(stimulus.directionsShown'==stimulus.directions(d)) ; % trials with that direction
   spkCount_direction_mean(:,d) = mean(spkCount(:,di),2) ;
   spkCount_direction_std(:,d) = std(spkCount(:,di),[],2) ;
end
 
% taylor/vaney dsi
for cells=1:length(dataRun.spikes) ; % for each cell
    PV = PolarVectorAddition([stimulus.directions',spkCount_direction_mean(cells,:)']) ;
    PreferedDirection(cells) = PV(1) ;
    DSindex(cells) = PV(2)/sum(spkCount_direction_mean(cells,:)) ;
end
    
% psth
psthTime = [0:StimFrameRate:StimDuration] ; 
for cells=1:length(dataRun.spikes) ; % for each cell
    for st = 1:length(stimulus.directions) ; % for each stimulus direction
        ti = find(dirShownTranspose(:)==stimulus.directions(st)) ; % index of triggers for that stim
        for t=1:NumTrialsPerDir ; % for each trial
            spk = dataRun.spikes{cells}-triggs(ti(t)) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;
            psth{cells}{st}(t,:) = histc(spk,psthTime) ;
        end
        psth_mean{cells}(st,:)=mean(psth{cells}{st},1) ;
        psth_var{cells}(st,:)=var(psth{cells}{st},[],1) ;
    end
end

%% noise correlations
% psth residuals
for cells=1:length(dataRun.spikes) ; % for each cell
    for st = 1:length(stimulus.directions) ; % for each stimulus direction
        for t=1:NumTrialsPerDir ; % for each trial
            psth_residual{cells}{st}(t,:) = psth{cells}{st}(t,:)-psth_mean{cells}(st,:) ;
        end
    end
end

% noise correlations as function of time (CAUTION- not local residuals)
for st = 1:length(stimulus.directions) ; % for each stimulus direction
    corr{st} = nan(length(stimulus.directions),length(stimulus.directions),length(psthTime)) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        for cells2=1:length(dataRun.spikes) ; % for each cell
            for pt = 1:length(psthTime) ; % for each time point of psth
                temp = corrcoef(psth_residual{cells}{st}(:,pt)',psth_residual{cells2}{st}(:,pt)') ;
                corr{st}(cells,cells2,pt)= temp(1,2) ;
            end
        end
    end
end
                
% corr histograms 
for st = 1:length(stimulus.directions) ; % for each stimulus direction                
    corrhist(st,:) = hist(corr{st}(:),[-1:.1:1]) ; 
end

%% drifting grating data

% stim path
TrialTrigInterval = 10 ;
slashi = strfind(DataBlock(DB).DsPath{DsPathNum},'/') ; % find the /
dataRunDg.triggers = dataRun.triggers(dataRun.triggers<=DgStimTime) ;
dataRunDg.names.stimulus_path = [DataBlock(DB).DsPath{Params.DsPathNum}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).DsPath{DsPathNum}(end-1:end),'.txt'] ;
dataRunDg = load_stim(dataRunDg,'user_defined_trigger_interval', TrialTrigInterval) ;

% psth
psthTimeDg = [0:StimFrameRate:TrialTrigInterval] ; 
for cells=1:length(dataRun.spikes) ; % for each cell
    for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
        ti = find(dataRunDg.stimulus.trial_list==st) ; % index of triggers for that stim
        for t=1:length(ti) ; % for each trial
            spk = dataRun.spikes{cells}-dataRunDg.stimulus.triggers(ti(t)) ;
            spk = spk(spk>=0 & spk<=TrialTrigInterval) ;
            psthDg{cells}{st}(t,:) = histc(spk,psthTimeDg) ;
        end
        psthDg_mean{cells}(st,:)=mean(psthDg{cells}{st},1) ;
        psthDg_var{cells}(st,:)=var(psthDg{cells}{st},[],1) ;
    end
end

%% figures

figure % plot tuning curves (all cells)
for cells=1:length(dataRun.spikes) ; % for each cell
    errorbar(stimulus.directions,spkCount_direction_mean(cells,:),spkCount_direction_std(cells,:))
    hold on
end
   
figure % plot tuning curves (ds cells only)
for DsType = 1:length(DsTypeName) ;
    subplot(length(DsTypeName),1,DsType)
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        errorbar(stimulus.directions,spkCount_direction_mean(DsCelli{DsType}(cells),:),spkCount_direction_std(DsCelli{DsType}(cells),:))
        hold on
    end
end
   
figure % dsi
plot(DSindex,'k*')
hold on
for DsType = 1:length(DsTypeName) ;
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        plot(DsCelli{DsType}(cells),DSindex(DsCelli{DsType}(cells)),'ro')
    end
end

figure % dsi and directions (polar)
polar([0,0],[0,1])
hold on
for cells=1:length(dataRun.spikes) ; % for each cell
    polar([PreferedDirection(cells),PreferedDirection(cells)]*pi/180,[DSindex(cells),DSindex(cells)],'*')
    hold on
end

figure % dsi and direction (linear)
plot(DSindex,PreferedDirection,'*')

figure % plot raster (ordered by stimulus direction)
NumTrialsPerDir = size(stimulus.directionsShown,1) ; % number of trial for each direction
for cells=1:length(dataRun.spikes) ; % for each cell
    temp = stimulus.directionsShown' ;
    [temp,ti] = sort(temp(:)) ; % get direction
    for t=1:length(triggs) ;
        spk = dataRun.spikes{cells}-triggs(ti(t)) ;
        spk = spk(spk>=0 & spk<=StimDuration) ;
        plot(spk,ones(1,length(spk))*t,'k.')
        axis([0,4,0,length(triggs)]) 
        hold on
    end
    title(num2str(cells))
    hold off
    pause
end

figure % plot raster (ordered by stimulus direction) - DS cells only
NumTrialsPerDir = size(stimulus.directionsShown,1) ; % number of trial for each direction
for DsType = 1:length(DsTypeName) ;
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        temp = stimulus.directionsShown' ;
        [temp,ti] = sort(temp(:)) ; % get direction
        for t=1:length(triggs) ;
            spk = dataRun.spikes{DsCelli{DsType}(cells)}-triggs(ti(t)) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;
            plot(spk,ones(1,length(spk))*t,'ro')
            axis([0,4,0,length(triggs)]) 
            hold on
        end
        title(num2str(DsCelli{DsType}(cells)))
        
        pause
        clf
    end
end


figure % plot raster of all cells for each trial - DS cells only
for t=1:length(triggs) ;
    clf
    c=0 ;
    for DsType = 1:length(DsTypeName) ;
        for cells=1:length(DsCelli{DsType}) ; % for each cell
            c=c+1 ;
            spk = dataRun.spikes{DsCelli{DsType}(cells)}-triggs(t) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;
            plot(spk,ones(1,length(spk))*c,[Color_list{DsType},'.']) 
            hold on
        end
    end
    axis([0,4,0,c])
    title(num2str(t))
    
    pause
end

figure % plot raster of all cells for each trial
for t=1:length(triggs) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        spk = dataRun.spikes{cells}-triggs(t) ;
        spk = spk(spk>=0 & spk<=StimDuration) ;
        plot(spk,ones(1,length(spk))*cells,color)
        axis([0,4,0,length(triggs)]) 
        hold on
    end
    title(num2str(cells))
    hold off
    pause
end

figure % grading ds index vs image slide dsi
for cells=1:length(dataRun.spikes) ; % for each cell
    if ~isnan(cell_list_mat(cells)) ; % if there is a mapped cell
        ci = get_cell_indices(dataRunDs,cell_list_mat(cells)) ;
        plot(mean([ForIgor.dsi_all{1}(ci),ForIgor.dsi_all{2}(ci)]), DSindex(cells),'*')
        hold on
    else
        plot(-.1, DSindex(cells),'*')
        
    end
end
xlabel('grading dsi')
ylabel('image dsi')

figure % noise correlations
for st = 1:length(stimulus.directions) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        for cells2=1:length(dataRun.spikes) ; % for each cell
            plot(sqrt(psth_mean{cells}(st,:).*psth_mean{cells2}(st,:)),squeeze(corr{st}(cells,cells2,:)),'*')
            pause
        end
    end
end
   
figure % noise correlations
for st = 1:length(stimulus.directions) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        for cells2=1:length(dataRun.spikes) ; % for each cell
            plot(psthTime,psth_mean{cells}(st,:),'b')
            hold on
            plot(psthTime,psth_mean{cells2}(st,:),'r')
            plot(psthTime,squeeze(corr{st}(cells,cells2,:)),'k')
            pause
            clf
        end
    end
end

figure % DS cell EI positions  
for DsType = 1:length(DsTypeName) ;
    subplot(length(DsTypeName),1,DsType)
    %figure
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        plot(EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),2),[Color_list{DsType},'+'])
        drawCircle(EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),2),50,'color',Color_list{DsType})
        hold on
    end
end

figure % DS EI positions firing rate movie
TypeAngle = [270,0,90,180] ;
EiSpikeRateFactor = 50 ; % (um/spike)
for st = 1:length(stimulus.directions) ; 
    title(num2str(st))
    for PsthPnt=1:length(psthTime) ; % for each point in psth
        for DsType = 1:length(DsTypeName) ;
            for cells=1:length(DsCelli{DsType}) ; % for each cell
                plot(EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),2),[Color_list{DsType},'o'])
                hold on
                
                spikeNorm = EiSpikeRateFactor*psth_mean{DsCelli{DsType}(cells)}(st,PsthPnt) ; % normalized spike rate

                deltaX = spikeNorm*cosd(TypeAngle(DsType)) ;
                deltaY = spikeNorm*sind(TypeAngle(DsType)) ;

                X=[EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),1)+deltaX] ;
                Y=[EiCtr(DsCelli{DsType}(cells),2),EiCtr(DsCelli{DsType}(cells),2)+deltaY] ;

                plot(X,Y,'Color',Color_list{DsType})
            end
        end
        pause
        clf
    end 
end

% for joel

ForJoel(DB,ImPathNum,DsPathNum).stimulus = struct(stimulus) ;
ForJoel(DB,ImPathNum,DsPathNum).psth = psth ; % {cell}{stimulus}(trial,:)
ForJoel(DB,ImPathNum,DsPathNum).psth_mean = psth_mean ; % {cell}(stimulus,:)
ForJoel(DB,ImPathNum,DsPathNum).psth_timeBins = psthTime ; % 
ForJoel(DB,ImPathNum,DsPathNum).EiCenters = EiCtr ; % (cellX,cellY)
ForJoel(DB,ImPathNum,DsPathNum).DsCellsIndices = DsCelli ; % {DsType}(cell)
ForJoel(DB,ImPathNum,DsPathNum).GratingDsi = ForIgor.dsi_all ; % {grating speed}(cell)
ForJoel(DB,ImPathNum,DsPathNum).psthDg_timeBins = psthTimeDg ; % 
ForJoel(DB,ImPathNum,DsPathNum).psthDg = psthDg ; % {cell}{stimulus}(trial,:)
ForJoel(DB,ImPathNum,DsPathNum).psthDg_mean = psthDg_mean ; % {cell}(stimulus,:)
ForJoel(DB,ImPathNum,DsPathNum).stimulusDg = dataRunDg.stimulus.combinations ; %


