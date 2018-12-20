function ForIgor = MbDrugAnalysis(DataBlock, DB, Params)

% this function will examine responses to moving bars across time and drug
% conditions

% JC 2018-01-10

RunAsScript = true ;
if RunAsScript ;
    DB = 51 ;
    [DataBlock,Params] = DataBlocks_KO ;
end

MbPathNum = 1 ; 

psthBinTime = .001 ; % (sec) size of psthNi Bin step
psthSmoothTime = .02 ; % (sec) size of psthNi Bin (smooth sliding window)

% load data
dataRun = load_data(DataBlock(DB).MbPath{MbPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION
dataRun = load_ei(dataRun, 'all') ;
NumCells = length(dataRun.spikes) ; % number of cells

% load stim 
slashi = strfind(DataBlock(DB).MbPath{MbPathNum},'/') ; % find the /
dataRun.names.stimulus_path = [DataBlock(DB).MbPath{MbPathNum}(1:slashi(end-1)),...
    'stimuli/s',DataBlock(DB).MbPath{MbPathNum}(end-1:end),'.txt'] ;
dataRun = load_stim(dataRun, 'user_defined_trigger_set', [1:2:length(dataRun.triggers)]) ;

TrialTrigInterval = min(diff(dataRun.stimulus.triggers)) ;

NumDirMb = length(dataRun.stimulus.params.DIRECTION) ;
Contrasts = unique(dataRun.stimulus.params.RGB) ; % assumes only 1 bg value
NumCntrst = length(Contrasts) ; 

% stim trials
for t=1:length(dataRun.stimulus.trials) ; % for each DG trial
    MbParams(t,1) = dataRun.stimulus.trials(t).DIRECTION ; % direction
    MbParams(t,2) = dataRun.stimulus.trials(t).RGB(1); % temporal period
end

% 
psthSmoothPnts = floor(psthSmoothTime/psthBinTime) ;

psthTimeMb = [0:psthBinTime:TrialTrigInterval] ;
for cells=1:NumCells ; % for each cell
    for cntrst = 1:NumCntrst ; % for each temp period
        for dir = 1:NumDirMb ; % for each direction 
            prm = [dataRun.stimulus.params.DIRECTION(dir),Contrasts(cntrst)] ; % params set
            ti = find((MbParams(:,1)==prm(1)).*(MbParams(:,2)==prm(2))==1) ; % index of triggers for that stim
            
            spikeTrain{cells}{cntrst}{dir} = zeros(length(ti),length(psthTimeMb)) ; % empty spike train
            for t=1:length(ti) ; % for each trial
                trigTimes{cntrst}{dir}(t) = dataRun.stimulus.triggers(ti(t)) ; % time of trigger
                spk = dataRun.spikes{cells}-dataRun.stimulus.triggers(ti(t)) ;
                spk = spk(spk>=0 & spk<=TrialTrigInterval) ;

                spkPnts = 1+floor(spk/psthBinTime) ; % spike points 
                spikeTrain{cells}{cntrst}{dir}(t,spkPnts) = 1 ; % spike train
                
                psthMb{cells}{cntrst}{dir}(t,:) = smooth(spikeTrain{cells}{cntrst}{dir}(t,:), psthSmoothPnts) ; % psthNi
                spkCountMb{cells}{cntrst}{dir}(t) = sum(spikeTrain{cells}{cntrst}{dir}(t,:)) ; % spike count
            end
            psthMb_mean{cells}{cntrst}(dir,:)=mean(psthMb{cells}{cntrst}{dir},1) ;
            psthMb_var{cells}{cntrst}(dir,:)=var(psthMb{cells}{cntrst}{dir},[],1) ;
            
            spkCountMb_mean{cells}{cntrst}(dir)=mean(spkCountMb{cells}{cntrst}{dir}) ;
            spkCountMb_var{cells}{cntrst}(dir)=var(spkCountMb{cells}{cntrst}{dir}) ;
        end
    end
end

% id drug switch trials from trig times
for cntrst = 1:NumCntrst ; % for each temp period
    for dir = 1:NumDirMb ; % for each direction
        for sw = 1:length(DataBlock(DB).MbPathDrugTimes{MbPathNum}) ; % for each drug switch
            trigTrials{cntrst}{dir}(sw) = find(trigTimes{cntrst}{dir}>DataBlock(DB).MbPathDrugTimes{MbPathNum}(sw),1,'first') ;
        end
    end
end

%% tuning curve avearged over X# trials before drug switch
trialNumAv = 4 ; % number of trial to average over

% find trials for each drug condition
cnd = 1 ;
cndBlocks{cnd} = [1:trigTrials{cntrst}{dir}(cnd)-1] ;

for sw = 1:length(DataBlock(DB).MbPathDrugTimes{MbPathNum}) ; % for each drug switch
    cnd = cnd+1 ;
    if sw~=length(DataBlock(DB).MbPathDrugTimes{MbPathNum}) ;
        cndBlocks{cnd} = [trigTrials{cntrst}{dir}(sw):trigTrials{cntrst}{dir}(sw+1)-1] ;
    else
        cndBlocks{cnd} = [trigTrials{cntrst}{dir}(sw):size(psthMb{cells}{cntrst}{dir},1)] ;
    end
end
      
% tuning curve avearged over X# trials before drug switch
for cells=1:NumCells ; % for each cell
    for cnd = 1:length(cndBlocks); % for each drug condition
        for cntrst = 1:NumCntrst ; % for each contrast
            for dir = 1:NumDirMb ; % for each direction
                CndRspMean{cells}{cnd}(cntrst,dir) = mean(spkCountMb{cells}{cntrst}{dir}(cndBlocks{cnd}(end-trialNumAv:end))) ;
                CndRspStd{cells}{cnd}(cntrst,dir) = std(spkCountMb{cells}{cntrst}{dir}(cndBlocks{cnd}(end-trialNumAv:end))) ;
            end
        end
    end
end

%% mse between between condition X and previous condition

for cnd = 1:length(cndBlocks); % for each drug condition
    for cells=1:NumCells ; % for each cell
        for cntrst = 1:NumCntrst ; % for each temp period
            for dir = 1:NumDirMb ; % for each direction
                 psthCndMean(cnd,cells,cntrst,dir,:)  = mean(psthMb{cells}{cntrst}{dir}(cndBlocks{cnd},:)) ;
                 psthCndVar(cnd,cells,cntrst,dir,:)  = var(psthMb{cells}{cntrst}{dir}(cndBlocks{cnd},:)) ;
                if cnd==1 ;
                    psthCndMse(cnd,cells,cntrst,dir) = 0 ;
                else
                    sig = max(mean(psthCndMean(cnd,cells,cntrst,dir,:)),mean(psthCndMean(cnd-1,cells,cntrst,dir,:))) ;
                    if sig>0 ; % if either condition has spikes
                        psthCndMse(cnd,cells,cntrst,dir) = mean((psthCndMean(cnd,cells,cntrst,dir,:) - ...
                            psthCndMean(cnd-1,cells,cntrst,dir,:)).^2)/sig;
                    else
                        psthCndMse(cnd,cells,cntrst,dir) = 0 ;
                    end
                end
            end
        end
    end
    psthCndMseAv(cnd,:) = squeeze(mean(mean(psthCndMse(cnd,:,:,:),3),4)) ; % average across cntrst and direction
end

%% map location of cells on epiflourescence image
numElectrodeLayers = 2 ;

for cells=1:NumCells ;
    ctr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end
Tform = map_Ei_to_camera(DataBlock(DB).EpiImagePath, DataBlock(DB).MbPath{MbPathNum}) ; % get transform to camera

%figure
im_array = imread(DataBlock(DB).EpiImagePath);

figure
imshow(im_array);
hold on
for cells=1:NumCells ;
    if ismember(cells,DataBlock(DB).PutativeDs)
        col = 'r' ;
    else
        col = 'c' ;
    end
    if ismember(cells,checkOn)
        col = 'b' ;
    end
    
    ctrImage = tformfwd(Tform,ctr(cells,1),ctr(cells,2)) ;
    plot(ctrImage(1),ctrImage(2),'o','MarkerSize',10*psthCndMseAv(2,cells)/mean(psthCndMseAv(2,:)),'Color',col)
end

%% figures
FigSavePath = ['/Users/jcafaro/Documents/AnalysisFigures/MbDrugAnalyis/Db',num2str(DB),'/'] ;  
mkdir(FigSavePath)
color_list = {'b','r','g','k','c'} ;

figure % spike rasters
for cells=1:NumCells ; % for each cell
    r=1 ;
    for cntrst = 1:NumCntrst ; % for each temp period
        for dir = 1:NumDirMb ; % for each direction
            subplot(NumCntrst,NumDirMb,r)
            r=r+1 ;
            
            for t=1:length(ti) ; % for each trial
                plot(psthTimeMb,spikeTrain{cells}{cntrst}{dir}(t,:)*t,'k.') 
                hold on
            end

            axis([0 psthTimeMb(end) 0 length(ti)])
            line([zeros(1,length(trigTrials{cntrst}{dir}));zeros(1,length(trigTrials{cntrst}{dir}))+psthTimeMb(end)],...
                [trigTrials{cntrst}{dir};trigTrials{cntrst}{dir}])

            hold off
        end
    end
    subplot(NumCntrst,NumDirMb,1)
    title(['celli: ',num2str(cells)])
    hgsave([FigSavePath,'rasterCelli',num2str(cells)])
    print(gcf, '-dtiff', [FigSavePath,'rasterCelli',num2str(cells)])
    %pause
end

figure % spike number tuning curves and max and min measures
%for cells=1:NumCells ; % for each cell
for cells= DataBlock(DB).PutativeDs ;
    % mean tuning curves
    subplot(5,1,1)
    for cntrst = 1:NumCntrst ; % for each contrast
        plot(spkCountMb_mean{cells}{cntrst})
        hold on
        [temp,mini(cntrst)] = min(spkCountMb_mean{cells}{cntrst}) ;
        [temp,maxi(cntrst)] = max(spkCountMb_mean{cells}{cntrst}) ;
    end
    hold off
    
    % tuning curves in time
    for cntrst = 1:NumCntrst ; % for each contrast
        subplot(5,1,cntrst+1)
        for t=1:length(ti) ; % for each trial
            plot(t,spkCountMb{cells}{cntrst}{mini(cntrst)}(t),'k*')
            hold on
            plot(t,spkCountMb{cells}{cntrst}{maxi(cntrst)}(t),'ro')
        end
        line([trigTrials{cntrst}{mini(cntrst)};trigTrials{cntrst}{mini(cntrst)}],...
            [zeros(1,length(trigTrials{cntrst}{dir}));zeros(1,length(trigTrials{cntrst}{mini(cntrst)}))+20])

        hold off
    end
    pause
        
end

figure % tuning curve avearged over X# trials before drug switch
for cells=DataBlock(DB).PutativeDs ;
    for cntrst = 1:NumCntrst ; % for each contrast
        subplot(NumCntrst,1,cntrst)
        %for cnd = 1:length(cndBlocks); % for each drug condition
        for cnd = [1,2,4]; % for each drug condition    
            errorbar([1:NumDirMb],CndRspMean{cells}{cnd}(cntrst,:),CndRspStd{cells}{cnd}(cntrst,:))
            hold on
        end
        hold off
    end
    title(num2str(cells))
    pause
end

figure ; % compare max and min for two drug conditions
cnd1=2 ; % use as control 
cnd2=4 ; % 
for cells=DataBlock(DB).PutativeDs ; % for each putative DS cell
    for cntrst = 1:NumCntrst ; % for each contrast
        [Cnd1Min,Cnd1Mini] = min(CndRspMean{cells}{cnd1}(cntrst,:)) ; % min for condition 1
        [Cnd1Max,Cnd1Maxi] = max(CndRspMean{cells}{cnd1}(cntrst,:)) ; % max for condition 1
        Cnd2Min = CndRspMean{cells}{cnd2}(cntrst,Cnd1Mini) ; % value at condtion 1 min
        Cnd2Max = CndRspMean{cells}{cnd2}(cntrst,Cnd1Maxi) ;
        
        subplot(2,NumCntrst,cntrst)
        %plot(Cnd1Min,Cnd2Min,'k*')
        plot(Cnd1Min/Cnd1Max,Cnd2Min/Cnd1Max,'K*')
        hold on
        %plot(Cnd1Min/Cnd1Max,Cnd1Min/Cnd1Max,'r.')
        %plot(Cnd2Min/Cnd1Max,Cnd2Min/Cnd1Max,'r.')
        plot([0,1.5],[0,1.5],'r')
        %set(gca,'xscale','log')
        %set(gca,'yscale','log')
        axis([0,1.5,0,1.5])
        xlabel('cond 1 min')
        ylabel('cond 2 min')
        
        subplot(2,NumCntrst,NumCntrst+cntrst)
        %plot(Cnd1Max,Cnd2Max,'k*')
        plot(Cnd1Max/Cnd1Max,Cnd2Max/Cnd1Max,'k*')
        hold on
        %plot(Cnd1Max,Cnd1Max,'r.')
        %plot(Cnd2Max,Cnd2Max,'r.')
        plot([0,1.5],[0,1.5],'r')
        %set(gca,'xscale','log')
        %set(gca,'yscale','log')
        axis([0,1.5,.5,1.5])
        
        xlabel('cond 1 max')
        ylabel('cond 2 max')
        
    end
end


figure % spike number tuning curves and max and min measures
%for cells=1:NumCells ; % for each cell
for cells= DataBlock(DB).PutativeDs ;

    % tuning curves in time
    for cntrst = 1:NumCntrst ; % for each contrast
        subplot(NumCntrst,1,cntrst)
        
        [temp,mini(cntrst)] = min(spkCountMb_mean{cells}{cntrst}) ;
        [temp,maxi(cntrst)] = max(spkCountMb_mean{cells}{cntrst}) ;
 
        plot(spkCountMb{cells}{cntrst}{mini(cntrst)}./spkCountMb{cells}{cntrst}{maxi(cntrst)}(1),'k-')
        hold on
   
        line([trigTrials{cntrst}{mini(cntrst)};trigTrials{cntrst}{mini(cntrst)}],...
            [zeros(1,length(trigTrials{cntrst}{dir}));zeros(1,length(trigTrials{cntrst}{mini(cntrst)}))+1])
    end     
end


figure % figure for Greg- grant
cells = 153 ;
cntrst = 4 ;
for dir = 1:NumDirMb ; % for each direction
    subplot(2,8,dir)
    r=0 ;
    for cnd = [2,4] ; % for each condition
        
        for t=cndBlocks{cnd} ; % for each trial
            r=r+1 ;
            tempi = find(spikeTrain{cells}{cntrst}{dir}(t,:)==1) ;
            for s = 1:sum(spikeTrain{cells}{cntrst}{dir}(t,:)) ; % for each spike
                temp = psthTimeMb(tempi(s)) ;
                plot([temp,temp],[r-.5,r+.5],color_list{cnd-1}) 
                hold on
            end
            
        end
        axis([0 4 0 25]) 
    end
end

for cnd = [2,4] ; % for each condition
    %subplot(2,8,[9:16])
    %errorbar([1:NumDirMb],CndRspMean{cells}{cnd}(cntrst,:),CndRspStd{cells}{cnd}(cntrst,:),color_list{cnd-1})
    hold on
    polar([dataRun.stimulus.params.DIRECTION*pi/180,dataRun.stimulus.params.DIRECTION(1)*pi/180],...
        [CndRspMean{cells}{cnd}(cntrst,:),CndRspMean{cells}{cnd}(cntrst,1)])
    
end



figure ; % for Greg/Grant
cnd1 = 2 ;
cnd2 = 4 ;
for cells=DataBlock(DB).PutativeDs ; % for each putative DS cell
    [Cnd1Max,Cnd1Maxi] = max(CndRspMean{cells}{cnd1}(cntrst,:)) ; % max for condition 1
    [Cnd1Min,Cnd1Mini] = min(CndRspMean{cells}{cnd1}(cntrst,:)) ; % min for condition 1
    Cnd2Min = CndRspMean{cells}{cnd2}(cntrst,Cnd1Mini) ; % value at condtion 1 min
    Cnd2Max = CndRspMean{cells}{cnd2}(cntrst,Cnd1Maxi) ;
    plot([1,2],[Cnd1Min,Cnd2Min]/Cnd1Max,'k')
    hold on
end

figure ; % compare max and min for two drug conditions
cnd1=2 ; % use as control 
cnd2=4 ; % 
cntrst = 3 ;
%for cells=DataBlock(DB).PutativeDs ; % for each putative DS cell
for cells=1:NumCells ; % for each putative DS cell
    
        [Cnd1Min,Cnd1Mini] = min(CndRspMean{cells}{cnd1}(cntrst,:)) ; % min for condition 1
        [Cnd1Max,Cnd1Maxi] = max(CndRspMean{cells}{cnd1}(cntrst,:)) ; % max for condition 1
        Cnd2Min = CndRspMean{cells}{cnd2}(cntrst,Cnd1Mini) ; % value at condtion 1 min
        Cnd2Max = CndRspMean{cells}{cnd2}(cntrst,Cnd1Maxi) ;
        
        if ismember(cells,DataBlock(DB).PutativeDs)
            %if Cnd1Min/Cnd1Max<=0.4
                plot(Cnd1Min/Cnd1Max,Cnd2Min/Cnd1Max,'k*')
                hold on
            %end
        else
            %plot(Cnd1Min/Cnd1Max,Cnd2Min/Cnd1Max,'k*')
        end
            
    plot([0, 1],[0,1],'r')
end
xlabel('control (min/max control)')
ylabel('+Drug (min/max control)')
print(gcf,'-dpdf','FigForGregPop4')


