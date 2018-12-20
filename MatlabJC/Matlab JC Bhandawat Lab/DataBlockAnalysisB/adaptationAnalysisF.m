function ForIgor = adaptationAnalysisF(Input,A,AnalysisFlag)
% adaptation of "adaptationAnalysisE" 
% added flag to analyze low pass filtered data (ORN lfp and PN subthreshold) instead of spikes
% added absolute concentration (pulse+bg) axis
% added export of pulse - rest data (not just pulse - bg data) 
% added fraction change metric
% removed adaptation factors for bg data
% removed kinetic analysis
% added PID import

% Input.id1 should be a cell array with trial numbers grouped according to vector specified in Input.id2.flag 
% JC 12/18/13

% edited 11/9/15 to make Mac path compatible

for DB = 1:length(A) ; % for each cell

    %% parameters
    PopData.sampleRate = 10000 ; % (hz) temp hard coded - should be saved in file 
    driftCheckTime = 0.25 ; %(sec) time at begining and end which current injected is inspected for changes
    bgTransTime = 5 ; %(sec) time after background odor starts to avoid transient response
    PsthBinTime = 0.1 ; % (sec) time width of psth bin
    OdorRspTime = .05 ; % (sec) time before and after peak of mean psth which odor response is assessed
    PopData.RspTime = 1.5 ; % (sec) time after odor pulse onset over which odor respone ensemble is defined
    
    PopData.SinglePeakTime = .45 ; % (sec) time after odor pulse used as peak time for all data when SinglePeakTimeFlag=1 
    PopData.SinglePeakTimeFlag = false ;
    
    PopData.LpfFlag = false ; % if true than ignore spikes and deal with raw data lowpass filtered 
    PopData.LpfsmoothTime = .05 ; % (sec) time over which data is averaged   

    PopData.PidFlag = true ; % if true than import pid values values and make plots with PID and dilutation values
    PidDataPath = ['/Volumes/bio-home/jc401:/Cafaro Documents/Analysis/PidData/'] ;
    
    ExampleDataSaveFlag = false ; % if true than save example data from all data blocks
    
    PopData.AnalysisFlag = AnalysisFlag ;
    
    if strcmp(PopData.AnalysisFlag,'ORNdata')
        spikeDetectionParameters.numClusters = 1 ;
        spikeDetectionParameters.numPCs = 2 ;
        spikeDetectionParameters.NegDiffThreshStd = 1.25 ; 
    elseif strcmp(PopData.AnalysisFlag,'PNdata')
        spikeDetectionParameters.absRefTime = 0.002 ; % sec (.002)
        spikeDetectionParameters.minRise = 0.2 ; % mV (.3)
        spikeDetectionParameters.minFall = 0.1 ; % mV (.15)
        spikeDetectionParameters.filterOrder = 1 ; % (1)
        spikeDetectionParameters.lpfCutOff = 4000 ; % hz (4000)
        spikeDetectionParameters.minRiseTime = 0.001 ; % sec (.001)
        spikeDetectionParameters.minFallTime = 0.001 ; % sec (.001)
        spikeDetectionParameters.maxPlatTime = 0.001 ; % sec (.001)
    end

    spikeDataPath = ['/Volumes/bio-home/jc401/Cafaro Documents/Analysis/DetectedSpikes/'] ;

    id1 = 'OdorRsp' ;
    id2 = 'OdorConcentration' ;
    id3 = 'BgConcentration' ;
  
    %% stim values
    if strcmp(PopData.AnalysisFlag,'ORNdata')
        PopData.ConcentrationRange  = 10.^[-6:-2] ;  
    elseif strcmp(PopData.AnalysisFlag,'PNdata')
        PopData.ConcentrationRange  = 10.^[-7:-3] ;
    end        
    PopData.BgConcentrationRangeNum = [10^-9,10.^[-7:-4]] ;
    PopData.BgConcentrationRange = {'0','Log7','Log6','Log5','Log4'} ;
    
    if PopData.PidFlag ;
        PIDdata = load([PidDataPath,'PidData']) ;
        if strcmp(PopData.AnalysisFlag,'ORNdata')
            PopData.PidConcentrationRange = PIDdata.PulsePID(2:end) ;  
        elseif strcmp(PopData.AnalysisFlag,'PNdata')
            PopData.PidConcentrationRange = PIDdata.PulsePID(1:end-1) ;
        end        
        PopData.PidBgConcentrationRangeNum = [10^-9,PIDdata.BgPID(2:end)] ; % numerical version, control changed to 10^-x and wash not included
        %PopData.PidBgConcentrationRangeNum = [10^-9,PIDdata.PulsePID(2:end)] ; % if using the pulse peak instead of steady-state
    end     
        
    %% load data in matricies
    rootdir = ['/Volumes/bio-home/jc401/Cafaro Data Backup/', Input(A(DB)).cellname(1:6),'Data'];

    Concentrations = str2num(Input(A(DB)).(id2)) ;
    NumConcentrations = length(Concentrations) ;

    BgConcentrations = Input(A(DB)).(id3) ;
    NumBackgrounds = length(BgConcentrations) ;

    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration

            odorRspTrials{a}{b} = str2num(Input(A(DB)).(id1){a}{b}) ;
            NumTrials(a,b) = length(odorRspTrials{a}{b}) ;
            loopNum = 0 ;
            for c = 1:NumTrials(a,b) ; % for each trial
                loopNum = loopNum+1 ;

                temp = load([rootdir,'/',Input(A(DB)).cellname,'/','voltage_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                vData{a}{b}(loopNum,:) = temp.voltage ;

                temp = load([rootdir,'/',Input(A(DB)).cellname,'/','current_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                iData{a}{b}(loopNum,:) = temp.current ;

                temp = load([rootdir,'/',Input(A(DB)).cellname,'/','Ao0_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                ao0Data{a}{b}(loopNum,:) = temp.Ao0 ; % odor

                temp = load([rootdir,'/',Input(A(DB)).cellname,'/','Ao1_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                ao1Data{a}{b}(loopNum,:) = temp.Ao1 ; % background valve   

                temp = load([rootdir,'/',Input(A(DB)).cellname,'/','TrigTime_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                tData{a}{b}(loopNum) = temp.Trigtime ;
            end
        end
    end
    
    %% time vector
    time = [1:length(vData{1}{1})]/PopData.sampleRate ;

    FirstTime = [] ;
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration 
            if NumTrials(a,b)>0 ;
                FirstTime = min([FirstTime,tData{a}{b}]) ; % earliest time stamp in data used
            end
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration 
            if NumTrials(a,b)>0 ;
                tDataN{a}{b} = (tData{a}{b} - FirstTime)*24*60^2 ; % convert to seconds since experiment began
            end
        end
    end

    %% make odor valve pulse binary
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration 
            if NumTrials(a,b)>0 ;
                ao0DataB{a}{b} = ao0Data{a}{b} ; 
                ao0DataB{a}{b}(ao0Data{a}{b}>=5) = 1 ; 
                ao0DataB{a}{b}(ao0Data{a}{b}<5) = 0 ;
            end
        end
    end

    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration 
            if NumTrials(a,b)>0 ;
                ao1DataB{a}{b} = ao1Data{a}{b} ; 
                ao1DataB{a}{b}(ao1Data{a}{b}>=5) = 1 ; 
                ao1DataB{a}{b}(ao1Data{a}{b}<5) = 0 ;
            end
        end
    end

    %% make sure odor pulse was the same time
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                ao0DataBdiff = ao0DataB{a}{b} - repmat(ao0DataB{1}{1}(1,:),NumTrials(a,b),1) ;
                if sum(abs(ao0DataBdiff(:)))~=0 ;
                    disp('odor pulse discrepancy') ;
                end
            end
        end
    end

    %% check that input current is not changing substantially during any of the trials
    driftCheckPnts = driftCheckTime*PopData.sampleRate ;
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            for c = 1:NumTrials(a,b) ; % for each trial
                Idrift(c) = mean(iData{a}{b}(c,1:driftCheckPnts)) - mean(iData{a}{b}(c,end-driftCheckPnts:end)) ;
                if Idrift(c)>1 ;
                    disp(['significant I drift in trial ',num2str(c)]) ;
                end
            end
        end
    end
    
    %% detect spikes in voltage data
    if ~PopData.LpfFlag ; % 
    
        if strcmp(PopData.AnalysisFlag,'ORNdata')
            try PreviousDetect = load([spikeDataPath,'cell',num2str(A(DB))]) ; 
                if ~isequal(PreviousDetect.spikeDetectionParameters,spikeDetectionParameters) ; % if your previous spike parameters were not the same
                    disp('spike detection parameters assumed are not as specified above')
                end
                spikePnt = PreviousDetect.spikePnt ;

            catch

                % detect and seperate spikes
                for a = 1:NumBackgrounds ; % for each background
                    for b = 1:NumConcentrations ; % for concentration
                        if NumTrials(a,b)>0 ;
                            for c = 1:NumTrials(a,b) ;
                                [spikePntGroup] = spikeSorter(vData{a}{b}(c,:), PopData.sampleRate, spikeDetectionParameters.NegDiffThreshStd, spikeDetectionParameters.numClusters, spikeDetectionParameters.numPCs, false) ;
                                spikePnt{a}{b}{c}= spikePntGroup{1} ;
                            end
                        end
                    end
                end

                save([spikeDataPath,'cell',num2str(A(DB))], 'spikePnt','spikeDetectionParameters') ; % save spike data
            end
        elseif strcmp(PopData.AnalysisFlag,'PNdata')
            try PreviousDetect = load([spikeDataPath,'cell',num2str(A(DB))]) ; 
                if ~isequal(PreviousDetect.spikeDetectionParameters,spikeDetectionParameters) ; % if your previous spike parameters were not the same
                    disp('spike detection parameters assumed are not as specified above')
                end
                spikePnt = PreviousDetect.spikePnt ;

            catch

                for a = 1:NumBackgrounds ; % for each background
                    for b = 1:NumConcentrations ; % for concentration
                        if NumTrials(a,b)>0 ;
                            [TempSpikePnt,SpikeData,NonSpikeData] = spikeFinder(vData{a}{b},PopData.sampleRate,spikeDetectionParameters) ;
                            spikePnt{a}{b}= TempSpikePnt ;
                        end
                    end
                end

                save([spikeDataPath,'cell',num2str(A(DB))], 'spikePnt','spikeDetectionParameters') ; % save spike data
            end
        end

        %% spike trains and psth
        PsthBinPnts = PsthBinTime*PopData.sampleRate ;
        for a = 1:NumBackgrounds ; % for each background
            for b = 1:NumConcentrations ; % for concentration
                if NumTrials(a,b)>0 ;
                    for c = 1:NumTrials(a,b) ;
                        SpikeTrain{a}{b}(c,:)= zeros(1,length(vData{a}{b}(c,:))) ;
                        SpikeTrain{a}{b}(c,spikePnt{a}{b}{c}) = 1 ;
                        SpikeTrainSmooth{a}{b}(c,:) = smooth(SpikeTrain{a}{b}(c,:),PsthBinPnts) ;
                        Psth{a}{b}(c,:) = SpikeTrainSmooth{a}{b}(c,:)*PopData.sampleRate ;

                        Psth{a}{b}(c,1:floor(PsthBinPnts/2),:) = mean(Psth{a}{b}(c,floor(PsthBinPnts/2)+1:PsthBinPnts)) ; % avoid small bin artifacts by ignoring early psth values
                        Psth{a}{b}(c,end-floor(PsthBinPnts/2):end,:) = mean(Psth{a}{b}(c,end-PsthBinPnts:end-floor(PsthBinPnts/2)-1)) ; % avoid small bin artifacts by ignoring late psth values 
                    end
                    Psth_mean{a}{b} = mean(Psth{a}{b},1);
                    Psth_std{a}{b} = std(Psth{a}{b},[],1) ;
                    Psth_sem{a}{b} = Psth_std{a}{b}/sqrt(NumTrials(a,b)) ;
                end
            end
        end
        
        %% lpf data
    elseif PopData.LpfFlag ; % if true than ignore spikes and deal with raw data lowpass filtered 
        LpfsmoothPnts = PopData.LpfsmoothTime*PopData.sampleRate ;
        for a = 1:NumBackgrounds ; % for each background
            for b = 1:NumConcentrations ; % for concentration
                if NumTrials(a,b)>0 ;
                    for c = 1:NumTrials(a,b) ;
                        Psth{a}{b}(c,:) = smooth(vData{a}{b}(c,:), LpfsmoothPnts) ;
                        Psth{a}{b}(c,1:floor(LpfsmoothPnts/2),:) = mean(Psth{a}{b}(c,floor(LpfsmoothPnts/2)+1:LpfsmoothPnts)) ; % avoid small bin artifacts by ignoring early psth values
                        Psth{a}{b}(c,end-floor(LpfsmoothPnts/2):end,:) = mean(Psth{a}{b}(c,end-LpfsmoothPnts:end-floor(LpfsmoothPnts/2)-1)) ; % avoid small bin artifacts by ignoring late psth values 
                    
                        if strcmp(PopData.AnalysisFlag,'ORNdata') ; % if ORN lfp 
                            Psth{a}{b}(c,:) = -Psth{a}{b}(c,:) ; % flip data
                        end
                    end
                    Psth_mean{a}{b} = mean(Psth{a}{b},1);
                    Psth_std{a}{b} = std(Psth{a}{b},[],1) ;
                    Psth_sem{a}{b} = Psth_std{a}{b}/sqrt(NumTrials(a,b)) ;
                end
            end
        end
    end
    
    %% index of background odor and odor pulse
    iopb = find(ao0DataB{2}{1}(1,:)~=0,1,'first')-1 ; % odor pulse begining
    iope = find(ao0DataB{2}{1}(1,:)~=0,1,'last')-1 ; % odor pulse ending

    ibpb = find(ao1DataB{2}{1}(1,:)~=0,1,'first')-1 ; % background pulse beginging
    ibpe = find(ao1DataB{2}{1}(1,:)~=0,1,'last')-1 ; % background pulse end
    
    %% odor response time
    OdorRspPnts = OdorRspTime*PopData.sampleRate ;
    
    for a = 1:NumBackgrounds ;
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ; 
                if PopData.SinglePeakTimeFlag ;
                    mi = PopData.SinglePeakTime * PopData.sampleRate ;
                else
                    [m,mi] = max(Psth_mean{a}{b}(iopb:iopb+PopData.sampleRate)) ; % max point of mean psth within a second of odor pulse onset
                end
                SRiorb(a,b) = mi-1+iopb - OdorRspPnts ; % point of odor pulse begining
                SRiore(a,b) = mi-1+iopb +  OdorRspPnts ; % point of odor pulse end
            end
        end
    end
    
    %% background transient response time
    for a = 1:NumBackgrounds ;
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;         
                [m,mi] = max(Psth_mean{a}{b}(ibpb:ibpb+PopData.sampleRate)) ; % max point of mean psth within a second of odor background onset
                SRibgtb(a,b) = mi-1+ibpb - OdorRspPnts ; % point of odor pulse begining
                SRibgte(a,b) = mi-1+ibpb +  OdorRspPnts ; % point of odor pulse end
            end
        end
    end
    
    %% assess rest spike rate
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration   
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRrestMean{a}{b}(c) = mean(Psth{a}{b}(c,1:ibpb)) ; % first point: background pulse start
                end
                SRrestMean_mean(a,b) = mean(SRrestMean{a}{b}) ;
                SRrestMean_sem(a,b) = std(SRrestMean{a}{b})/sqrt(NumTrials(a,b)) ;
            else
                SRrestMean_mean(a,b) = nan ;
                SRrestMean_sem(a,b) = nan ;
            end
        end
    end

    %% assess background spike rate (during background odor after transient response)
    bgTransPnts = bgTransTime*PopData.sampleRate ;
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRbgMean{a}{b}(c) = mean(Psth{a}{b}(c,ibpb+bgTransPnts:iopb)) ; % (start of background + transient time: begining of odor pulse)
                    SRbgStd{a}{b}(c) = std(Psth{a}{b}(c,ibpb+bgTransPnts:iopb)) ; % variance during background
                end
                SRbgMean_mean(a,b) = mean(SRbgMean{a}{b}) ;
                SRbgStd_mean(a,b) = mean(SRbgStd{a}{b}) ;
                
                SRbgMean_std(a,b) = std(SRbgMean{a}{b}) ;
                SRbgStd_std(a,b) = std(SRbgStd{a}{b}) ;
                
                SRbgMean_sem(a,b) = SRbgMean_std(a,b)/sqrt(NumTrials(a,b)) ;
                SRbgStd_sem(a,b) = SRbgStd_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRbgMean_mean(a,b) = nan ;
                SRbgStd_mean(a,b) = nan ;
                
                SRbgMean_std(a,b) = nan ;
                SRbgStd_std(a,b) = nan ;
                
                SRbgMean_sem(a,b) = nan ;
                SRbgStd_sem(a,b) = nan ;
            end
        end
    end
    
    %% assess background transient spike rate (spike rate peak shortly after background begins)
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRbgTransMean{a}{b}(c) = mean(Psth{a}{b}(c,SRibgtb(a,b):SRibgte(a,b))); %mV (odor response depol start: odor response depol end)
                end
                SRbgTransMean_mean(a,b) = mean(SRbgTransMean{a}{b}) ;
                SRbgTransMean_std(a,b) = std(SRbgTransMean{a}{b}) ;
                SRbgTransMean_sem(a,b) = SRbgTransMean_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRbgTransMean_mean(a,b) = nan ;
                SRbgTransMean_std(a,b) = nan ;
                SRbgTransMean_sem(a,b) = nan ;
            end       
        end
    end
    
    %% odor pulse response spike rate
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRpulseMean{a}{b}(c) = mean(Psth{a}{b}(c,SRiorb(a,b):SRiore(a,b))); %mV (odor response depol start: odor response depol end)
                end
                SRpulseMean_mean(a,b) = mean(SRpulseMean{a}{b}) ;
                SRpulseMean_std(a,b) = std(SRpulseMean{a}{b}) ;
                SRpulseMean_sem(a,b) = SRpulseMean_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRpulseMean_mean(a,b) = nan ;
                SRpulseMean_std(a,b) = nan ;
                SRpulseMean_sem(a,b) = nan ;
            end       
        end
    end

    %% pulse minus rest
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMeanMinRest{a}{b} = SRpulseMean{a}{b} - SRrestMean{a}{b} ;
                SRpulseMeanMinRest_mean(a,b) = mean(SRpulseMeanMinRest{a}{b}) ;
                SRpulseMeanMinRest_std(a,b) = std(SRpulseMeanMinRest{a}{b}) ;
                SRpulseMeanMinRest_std_unc(a,b) = SRpulseMeanMinRest_std(a,b)/sqrt(2*(NumTrials(a,b)-1)) ; % uncertanty in the std measure (std/sqrt(2*(N-1)) p.298 of J. Taylor " Intro to Error Analysis"
                SRpulseMeanMinRest_sem(a,b) = SRpulseMeanMinRest_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRpulseMeanMinRest_mean(a,b) = nan ;
                SRpulseMeanMinRest_std(a,b) = nan ;
                SRpulseMeanMinRest_std_unc(a,b) = nan ;
                SRpulseMeanMinRest_sem(a,b) = nan ;
            end
        end
    end
    
    %% backgound minus rest
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRbgMeanMinRest{a}{b} = SRbgMean{a}{b} - SRrestMean{a}{b} ;
                SRbgMeanMinRest_mean(a,b) = mean(SRbgMeanMinRest{a}{b}) ;
                SRbgMeanMinRest_std(a,b) = std(SRbgMeanMinRest{a}{b}) ;
                SRbgMeanMinRest_sem(a,b) = SRbgMeanMinRest_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRbgMeanMinRest_mean(a,b) = nan ;
                SRbgMeanMinRest_std(a,b) = nan ;
                SRbgMeanMinRest_sem(a,b) = nan ;
            end
        end
        SRbgMeanMinRest_ApMean(a) = mean(cell2mat(SRbgMeanMinRest{a})) ; % mean across of all pulse (Ap) trials
        SRbgMeanMinRest_ApSem(a) = std(cell2mat(SRbgMeanMinRest{a}))/sqrt(sum(NumTrials(a,:))) ; % sem across of all pulse (Ap) trials
        
        SRbgStd_ApMean(a) = mean(cell2mat(SRbgStd{a})) ;
        SRbgStd_ApSem(a) = std(cell2mat(SRbgStd{a}))/sqrt(sum(NumTrials(a,:))) ;
    end

    %% pulse minus bg
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMeanMinBg{a}{b} = SRpulseMean{a}{b} - SRbgMean{a}{b} ;
                SRpulseMeanMinBg_mean(a,b) = mean(SRpulseMeanMinBg{a}{b}) ;
                SRpulseMeanMinBg_std(a,b) = std(SRpulseMeanMinBg{a}{b}) ;
                SRpulseMeanMinBg_std_unc(a,b) = SRpulseMeanMinBg_std(a,b)/sqrt(2*(NumTrials(a,b)-1)) ; % uncertanty in the std measure (std/sqrt(2*(N-1)) p.298 of J. Taylor " Intro to Error Analysis"
                SRpulseMeanMinBg_sem(a,b) = SRpulseMeanMinBg_std(a,b)/sqrt(NumTrials(a,b)) ; % uncertanty in the mean measure
            else
                SRpulseMeanMinBg_mean(a,b) = nan ;
                SRpulseMeanMinBg_std(a,b) = nan ;
                SRpulseMeanMinBg_std_unc(a,b) = nan ;
                SRpulseMeanMinBg_sem(a,b) = nan ;
            end
        end
    end
    
    %% bg transient minus rest
     for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRbgTransMeanMinRest{a}{b} =  SRbgTransMean{a}{b} - SRrestMean{a}{b} ;
                SRbgTransMeanMinRest_mean(a,b) = mean(SRbgTransMeanMinRest{a}{b}) ;
                SRbgTransMeanMinRest_std(a,b) = std(SRbgTransMeanMinRest{a}{b}) ;
                SRbgTransMeanMinRest_sem(a,b) = SRbgTransMeanMinRest_std(a,b)/sqrt(NumTrials(a,b)) ; % uncertanty in the mean measure
            else
                SRbgTransMeanMinRest_mean(a,b) = nan ;
                SRbgTransMeanMinRest_std(a,b) = nan ;
                SRbgTransMeanMinRest_sem(a,b) = nan ; % uncertanty in the mean measure
            end
        end
        SRbgTransMeanMinRest_ApMean(a) = mean(cell2mat(SRbgTransMeanMinRest{a})) ; % mean across of all pulse (Ap) trials
        SRbgTransMeanMinRest_ApSem(a) = std(cell2mat(SRbgTransMeanMinRest{a}))/sqrt(sum(NumTrials(a,:))) ; % sem across of all pulse (Ap) trials
    end
    
    %% (bg+bgt(0))/background transient (adpatation factor for background, bg trans control prevents little adaptating backgrounds from having large adaptation factors)
     for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRbgMeanMinRest_DivBgTrans_Mean(a,b) = (SRbgMeanMinRest_mean(a,b)+SRbgTransMeanMinRest_mean(1,b))./SRbgTransMeanMinRest_mean(a,b) ;
                SRbgMeanMinRest_DivBgTrans_Sem(a,b) = sqrt((sqrt(SRbgMeanMinRest_sem(a,b).^2 +SRbgTransMeanMinRest_sem(1,b).^2)/(SRbgMeanMinRest_mean(a,b)+SRbgTransMeanMinRest_mean(1,b))).^2+...
                   (SRbgTransMeanMinRest_sem(a,b)./SRbgTransMeanMinRest_mean(a,b)).^2)*SRbgMeanMinRest_DivBgTrans_Mean(a,b) ; % sem (has both additive and divisive error propagation
            else
                SRbgMeanMinRest_DivBgTrans_Mean(a,b) = nan ;
                SRbgMeanMinRest_DivBgTrans_Sem(a,b) = nan ;
            end
        end
     end
    
    %% adaptation factors (+bg/control)
     for a = 1:NumBackgrounds ; % for each background

        SRbgStd_ApMean_DivCon(a) = SRbgStd_ApMean(a)/SRbgStd_ApMean(1) ;
        SRbgStd_ApSem_DivCon(a) = SRbgStd_ApMean_DivCon(a)*sqrt((SRbgStd_ApSem(a)/SRbgStd_ApMean(a))^2 + (SRbgStd_ApSem(1)/SRbgStd_ApMean(1))^2) ; % fractional uncertanties add

        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMeanMinBg_mean_DivCon(a,b) = SRpulseMeanMinBg_mean(a,b)/SRpulseMeanMinBg_mean(1,b) ;
                SRpulseMeanMinBg_sem_DivCon(a,b) = SRpulseMeanMinBg_mean_DivCon(a,b)*sqrt((SRpulseMeanMinBg_sem(a,b)/SRpulseMeanMinBg_mean(a,b))^2 + (SRpulseMeanMinBg_sem(1,b)/SRpulseMeanMinBg_mean(1,b))^2) ; % fractional uncertanties add
                
                SRpulseMeanMinBg_std_DivCon(a,b) = SRpulseMeanMinBg_std(a,b)/SRpulseMeanMinBg_std(1,b) ;
                SRpulseMeanMinBg_std_unc_DivCon(a,b) = SRpulseMeanMinBg_std_DivCon(a,b)*sqrt((SRpulseMeanMinBg_std_unc(a,b)/SRpulseMeanMinBg_std(a,b))^2 + (SRpulseMeanMinBg_std_unc(1,b)/SRpulseMeanMinBg_std(1,b))^2) ; % fractional uncertanties add                 
            else
                SRpulseMeanMinBg_mean_DivCon(a,b) = nan ;
                SRpulseMeanMinBg_sem_DivCon(a,b) = nan ;
                
                SRpulseMeanMinBg_std_DivCon(a,b) = nan ;
                SRpulseMeanMinBg_std_unc_DivCon(a,b) = nan ;
            end
        end
        
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMeanMinRest_mean_DivCon(a,b) = SRpulseMeanMinRest_mean(a,b)/SRpulseMeanMinRest_mean(1,b) ;
                SRpulseMeanMinRest_sem_DivCon(a,b) = SRpulseMeanMinRest_mean_DivCon(a,b)*sqrt((SRpulseMeanMinRest_sem(a,b)/SRpulseMeanMinRest_mean(a,b))^2 + (SRpulseMeanMinRest_sem(1,b)/SRpulseMeanMinRest_mean(1,b))^2) ; % fractional uncertanties add
                
                SRpulseMeanMinRest_std_DivCon(a,b) = SRpulseMeanMinRest_std(a,b)/SRpulseMeanMinRest_std(1,b) ;
                SRpulseMeanMinRest_std_unc_DivCon(a,b) = SRpulseMeanMinRest_std_DivCon(a,b)*sqrt((SRpulseMeanMinRest_std_unc(a,b)/SRpulseMeanMinRest_std(a,b))^2 + (SRpulseMeanMinRest_std_unc(1,b)/SRpulseMeanMinRest_std(1,b))^2) ; % fractional uncertanties add                 
            else
                SRpulseMeanMinRest_mean_DivCon(a,b) = nan ;
                SRpulseMeanMinRest_sem_DivCon(a,b) = nan ;
                
                SRpulseMeanMinRest_std_DivCon(a,b) = nan ;
                SRpulseMeanMinRest_std_unc_DivCon(a,b) = nan ;
            end
        end
    end
    
    %% rest and bg activity subtracted psth vectors
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    PsthMinRest{a}{b}(c,:) = Psth{a}{b}(c,:) - SRrestMean{a}{b}(c) ;
                    PsthMinBg{a}{b}(c,:) = Psth{a}{b}(c,:) - SRbgMean{a}{b}(c) ;
                end
            PsthMinRest_mean{a}{b} = mean(PsthMinRest{a}{b},1) ;
            PsthMinBg_mean{a}{b} = mean(PsthMinBg{a}{b},1) ;
            
            PsthMinRest_std{a}{b} = std(PsthMinRest{a}{b},[],1) ;
            PsthMinBg_std{a}{b} = std(PsthMinBg{a}{b},[],1) ;
            
            PsthMinRest_sem{a}{b} = PsthMinRest_std{a}{b}/sqrt(NumTrials(a,b)) ;
            PsthMinBg_sem{a}{b} = PsthMinBg_std{a}{b}/sqrt(NumTrials(a,b)) ;
            end
        end
    end

    %% repository status
    temp = getGitInfo ; 
    RepVer = temp.hash ; % repository version 

    [tempS,tempR] = system('git status') ;
    if length(tempR)==62 ;
        RepStat = ['GitHub up to date ',RepVer] ; % if no unsynced files
    else
        RepStat = ['GitHub not up to date ',RepVer] ;
    end

    %% figures
    Conc(1,:) = [.6,.5,.4,0,0,0,0] ;
    Conc(2,:) = [1,1,1,1,.6,.5,.4] ;
    Conc(3,:) = [.6,.5,.4,.3,.2,.1,0] ;
    for a=1:NumConcentrations ; % each concentration is a matrix and each background is a row within that matrix
        colorMat{a} = [Conc(3,a),Conc(3,a),Conc(3,a); Conc(2,a),Conc(1,a),Conc(1,a); Conc(1,a),Conc(1,a),Conc(2,a);... 
            Conc(1,a),Conc(2,a),Conc(1,a);Conc(2,a),Conc(2,a),Conc(1,a);...
            Conc(2,a),Conc(1,a),Conc(2,a);Conc(1,a),Conc(2,a),Conc(2,a)] ;
    end
        
    % spike detection
    figure 
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            for c = 1:NumTrials(a,b) ;
                plot(time,vData{a}{b}(c,:)) 
                hold on
                plot(time(spikePnt{a}{b}{c}),vData{a}{b}(c,spikePnt{a}{b}{c}),'r*')
                title(num2str(odorRspTrials{a}{b}(c)))
                hold off
                pause
            end
        end
    end
%           
%     % spike raster
%     figure
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             for c = 1:NumTrials(a,b) ;
%                 for d=1:length(spikePnt{a}{b}{c}) ;
%                     plot([1,1]*spikePnt{a}{b}{c}(d),[odorRspTrials{a}{b}(c)-1,odorRspTrials{a}{b}(c)],'Color',colorMat{b}(a,:))
%                     hold on
%                 end
%             end
%         end
%     end
% %     
%     % mean psth for different concentrations and backgrounds and time (raw data)
%     figure 
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(NumConcentrations+2,1,b)
%                 plot(time,Psth_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%                 hold on              
%             end
%         end
%     end
%     axis tight
%     
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(NumConcentrations+2,1,NumConcentrations+1)
%                 plot(tDataN{a}{b},SRbgMean{a}{b},'.','Color',colorMat{b}(a,:))
%                 hold on
%                 plot(tDataN{a}{b},SRpulseMean{a}{b},'o','Color',colorMat{b}(a,:))            
%                 axis tight
%                 
%                 subplot(NumConcentrations+2,1,NumConcentrations+2)
%                 plot(tDataN{a}{b},SRbgMean{a}{b},'.','Color',colorMat{b}(a,:))
%                 hold on
%                 plot(tDataN{a}{b},SRrestMean{a}{b},'*','Color',colorMat{b}(a,:))
%                 axis tight
%             end
%         end
%     end
%     xlabel('trig time (sec)')
%     ylabel('spike rate (hz)')
%     
%     axes('Position',[0 0 .02 1],'Visible','off');
%     text(0,.01,RepStat,'FontSize',5)
%     
% %     % mean psth for different concentrations and backgrounds and time (rest subtracted data)
% %     figure 
% %     for a = 1:NumBackgrounds ; % for each background
% %         for b = 1:NumConcentrations ; % for concentration
% %             if NumTrials(a,b)>0 ;
% %                 subplot(NumConcentrations+2,1,b)
% %                 plot(time,PsthMinRest_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
% %                 hold on   
% %                 axis tight
% %             end
% %         end
% %     end
% %     
% %     for a = 1:NumBackgrounds ; % for each background
% %         for b = 1:NumConcentrations ; % for concentration
% %             if NumTrials(a,b)>0 ;
% %                 subplot(NumConcentrations+2,1,NumConcentrations+2)
% %                 plot(tDataN{a}{b},SRbgMeanMinRest{a}{b},'.','Color',colorMat{b}(a,:))
% %                 hold on
% %                 axis tight
% %                 xlabel('trig time (sec)')
% %                 ylabel('bg sr - rest sr (hz)')
% %                 
% %                 subplot(NumConcentrations+2,1,NumConcentrations+1)
% %                 plot(tDataN{a}{b},SRpulseMeanMinRest{a}{b},'o','Color',colorMat{b}(a,:))
% %                 hold on
% %                 axis tight
% %                 xlabel('trig time (sec)')
% %                 ylabel('pulse sr - rest sr (hz)')
% %             end
% %         end
% %     end
% %     
% %      % mean psth for different concentrations and backgrounds and time (bg subtracted data)
% %     figure 
% %     for a = 1:NumBackgrounds ; % for each background
% %         for b = 1:NumConcentrations ; % for concentration
% %             if NumTrials(a,b)>0 ;
% %                 subplot(NumConcentrations+1,1,b)
% %                 plot(time,PsthMinBg_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
% %                 hold on   
% %                 axis tight
% %             end
% %         end
% %     end
% %     
% %     for a = 1:NumBackgrounds ; % for each background
% %         for b = 1:NumConcentrations ; % for concentration
% %             if NumTrials(a,b)>0 ;
% %                 subplot(NumConcentrations+1,1,NumConcentrations+1)
% %                 plot(tDataN{a}{b},SRpulseMeanMinBg{a}{b},'o','Color',colorMat{b}(a,:))
% %                 hold on
% %                 xlabel('trig time (sec)')
% %                 ylabel('pulse sr - bg sr (hz)')
% %                 axis tight
% %             end
% %         end
% %     end
% % 
% %     % comparing rest and background data across concentrations
% %     figure 
% %     for a = 1:NumBackgrounds ; % for each background
% %         for b = 1:NumConcentrations ; % for concentration
% %             if NumTrials(a,b)>0 ;
% %                 subplot(1,2,1)
% %                 plot(time,Psth_mean{a}{b},'Color',colorMat{b}(a,:))            
% %                 hold on     
% %     
% %                 subplot(1,2,2)
% %                 plot(time,PsthMinRest_mean{a}{b},'Color',colorMat{b}(a,:))            
% %                 hold on                  
% %             end
% %         end
% %     end
% %    
% %     % mean and std spike rate as a function of background
% %     figure
% %     for a = 1:NumBackgrounds ;
% %         subplot(2,1,1)
% %         errorbar(a, SRbgMeanMinRest_ApMean(a),SRbgMeanMinRest_ApSem(a),SRbgMeanMinRest_ApSem(a),'Color',colorMat{NumConcentrations}(a,:))
% %         hold on
% %         xlabel('background concentration')
% %         ylabel('mean sr bg-rest')
% % 
% %         subplot(2,1,2)
% %         errorbar(a, SRbgStd_ApMean(a),SRbgStd_ApSem(a),SRbgStd_ApSem(a),'Color',colorMat{NumConcentrations}(a,:))
% %         hold on
% %         xlabel('background concentration')
% %         ylabel('std sr bg')
% %     end
% %     
% %    % mean and std spike rate as a function of pulse and background
% %     figure
% %     for a = 1:NumBackgrounds ;
% %         subplot(2,1,1)
% %         errorbar(log10(Concentrations), SRpulseMeanMinBg_mean(a,:),SRpulseMeanMinBg_sem(a,:),SRpulseMeanMinBg_sem(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
% %         hold on
% %         xlabel('pulse concentration')
% %         ylabel('mean sr pulse-bg')
% % 
% %         subplot(2,1,2)
% %         plot(log10(Concentrations), SRpulseMeanMinBg_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
% %         hold on
% %         xlabel('background concentration')
% %         ylabel('std sr bg')  
% %     end
% % 
% %     %  spikes as a function of concentration
% %     figure
% %     subplot(2,1,1)
% %     for a = 1:NumBackgrounds ;
% %         errorbar([-10,log10(Concentrations)],[SRbgMeanMinRest_mean_DBwMean(a), SRpulseMeanMinRest_mean(a,:)],[SRbgStd_mean_DBwMean(a),SRpulseMean_std(a,:)],[SRbgStd_mean_DBwMean(a),SRpulseMean_std(a,:)],'*-','Color',colorMat{NumConcentrations}(a,:))
% %         hold on
% %     end
% %     ylabel('spike rate')
% %     xlabel('log concentration')
% % 
% %     subplot(2,1,1)
% %     for a = 1:NumBackgrounds ;
% %          errorbar([-10,log10(Concentrations)],[SRbgMeanMinBg_mean_DBwMean(a), SRpulseMeanMinBg_mean(a,:)],[SRbgStd_mean_DBwMean(a),SRpulseMeanMinBg_std(a,:)],[SRbgStd_mean_DBwMean(a),SRpulseMeanMinBg_std(a,:)],'*-','Color',colorMat{NumConcentrations}(a,:))
% %         hold on
% %     end
% %     ylabel('pulse min bg spike rate')
% %     xlabel('log concentration')
% % 
% % graphs for single page focus
%     MaxNumConcentrations = 6 ;
%     figure
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(MaxNumConcentrations+2,2,b*2-1)
%                 plot(time,Psth_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%                 axis tight
%                 ylabel(num2str(Concentrations(b)))
%                 hold on 
%             end
%         end
%     end
%     
%     subplot(MaxNumConcentrations+2,2,2:2:4)
%     plot(time(1,iopb:SRiore(1,NumConcentrations)),vData{1}{NumConcentrations}(1,iopb:SRiore(1,NumConcentrations)))
%     axis tight
%     
%     subplot(MaxNumConcentrations+2,2,6:2:8)
%     for a = 1:NumBackgrounds ;
%         errorbar(log10(Concentrations), SRpulseMeanMinBg_mean(a,:),SRpulseMeanMinBg_sem(a,:),SRpulseMeanMinBg_sem(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%     end
%     ylabel('pulse min rest spike rate')
%     xlabel('log concentration')
%     axis tight
% 
%     subplot(MaxNumConcentrations+2,2,10:2:12)
%     for a = 1:NumBackgrounds ;
%         errorbar(log10(Concentrations), SRpulseMeanMinBg_mean_DivCon(a,:),SRpulseMeanMinBg_sem_DivCon(a,:),SRpulseMeanMinBg_sem_DivCon(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%     end
%     ylabel('pulse +bg/-bg')
%     xlabel('log concentration')
%     axis tight
%     
%     subplot(MaxNumConcentrations+2,2,2*MaxNumConcentrations+1:2*MaxNumConcentrations+2)
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 plot(tDataN{a}{b},SRpulseMean{a}{b},'o','Color',colorMat{b}(a,:))
%                 hold on
%             end
%         end
%     end
%     xlabel('trig time (sec)')
%     ylabel('spike rate (hz)')
%     axis tight
%     
%     subplot(MaxNumConcentrations+2,2,2*MaxNumConcentrations+3:2*MaxNumConcentrations+4)
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 plot(tDataN{a}{b},SRbgMean{a}{b},'.','Color',colorMat{b}(a,:))
%                 hold on
%                 plot(tDataN{a}{b},SRrestMean{a}{b},'*','Color',colorMat{b}(a,:))
%             end
%         end
%     end
%     xlabel('trig time (sec)')
%     ylabel('spike rate (hz)')
%     axis tight    
%     
%      axes('Position',[0 0 .02 1],'Visible','off');
%      text(0,.02,['cell ',num2str(A(DB)),'   ',PopData.AnalysisFlag,' date', Input(A(DB)).cellname],'FontSize',10)

    %% temp results for population analysis
    
    % range of background and pulse concentrations assessed across
    % population and indicies to compartmentalize this DB correctly
    PopData.time = time ;  
    
    for a = 1:length(Concentrations) ; % for each concentration test for this cell
        if ismember(Concentrations(a),PopData.ConcentrationRange) ; % if that concentration was a selected in PopData
            [c,PopPulsei(a)] = intersect(PopData.ConcentrationRange,Concentrations(a)) ;
        else
            PopPulsei(a) = nan ;
        end
    end
    
    for a = 1:NumBackgrounds ;
        for b = 1:length(PopData.BgConcentrationRange) ;
            if strcmp(BgConcentrations{a},PopData.BgConcentrationRange{b}) ;
                PopBgi(a) = b ;
            end
        end
    end
    
    % NumTrials {bg}(cell,pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            for b=1:length(PopData.ConcentrationRange) ;
                PopData.NumTrials{a} = nan(length(A),length(PopData.ConcentrationRange));
            end
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                PopData.NumTrials{PopBgi(a)}(DB,PopPulsei(b)) = NumTrials(a,b) ;
            end
        end
    end 
    
    % Psth {bg}{pulse}(cell,time)  
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            for b=1:length(PopData.ConcentrationRange) ;
                PopData.Psth_mean{a}{b} = single(nan(length(A),length(time))) ;
                PopData.PsthMinRest_mean{a}{b} = single(nan(length(A),length(time))) ;
                PopData.PsthMinBg_mean{a}{b} = single(nan(length(A),length(time))) ;
                
                PopData.Psth_sem{a}{b} = single(nan(length(A),length(time))) ;
                PopData.PsthMinRest_sem{a}{b} = single(nan(length(A),length(time))) ;
                PopData.PsthMinBg_sem{a}{b} = single(nan(length(A),length(time))) ;
            end
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                if NumTrials(a,b)>0 ;
                    PopData.Psth_mean{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(Psth_mean{a}{b}) ;
                    PopData.PsthMinRest_mean{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(PsthMinRest_mean{a}{b}) ;
                    PopData.PsthMinBg_mean{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(PsthMinBg_mean{a}{b}) ;

                    PopData.Psth_sem{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(Psth_sem{a}{b}) ;
                    PopData.PsthMinRest_sem{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(PsthMinRest_sem{a}{b}) ;
                    PopData.PsthMinBg_sem{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(PsthMinBg_sem{a}{b}) ;
                end
            end
        end
    end 
    
    % SR pulse-bg data {bg}(cell, pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRpulseMeanMinBg_mean{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinBg_sem{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinBg_std{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinBg_std_unc{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                PopData.SRpulseMeanMinBg_mean{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinBg_mean(a,b) ;
                PopData.SRpulseMeanMinBg_sem{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinBg_sem(a,b) ;
                PopData.SRpulseMeanMinBg_std{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinBg_std(a,b) ;
                PopData.SRpulseMeanMinBg_std_unc{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinBg_std_unc(a,b) ;
            end
        end
    end
    
    % SR pulse-rest data {bg}(cell, pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRpulseMeanMinRest_mean{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinRest_sem{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinRest_std{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinRest_std_unc{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                PopData.SRpulseMeanMinRest_mean{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinRest_mean(a,b) ;
                PopData.SRpulseMeanMinRest_sem{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinRest_sem(a,b) ;
                PopData.SRpulseMeanMinRest_std{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinRest_std(a,b) ;
                PopData.SRpulseMeanMinRest_std_unc{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinRest_std_unc(a,b) ;
            end
        end
    end
    
    % SR background data {bg} (cell)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRbgMeanMinRest_ApMean{a} = nan(1,length(A)) ;
            PopData.SRbgMeanMinRest_ApSem{a} = nan(1,length(A)) ;
            PopData.SRbgStd_ApMean{a} = nan(1,length(A)) ;
            PopData.SRbgStd_ApSem{a} = nan(1,length(A)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        PopData.SRbgMeanMinRest_ApMean{PopBgi(a)}(DB) = SRbgMeanMinRest_ApMean(a) ;
        PopData.SRbgMeanMinRest_ApSem{PopBgi(a)}(DB) = SRbgMeanMinRest_ApSem(a) ;
        PopData.SRbgStd_ApMean{PopBgi(a)}(DB) = SRbgStd_ApMean(a) ;
        PopData.SRbgStd_ApSem{PopBgi(a)}(DB) = SRbgStd_ApSem(a) ;
    end
    
    % normilation factors
    
    % SR pulse-bg data {bg}(cell, pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRpulseMeanMinBg_mean_DivCon{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinBg_sem_DivCon{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinBg_std_DivCon{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinBg_std_unc_DivCon{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                PopData.SRpulseMeanMinBg_mean_DivCon{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinBg_mean_DivCon(a,b) ;
                PopData.SRpulseMeanMinBg_sem_DivCon{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinBg_sem_DivCon(a,b) ;
                PopData.SRpulseMeanMinBg_std_DivCon{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinBg_std_DivCon(a,b) ;
                PopData.SRpulseMeanMinBg_std_unc_DivCon{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinBg_std_unc_DivCon(a,b) ;
            end
        end
    end
    
    % SR pulse-rest data {bg}(cell, pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRpulseMeanMinRest_mean_DivCon{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinRest_sem_DivCon{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinRest_std_DivCon{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinRest_std_unc_DivCon{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                PopData.SRpulseMeanMinRest_mean_DivCon{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinRest_mean_DivCon(a,b) ;
                PopData.SRpulseMeanMinRest_sem_DivCon{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinRest_sem_DivCon(a,b) ;
                PopData.SRpulseMeanMinRest_std_DivCon{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinRest_std_DivCon(a,b) ;
                PopData.SRpulseMeanMinRest_std_unc_DivCon{PopBgi(a)}(DB,PopPulsei(b)) = SRpulseMeanMinRest_std_unc_DivCon(a,b) ;
            end
        end
    end
    
    % SR std background data {bg} (cell)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRbgStd_ApMean_DivCon{a} = nan(1,length(A)) ;
            PopData.SRbgStd_ApSem_DivCon{a} = nan(1,length(A)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        PopData.SRbgStd_ApMean_DivCon{PopBgi(a)}(DB) = SRbgStd_ApMean_DivCon(a) ;
        PopData.SRbgStd_ApSem_DivCon{PopBgi(a)}(DB) = SRbgStd_ApSem_DivCon(a) ;
    end
    
    % SR background (for each pulse - all data is combined above) {bg}(cell,pulse)
     if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRbgMeanMinRest_mean{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
         for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                PopData.SRbgMeanMinRest_mean{PopBgi(a)}(DB,PopPulsei(b)) =  SRbgMeanMinRest_mean(a,b) ;
            end
         end
    end
    
    % SR bg transient adaptation factors{bg}(cell,pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRbgMeanMinRest_DivBgTrans_Mean{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRbgMeanMinRest_DivBgTrans_Sem{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
         for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                PopData.SRbgMeanMinRest_DivBgTrans_Mean{PopBgi(a)}(DB,PopPulsei(b)) =  SRbgMeanMinRest_DivBgTrans_Mean(a,b) ;
                PopData.SRbgMeanMinRest_DivBgTrans_Sem{PopBgi(a)}(DB,PopPulsei(b)) =  SRbgMeanMinRest_DivBgTrans_Sem(a,b) ;
            end
         end
    end
    
    % SR bg transient {bg}(cell,pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRbgTransMeanMinRest_mean{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
         for b = 1:NumConcentrations ; % for concentration
            if ~isnan(PopPulsei(b)) ;
                PopData.SRbgTransMeanMinRest_mean{PopBgi(a)}(DB,PopPulsei(b)) =  SRbgTransMeanMinRest_mean(a,b) ;
            end
         end
    end
     
    if ExampleDataSaveFlag ;
    % example data {bg}{pulse}(cell,time)  
        if DB==1 ;
            for a=1:length(PopData.BgConcentrationRange) ;
                for b=1:length(PopData.ConcentrationRange) ;
                    PopData.Psth_Ex{a}{b} = single(nan(length(A),length(time))) ;
                    PopData.PsthMinRest_Ex{a}{b} = single(nan(length(A),length(time))) ;
                    PopData.PsthMinBg_Ex{a}{b} = single(nan(length(A),length(time))) ;

                end
            end
        end
        for a = 1:NumBackgrounds ; % for each background
            for b = 1:NumConcentrations ; % for concentration
                if ~isnan(PopPulsei(b)) ;
                    if NumTrials(a,b)>0 ;
                        PopData.Psth_Ex{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(Psth{a}{b}(1,:)) ;
                        PopData.PsthMinRest_Ex{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(PsthMinRest{a}{b}(1,:)) ;
                        PopData.PsthMinBg_Ex{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(PsthMinBg{a}{b}(1,:)) ;
                    end
                end
            end
        end 

        % example data {bg}{pulse}(cell,time)  
        if DB==1 ;
            for a=1:length(PopData.BgConcentrationRange) ;
                for b=1:length(PopData.ConcentrationRange) ;
                    PopData.Vdata_Ex{a}{b} = single(nan(length(A),length(time))) ;                
                end
            end
        end
        for a = 1:NumBackgrounds ; % for each background
            for b = 1:NumConcentrations ; % for concentration
                if ~isnan(PopPulsei(b)) ;
                    if NumTrials(a,b)>0 ;
                        PopData.Vdata_Ex{PopBgi(a)}{PopPulsei(b)}(DB,:) = single(vData{a}{b}(1,:)) ;
                    end
                end
            end
        end
    end
    
    disp(DB)
    
    %% clear uneeded variables 
    clearvars -except Input A ForIgor PopData iopb AnalysisFlag  
    
 end % Data block loop

%% population analysis

%% psth
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.ConcentrationRange) ; % for each pulse strength
            PopData.PsthMinRest_mean_AdMean{a}(b,:) = nanmean(PopData.PsthMinRest_mean{a}{b});
            PopData.PsthMinBg_mean_AdMean{a}(b,:) = nanmean(PopData.PsthMinBg_mean{a}{b});
            PopData.Psth_mean_AdMean{a}(b,:) = nanmean(PopData.Psth_mean{a}{b}) ;
        end
    end
end

%% background transient normalized PSTH mean
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.ConcentrationRange) ; % for each pulse strength
            PopData.PsthMinRest_mean_AdMean_TransNorm{a}(b,:) = nanmean(PopData.PsthMinRest_mean{a}{b}./repmat(PopData.SRbgTransMeanMinRest_mean{a}(:,b),1,length(PopData.time)));
        end
    end
end

%% bg spike rate mean
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        PopData.SRbgMeanMinRest_ApMean_AdbWMean(a) = nansum(PopData.SRbgMeanMinRest_ApMean{a}./PopData.SRbgMeanMinRest_ApSem{a}.^2)./nansum(1./PopData.SRbgMeanMinRest_ApSem{a}.^2) ; % wieghted mean across data blocks
        %PopData.SRbgMeanMinRest_ApMean_AdbSem(a) = 1./sqrt(nansum(1./PopData.SRbgMeanMinRest_ApSem{a}.^2)) ; % error of wieghted mean
        PopData.SRbgMeanMinRest_ApMean_AdbSem(a) = nanstd(PopData.SRbgMeanMinRest_ApMean{a})./sqrt(sum(~isnan(PopData.SRbgMeanMinRest_ApMean{a}))) ; % sem for non-wieghted mean 
    end
end

%% bg spike rate variability 
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        PopData.SRbgStd_ApMean_AdbWMean(a) = nansum(PopData.SRbgStd_ApMean{a}./PopData.SRbgStd_ApSem{a}.^2)./nansum(1./PopData.SRbgStd_ApSem{a}.^2) ; % wieghted mean across data blocks
        %PopData.SRbgStd_ApMean_AdbSem(a) = 1./sqrt(nansum(1./PopData.SRbgStd_ApSem{a}.^2)) ; % error of wieghted mean
        PopData.SRbgStd_ApMean_AdbSem(a) = nanstd(PopData.SRbgStd_ApMean{a})./sqrt(sum(~isnan(PopData.SRbgStd_ApMean{a}))) ; % sem for non-wieghted mean 
    end
end

%% bg spike rate signal to noise (mean of background response/std of spike response in absensence of bg)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        PopData.SRbgSNR_mean{a} = PopData.SRbgMeanMinRest_ApMean{a}./PopData.SRbgStd_ApMean{a} ; % signal to noise for each cell (its a vector with nans)
        PopData.SRbgSNR_sem{a} = PopData.SRbgSNR_mean{a}.*sqrt((PopData.SRbgMeanMinRest_ApSem{a}./PopData.SRbgMeanMinRest_ApMean{a}).^2+(PopData.SRbgStd_ApSem{a}./PopData.SRbgStd_ApMean{a}).^2) ; % uncertainty for each cell
        
        PopData.SRbgSNR_mean_AdbWMean(a) =  nansum(PopData.SRbgSNR_mean{a}./PopData.SRbgSNR_sem{a}.^2)./nansum(1./PopData.SRbgSNR_sem{a}.^2) ; % wieghted mean across data blocks
        PopData.SRbgSNR_mean_AdbSem(a) = nanstd(PopData.SRbgSNR_mean{a})./sqrt(sum(~isnan(PopData.SRbgSNR_mean{a}))) ; % sem for non-wieghted mean
    end
end

%% pulse-bg spike rate mean
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinBg_mean_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinBg_mean{a}./PopData.SRpulseMeanMinBg_sem{a}.^2,1)./nansum(1./PopData.SRpulseMeanMinBg_sem{a}.^2,1) ; % wiehgted mean across db
        %PopData.SRpulseMeanMinBg_mean_AdbSem(a,:) = 1./sqrt(nansum(1./PopData.SRpulseMeanMinBg_sem{a}.^2,1)) ; % error of wieghted mean
        PopData.SRpulseMeanMinBg_mean_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinBg_mean{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinBg_mean{a}),1)) ; % sem for non-wieghted mean 
    end
end

%% pulse-rest spike rate mean
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinRest_mean_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinRest_mean{a}./PopData.SRpulseMeanMinRest_sem{a}.^2,1)./nansum(1./PopData.SRpulseMeanMinRest_sem{a}.^2,1) ; % wiehgted mean across db
        %PopData.SRpulseMeanMinRest_mean_AdbSem(a,:) = 1./sqrt(nansum(1./PopData.SRpulseMeanMinRest_sem{a}.^2,1)) ; % error of wieghted mean
        PopData.SRpulseMeanMinRest_mean_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinRest_mean{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinRest_mean{a}),1)) ; % sem for non-wieghted mean 
    end
end

%% pulse-bg spike rate variability
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinBg_std_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinBg_std{a}./PopData.SRpulseMeanMinBg_std_unc{a}.^2,1)./nansum(1./PopData.SRpulseMeanMinBg_std_unc{a}.^2,1) ; % wiehgted mean across db
        %PopData.SRpulseMeanMinBg_std_AdbSem(a,:) = 1./sqrt(nansum(1./PopData.SRpulseMeanMinBg_std_unc{a}.^2,1)) ; % error of wieghted mean
        PopData.SRpulseMeanMinBg_std_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinBg_std{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinBg_std{a}),1)) ; % sem for non-wieghted mean 
    end
end

%% pulse-rest spike rate variability
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinRest_std_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinRest_std{a}./PopData.SRpulseMeanMinRest_std_unc{a}.^2,1)./nansum(1./PopData.SRpulseMeanMinRest_std_unc{a}.^2,1) ; % wiehgted mean across db
        %PopData.SRpulseMeanMinRest_std_AdbSem(a,:) = 1./sqrt(nansum(1./PopData.SRpulseMeanMinRest_std_unc{a}.^2,1)) ; % error of wieghted mean
        PopData.SRpulseMeanMinRest_std_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinRest_std{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinRest_std{a}),1)) ; % sem for non-wieghted mean 
    end
end

%% normalization factors

%% bg spike rate variability 
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        PopData.SRbgStd_ApMean_DivCon_AdbWMean(a) = nansum(PopData.SRbgStd_ApMean_DivCon{a}./PopData.SRbgStd_ApSem_DivCon{a}.^2)./nansum(1./PopData.SRbgStd_ApSem_DivCon{a}.^2) ; % wieghted mean across data blocks
        %PopData.SRbgStd_ApMean_DivCon_AdbSem(a) = 1./sqrt(nansum(1./PopData.SRbgStd_ApSem_DivCon{a}.^2)) ; % error of wieghted mean
        PopData.SRbgStd_ApMean_DivCon_AdbSem(a) = nanstd(PopData.SRbgStd_ApMean_DivCon{a})./sqrt(sum(~isnan(PopData.SRbgStd_ApMean_DivCon{a}))) ; % sem for non-wieghted mean 
    end
end

%% pulse-bg spike rate mean 
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinBg_mean_DivCon{a}./PopData.SRpulseMeanMinBg_sem_DivCon{a}.^2,1)./nansum(1./PopData.SRpulseMeanMinBg_sem_DivCon{a}.^2,1) ; % wiehgted mean across db
        %PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem(a,:) = 1./sqrt(nansum(1./PopData.SRpulseMeanMinBg_sem_DivCon{a}.^2,1)) ; % error of wieghted mean
        PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinBg_mean_DivCon{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinBg_mean_DivCon{a}),1)) ; % sem for non-wieghted mean 
    end
end

%% pulse-rest spike rate mean 
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinRest_mean_DivCon{a}./PopData.SRpulseMeanMinRest_sem_DivCon{a}.^2,1)./nansum(1./PopData.SRpulseMeanMinRest_sem_DivCon{a}.^2,1) ; % wiehgted mean across db
        %PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem(a,:) = 1./sqrt(nansum(1./PopData.SRpulseMeanMinRest_sem_DivCon{a}.^2,1)) ; % error of wieghted mean
        PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinRest_mean_DivCon{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinRest_mean_DivCon{a}),1)) ; % sem for non-wieghted mean 
    end
end

%% pulse-bg spike rate variability 
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinBg_std_DivCon_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinBg_std_DivCon{a}./PopData.SRpulseMeanMinBg_std_unc_DivCon{a}.^2,1)./nansum(1./PopData.SRpulseMeanMinBg_std_unc_DivCon{a}.^2,1) ; % wiehgted mean across db
        %PopData.SRpulseMeanMinBg_std_DivCon_AdbSem(a,:) = 1./sqrt(nansum(1./PopData.SRpulseMeanMinBg_std_unc_DivCon{a}.^2,1)) ; % error of wieghted mean
        PopData.SRpulseMeanMinBg_std_DivCon_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinBg_std_DivCon{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinBg_std_DivCon{a}),1)) ; % sem for non-wieghted mean 
    end
end

%% pulse-rest spike rate variability 
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinRest_std_DivCon_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinRest_std_DivCon{a}./PopData.SRpulseMeanMinRest_std_unc_DivCon{a}.^2,1)./nansum(1./PopData.SRpulseMeanMinRest_std_unc_DivCon{a}.^2,1) ; % wiehgted mean across db
        %PopData.SRpulseMeanMinRest_std_DivCon_AdbSem(a,:) = 1./sqrt(nansum(1./PopData.SRpulseMeanMinRest_std_unc_DivCon{a}.^2,1)) ; % error of wieghted mean
        PopData.SRpulseMeanMinRest_std_DivCon_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinRest_std_DivCon{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinRest_std_DivCon{a}),1)) ; % sem for non-wieghted mean 
    end
end

%% applying normalization factors to mean control data
PopData.SRbgStd_ApMean_DivCon_AdbWMean_TimCon = PopData.SRbgStd_ApMean_AdbWMean(1) * PopData.SRbgStd_ApMean_DivCon_AdbWMean ; 
PopData.SRbgStd_ApMean_DivCon_AdbSem_TimCon = PopData.SRbgStd_ApMean_DivCon_AdbWMean_TimCon .* sqrt((PopData.SRbgStd_ApMean_AdbSem(1)/PopData.SRbgStd_ApMean_AdbWMean(1))^2 + (PopData.SRbgStd_ApMean_DivCon_AdbSem./PopData.SRbgStd_ApMean_DivCon_AdbWMean).^2) ;

for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:) = PopData.SRpulseMeanMinBg_mean_AdbWMean(1,:) .* PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean(a,:) ;
        PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:) = PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:) .* sqrt((PopData.SRpulseMeanMinBg_mean_AdbSem(1,:)./PopData.SRpulseMeanMinBg_mean_AdbWMean(1,:)).^2 + (PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem(a,:)./PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean(a,:)).^2) ;
    
        PopData.SRpulseMeanMinBg_std_DivCon_AdbWMean_TimCon(a,:) = PopData.SRpulseMeanMinBg_std_AdbWMean(1,:) .* PopData.SRpulseMeanMinBg_std_DivCon_AdbWMean(a,:) ;
        PopData.SRpulseMeanMinBg_std_DivCon_AdbSem_TimCon(a,:) = PopData.SRpulseMeanMinBg_std_DivCon_AdbWMean_TimCon(a,:) .* sqrt((PopData.SRpulseMeanMinBg_std_AdbSem(1,:)./PopData.SRpulseMeanMinBg_std_AdbWMean(1,:)).^2 + (PopData.SRpulseMeanMinBg_std_DivCon_AdbSem(a,:)./PopData.SRpulseMeanMinBg_std_DivCon_AdbWMean(a,:)).^2) ;
    end
end

for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:) = PopData.SRpulseMeanMinRest_mean_AdbWMean(1,:) .* PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean(a,:) ;
        PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:) = PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:) .* sqrt((PopData.SRpulseMeanMinRest_mean_AdbSem(1,:)./PopData.SRpulseMeanMinRest_mean_AdbWMean(1,:)).^2 + (PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem(a,:)./PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean(a,:)).^2) ;
    
        PopData.SRpulseMeanMinRest_std_DivCon_AdbWMean_TimCon(a,:) = PopData.SRpulseMeanMinRest_std_AdbWMean(1,:) .* PopData.SRpulseMeanMinRest_std_DivCon_AdbWMean(a,:) ;
        PopData.SRpulseMeanMinRest_std_DivCon_AdbSem_TimCon(a,:) = PopData.SRpulseMeanMinRest_std_DivCon_AdbWMean_TimCon(a,:) .* sqrt((PopData.SRpulseMeanMinRest_std_AdbSem(1,:)./PopData.SRpulseMeanMinRest_std_AdbWMean(1,:)).^2 + (PopData.SRpulseMeanMinRest_std_DivCon_AdbSem(a,:)./PopData.SRpulseMeanMinRest_std_DivCon_AdbWMean(a,:)).^2) ;
    end
end
        
%% fit dose response with sigmoid 
% for a = 1:length(PopData.BgConcentrationRange) ;
%     InitCoefs = [.5,10^-5] ; % initial estimated mean and variance
%     PopData.ConcentrationRangeFit = 10.^[min(log10(PopData.ConcentrationRange)):.1:max(log10(PopData.ConcentrationRange))] ; % interpolated data range
%     if sum(~isnan(PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:)))>1 ; % if there is data
%         Firsti = find(~isnan(PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:))==1,1,'first') ;
%         Lasti = find(~isnan(PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:))==1,1,'last') ;
%         FitParams.MaxY = (max(PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:))+PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,Lasti))/2 ; % amp of sig fit (max+last point)/2         
%         %FitParams.MinY = (min(PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:))+PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,Firsti))/2 ; % amp of sig fit (min+first point)/2
%         FitParams.MinY = 0 ; 
%         FitParams.x = PopData.ConcentrationRange ; % data x values
%         FitParams.Xinterp = PopData.ConcentrationRangeFit ; % interpolated data range
%         PopData.FitCoefs{a} = nlinfit(FitParams,PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:),@SatFun2free2fixed,InitCoefs) ; % find fit params
%         FitParams.x = FitParams.Xinterp ; % temporarally change x to get interpolated data
%         PopData.SigFit(a,:) = SatFun2free2fixed(PopData.FitCoefs{a},FitParams) ; % generate curve
%         MidPoint(a) = log10(interp1(PopData.SigFit(a,:),PopData.ConcentrationRangeFit,(FitParams.MaxY+FitParams.MinY)/2)) ; % interpolate to find midpoint
%     else
%         PopData.FitCoefs{a} = nan(1,2) ;
%         PopData.SigFit(a,:) = nan(1,length(PopData.ConcentrationRangeFit)) ;
%     end
% end

%% fit mean of sigmoids vs. background with a thresholding nonlinearity
% %FitCoefsMat = cell2mat(PopData.FitCoefs') ;
% InitCoefs = [-7,MidPoint(1),.7] ;
% PopData.BgConcentrationRangeNumLogFit = [-9:.1:-4] ;
% PopData.HockeyStick_FitCoefs = nlinfit(log10(PopData.BgConcentrationRangeNum), MidPoint,@HockeyStick,InitCoefs) ;
% PopData.HockeyStickfit = HockeyStick(PopData.HockeyStick_FitCoefs,PopData.BgConcentrationRangeNumLogFit) ;

%% pulse contrast axis
for a = 1:length(PopData.BgConcentrationRangeNum) ; % for each background
    PopData.ConcentrationRangeDivBg(a,:) = PopData.ConcentrationRange/PopData.BgConcentrationRangeNum(a) ;
    %PopData.ConcentrationRangeFitDivBg(a,:) = PopData.ConcentrationRangeFit/PopData.BgConcentrationRangeNum(a) ;
end

if PopData.PidFlag ;
    for a = 1:length(PopData.BgConcentrationRangeNum) ; % for each background
        PopData.PidConcentrationRangeDivBg(a,:) = PopData.PidConcentrationRange/PopData.PidBgConcentrationRangeNum(a) ;
    end
end
    
%% absolute concentration axis
for a = 1:length(PopData.BgConcentrationRangeNum) ; % for each background
    PopData.ConcentrationRangePlusBg(a,:) = PopData.ConcentrationRange+PopData.BgConcentrationRangeNum(a) ;
    %PopData.ConcentrationRangeFitPlusBg(a,:) = PopData.ConcentrationRangeFit+PopData.BgConcentrationRangeNum(a) ;
end

if PopData.PidFlag ;
    for a = 1:length(PopData.BgConcentrationRangeNum) ; % for each background
        PopData.PidConcentrationRangePlusBg(a,:) = PopData.PidConcentrationRange+PopData.PidBgConcentrationRangeNum(a) ;
    end
end
    
%% relative bg and bg transient adaptation
SRpulseMeanMindBg_mean_DivCon_Mat = cell2mat(PopData.SRpulseMeanMinBg_mean_DivCon') ; % matrix from cell array
SRbgMeanMinRest_DivBgTrans_Mean_Mat = cell2mat(PopData.SRbgMeanMinRest_DivBgTrans_Mean') ; % matrix from cell array
for a=1:length(PopData.ConcentrationRange) ; % for each pulse concentration
     SRpulseMeanMindBg_mean_DivCon_MatNorm(:,a)= SRpulseMeanMindBg_mean_DivCon_Mat(:,a)/max(SRpulseMeanMindBg_mean_DivCon_Mat(:,a)) ;
     SRbgMeanMinRest_DivBgTrans_Mean_MatNorm(:,a)= SRbgMeanMinRest_DivBgTrans_Mean_Mat(:,a)/max(SRbgMeanMinRest_DivBgTrans_Mean_Mat(:,a)) ;
end

%% bg transient minus bg
for a = 1:length(PopData.BgConcentrationRangeNum) ; % for each background
    PopData.SRbgTransMeanMinRest_mean_MinBg{a} = PopData.SRbgTransMeanMinRest_mean{a}-PopData.SRbgMeanMinRest_mean{a} ;
end

%% fraction change ((control-adpated)/range)
RangeFactor = max(PopData.SRpulseMeanMinRest_mean_AdbWMean(1,:)) ; % the max range of the cells response

% for each cell find fractional change
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
            PopData.SRpulseMeanMinBg_FracChange{a}(:,b) = (PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(1,b)-...
                (PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,b)*PopData.SRpulseMeanMinRest_mean_DivCon{a}(:,b)))/RangeFactor ;
        end
    end
end

% for the average
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinRest_AdbWMean_FracChange(a,:) = (PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(1,:)...
            -PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:))/RangeFactor ;
        
        PopData.SRpulseMeanMinRest_AdbSem_FracChange(a,:) = sqrt(PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(1,:).^2 ...
            +PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:).^2)/RangeFactor ;
    end
end

%% figures
Conc(1,:) = [.7,.6,.5,.4,0,0,0,0] ;
Conc(2,:) = [1,1,1,1,.7,.6,.5,.4] ;
Conc(3,:) = [.7,.6,.5,.4,.3,.2,.1,0] ;
for a=1:length(PopData.ConcentrationRange) ; % each concentration is a matrix and each background is a row within that matrix
    colorMat{a} = [Conc(3,a),Conc(3,a),Conc(3,a); Conc(2,a),Conc(1,a),Conc(1,a); Conc(1,a),Conc(1,a),Conc(2,a);... 
        Conc(1,a),Conc(2,a),Conc(1,a);Conc(2,a),Conc(2,a),Conc(1,a);...
        Conc(2,a),Conc(1,a),Conc(2,a);Conc(1,a),Conc(2,a),Conc(2,a);...
        Conc(3,a),Conc(3,a),Conc(1,a)] ;
end

figure % psth (-rest) by pulse
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        if ~isnan(sum(PopData.PsthMinRest_mean_AdMean{a}(b,:))) ; % if there is data
            subplot(length(PopData.ConcentrationRange),1,b) ;
%             plot(PopData.time,PopData.PsthMinRest_mean_AdWMean{a}(b,:),'Color',colorMat{5}(a,:))
%             hold on
            %plot(PopData.time(AllCellPeakPnt(a,b)),PopData.PsthMinRest_mean_AdWMean{a}(b,AllCellPeakPnt(a,b)),'o','Color',colorMat{5}(a,:))
            plot(PopData.time,PopData.PsthMinRest_mean_AdMean{a}(b,:),'Color',colorMat{5}(a,:))
            hold on
            
        end
    end
end

figure % psth (-rest) by pulse
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        if ~isnan(sum(PopData.PsthMinRest_mean_AdWMean{a}(b,:))) ; % if there is data
            subplot(length(PopData.ConcentrationRange),1,b) ;           
            plot(PopData.time,PopData.PsthMinBg_mean_AdMean{a}(b,:),'Color',colorMat{5}(a,:))
            hold on    
            ForVikas.ORN.PsthMinRest{a}(b,:) = PopData.PsthMinRest_mean_AdMean{a}(b,:) ;
            ForVikas.ORN.PsthMinBg{a}(b,:) = PopData.PsthMinBg_mean_AdMean{a}(b,:) ;
        end
    end
end

% psth (-rest) by pulse as requested by vikas
for a = 1:length(PopData.ConcentrationRange) ; % for each pulse
    %figure
    for b = 1:length(PopData.BgConcentrationRange) ; % for each bg concentration 
        for c = 1:length(PopData.NumTrials{b}(:,a)) ; 
            if nansum(PopData.NumTrials{b}(:,a))>0 ; % if there is data in at this pulse and background
                %subplot(length(PopData.BgConcentrationRange),1,b) ;
                plot(PopData.time,PopData.PsthMinRest_mean{b}{a}(c,:),'Color',colorMat{5}(b,:))
                hold on
            end
        end
    end
    axes('Position',[0 0 .02 1],'Visible','off');      
    text(0,.02,['log odor pulse= ',num2str(log10(PopData.ConcentrationRange(a)))],'FontSize',10)
end

figure %psth (-bg) pulse kinetics (with only valid pairwise comparisons averaged in)
for a = 2:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        for c = 1:size(PopData.NumTrials{a},1) ;
            if PopData.NumTrials{a}(c,b)>0 ; % if there is data in this cell at this background
                subplot(length(PopData.BgConcentrationRange)-1,length(PopData.ConcentrationRange),sub2ind([length(PopData.ConcentrationRange),length(PopData.BgConcentrationRange)-1],b,a-1))
                plot(PopData.time(150000:165000),PopData.PsthMinBg_mean{a}{b}(c,150000:165000),'Color',colorMat{5}(a,:))
                hold on
                plot(PopData.time(150000:165000),PopData.PsthMinBg_mean{1}{b}(c,150000:165000),'Color',colorMat{5}(1,:))
                axis tight
            end
        end
    end
end

figure % bg respones for each background
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background

            subplot(3,1,1)
            errorbar(a,PopData.SRbgMeanMinRest_ApMean{a}(b),PopData.SRbgMeanMinRest_ApSem{a}(b),PopData.SRbgMeanMinRest_ApSem{a}(b),'Color',colorMat{5}(a,:))
            hold on
            plot(a,PopData.SRbgMeanMinRest_ApMean_AdbWMean(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            ylabel('SRbgMeanMinRest')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end))

            subplot(3,1,2)
            errorbar(a,PopData.SRbgMeanMinRest_ApMean_AdbWMean(a),PopData.SRbgMeanMinRest_ApMean_AdbSem(a),PopData.SRbgMeanMinRest_ApMean_AdbSem(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            plot([1:length(PopData.BgConcentrationRange)],ones(1,length(PopData.BgConcentrationRange))*PopData.SRbgStd_ApMean_AdbWMean(1),'k--')
            ylabel('SRbgMeanMinRest')
            xlabel('background')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end))
            
            subplot(3,1,3)
            errorbar(a, PopData.SRbgSNR_mean_AdbWMean(a), PopData.SRbgSNR_mean_AdbSem(a),PopData.SRbgSNR_mean_AdbSem(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('bg mean / control std')
            xlabel('background')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end))
        end
    end
end

figure % bg variability for each background (non normalized)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background

            subplot(2,1,1)
            errorbar(a,PopData.SRbgStd_ApMean{a}(b),PopData.SRbgStd_ApSem{a}(b),PopData.SRbgStd_ApSem{a}(b),'*','Color',colorMat{5}(a,:))
            hold on
            plot(a,PopData.SRbgStd_ApMean_AdbWMean(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            ylabel('SRbgStd')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)-1])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end-1))

            subplot(2,1,2)
            errorbar(a,PopData.SRbgStd_ApMean_AdbWMean(a),PopData.SRbgStd_ApMean_AdbSem(a),PopData.SRbgStd_ApMean_AdbSem(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('SRbgStd')
            xlabel('background')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)-1])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end-1))
        end
    end
end

figure % dose response from avearges (non normalized)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background
            subplot(2,1,1)
            errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_mean{a}(b,:),PopData.SRpulseMeanMinBg_sem{a}(b,:),PopData.SRpulseMeanMinBg_sem{a}(b,:),'Color',colorMat{5}(a,:))
            hold on
            plot(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_mean_AdbWMean(a,:),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            ylabel('SRpulseMeanMinBg')

            subplot(2,1,2)
            errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_mean_AdbWMean(a,:),PopData.SRpulseMeanMinBg_mean_AdbSem(a,:),PopData.SRpulseMeanMinBg_mean_AdbSem(a,:),'LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('SRpulseMeanMinBg')
            xlabel('log(pulse)')
            
        end
    end
end

figure % pulse variability (non normalized)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_std{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background
            subplot(2,1,1)
            errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_std{a}(b,:),PopData.SRpulseMeanMinBg_std_unc{a}(b,:),PopData.SRpulseMeanMinBg_std_unc{a}(b,:),'Color',colorMat{5}(a,:))
            hold on
            plot(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_std_AdbWMean(a,:),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            ylabel('SRpulseMeanMinBg std')

            subplot(2,1,2)
            errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_std_AdbWMean(a,:),PopData.SRpulseMeanMinBg_std_AdbSem(a,:),PopData.SRpulseMeanMinBg_std_AdbSem(a,:),'LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('SRpulseMeanMinBg std')
            xlabel('log(pulse)')
            
        end
    end
end
 
figure % contrast tuning curves (non normalized)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background
            errorbar(log10(PopData.ConcentrationRangeDivBg(a,:)),PopData.SRpulseMeanMinBg_mean_AdbWMean(a,:),PopData.SRpulseMeanMinBg_mean_AdbSem(a,:),PopData.SRpulseMeanMinBg_mean_AdbSem(a,:),'LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('SRpulseMeanMinBg')
            xlabel('log(pulse/Bg)')    
        end
    end
end

figure % data peak comparisons
for a = 2:length(PopData.BgConcentrationRange); % for each background but not control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background
            for c = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
                subplot(1,length(PopData.ConcentrationRange),c)
                plot(PopData.SRpulseMeanMinBg_mean{1}(b,c),PopData.SRpulseMeanMinBg_mean{a}(b,c),'*','Color',colorMat{5}(a,:))
                hold on
            end    
        end
    end
end
for c = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
    subplot(1,length(PopData.ConcentrationRange),c)
    temp = axis ;
    plot([min(temp),max(temp)],[min(temp),max(temp)],'k')
    xlabel('control spike rate (Hz)')
    ylabel('+bg spike rate (Hz)')
    axis tight
end

% normalization figures

figure % pules-bg adaptation factors and tuning curve
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background

            subplot(3,1,1)
            errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_mean_DivCon{a}(b,:),PopData.SRpulseMeanMinBg_sem_DivCon{a}(b,:),PopData.SRpulseMeanMinBg_sem_DivCon{a}(b,:),'Color',colorMat{5}(a,:))
            plot(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_mean_DivCon{a}(b,:),'*','Color',colorMat{5}(a,:))
            hold on
            xlabel('log(pulse)')
            ylabel('bg/control')
            
            subplot(3,1,2)
            errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem(a,:),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            xlabel('log(pulse)')
            ylabel('bg/control')
            
            subplot(3,1,3)
            errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),'o-','LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            %plot(log10(PopData.ConcentrationRangeFit),PopData.SigFit(a,:),'Color',colorMat{5}(a,:),'LineWidth',1)
            ylabel('SRpulseMeanMinBg')
            xlabel('log(pulse)')
        end
    end
end

figure % background variance
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background

            subplot(2,1,1)
            errorbar(a,PopData.SRbgStd_ApMean_DivCon{a}(b),PopData.SRbgStd_ApSem_DivCon{a}(b),PopData.SRbgStd_ApSem_DivCon{a}(b),'*','Color',colorMat{5}(a,:))
            hold on
            errorbar(a,PopData.SRbgStd_ApMean_DivCon_AdbWMean(a),PopData.SRbgStd_ApMean_DivCon_AdbSem(a),PopData.SRbgStd_ApMean_DivCon_AdbSem(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            xlabel('log(pulse)')
            ylabel('bg/control')
            
            subplot(2,1,2)
            errorbar(a,PopData.SRbgStd_ApMean_DivCon_AdbWMean_TimCon(a),PopData.SRbgStd_ApMean_DivCon_AdbSem_TimCon(a),PopData.SRbgStd_ApMean_DivCon_AdbSem_TimCon(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('SRbg std')
            xlabel('log(pulse)')
        end
    end
end

figure % tuning curves (dilution)
subplot(1,6,1) % (pulse-rest vs absolute concentration)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.ConcentrationRangePlusBg(a,:)),PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),'*-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFitDivBg(a,:)),PopData.SigFit(a,:),'-','Color',colorMat{5}(a,:))
        ylabel('SR pulse-rest')
        xlabel('log(pulse+Bg)')    
    end
end

subplot(1,6,2) % (pulse-bg vs absolute concentration)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.ConcentrationRangePlusBg(a,:)),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),'o-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFit),PopData.SigFit(a,:),'Color',colorMat{5}(a,:),'LineWidth',1)
        ylabel('SR pulse-bg')
        xlabel('log(pulse+bg)')
    end
end

subplot(1,6,3)  % (pulse-rest vs pulse concentration)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),'*-','Color',colorMat{5}(a,:))
        hold on
        ylabel('SR pulse-rest')
        xlabel('log(pulse)')    
    end
end

subplot(1,6,4) % (pulse-bg vs pulse concentration)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),'o-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFit),PopData.SigFit(a,:),'Color',colorMat{5}(a,:),'LineWidth',1)
        ylabel('SR pulse-bg')
        xlabel('log(pulse)')
    end
end

subplot(1,6,5)  % contrast axis (pulse-rest vs contrast)
for a = 2:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.ConcentrationRangeDivBg(a,:)),PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),'*-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFitDivBg(a,:)),PopData.SigFit(a,:),'-','Color',colorMat{5}(a,:))
        ylabel('SR pulse - rest')
        xlabel('log(pulse/Bg)')    
    end
end

subplot(1,6,6)  % contrast axis (pulse-bg vs contrast)
for a = 2:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background   
        errorbar(log10(PopData.ConcentrationRangeDivBg(a,:)),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),'*-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFitDivBg(a,:)),PopData.SigFit(a,:),'-','Color',colorMat{5}(a,:))
        ylabel('SR pulse-bg')
        xlabel('log(pulse/Bg)')    
    end
end


figure % tuning curves (PID)
subplot(1,6,1) % (pulse-rest vs absolute concentration)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.PidConcentrationRangePlusBg(a,:)),PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),'*-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFitDivBg(a,:)),PopData.SigFit(a,:),'-','Color',colorMat{5}(a,:))
        ylabel('SR pulse-rest')
        xlabel('log(pulse+Bg)')    
    end
end

subplot(1,6,2) % (pulse-bg vs absolute concentration)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.PidConcentrationRangePlusBg(a,:)),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),'o-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFit),PopData.SigFit(a,:),'Color',colorMat{5}(a,:),'LineWidth',1)
        ylabel('SR pulse-bg')
        xlabel('log(pulse+bg)')
    end
end

subplot(1,6,3)  % (pulse-rest vs pulse concentration)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.PidConcentrationRange),PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),'*-','Color',colorMat{5}(a,:))
        hold on
        ylabel('SR pulse-rest')
        xlabel('log(pulse)')    
    end
end

subplot(1,6,4) % (pulse-bg vs pulse concentration)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.PidConcentrationRange),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),'o-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFit),PopData.SigFit(a,:),'Color',colorMat{5}(a,:),'LineWidth',1)
        ylabel('SR pulse-bg')
        xlabel('log(pulse)')
    end
end

subplot(1,6,5)  % contrast axis (pulse-rest vs contrast)
for a = 2:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinRest_mean{a}) ; % if there are DB at this background
        errorbar(log10(PopData.PidConcentrationRangeDivBg(a,:)),PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon(a,:),'*-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFitDivBg(a,:)),PopData.SigFit(a,:),'-','Color',colorMat{5}(a,:))
        ylabel('SR pulse - rest')
        xlabel('log(pulse/Bg)')    
    end
end

subplot(1,6,6)  % contrast axis (pulse-bg vs contrast)
for a = 2:length(PopData.BgConcentrationRange) ; % for each background but control
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background   
        errorbar(log10(PopData.PidConcentrationRangeDivBg(a,:)),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon(a,:),'*-','Color',colorMat{5}(a,:))
        hold on
        %plot(log10(PopData.ConcentrationRangeFitDivBg(a,:)),PopData.SigFit(a,:),'-','Color',colorMat{5}(a,:))
        ylabel('SR pulse-bg')
        xlabel('log(pulse/Bg)')    
    end
end


figure % sigmoid fits as function of background with persect contrast adaptation fits
plot(log10(PopData.BgConcentrationRangeNum), MidPoint,'*')
hold on
plot(PopData.BgConcentrationRangeNumLogFit,PopData.HockeyStickfit,'r','LineWidth',1)
text(.1,.8,['slope=',num2str(PopData.HockeyStick_FitCoefs(3))],'Units','Normalized')
text(.1,.7,['infection pnt=',num2str(PopData.HockeyStick_FitCoefs(1))],'Units','Normalized')

figure % background transient adaptation as a function of bg

figure % background transient adaptation vs. pulse adaptation
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        if nansum(PopData.NumTrials{a}(:,b))>0 ; % if there is data at this background and pulse concentration
            subplot(1,length(PopData.ConcentrationRange),b)
            plot(PopData.SRbgMeanMinRest_DivBgTrans_Mean{a}(:,b),PopData.SRpulseMeanMinBg_mean_DivCon{a}(:,b),'*','Color',colorMat{5}(a,:)) ;
            hold on
        end
    end
end

figure % background transient decrease vs pulse adaptation
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        if nansum(PopData.NumTrials{a}(:,b))>0 ; % if there is data at this background and pulse concentration
            subplot(1,length(PopData.ConcentrationRange),b)
            plot(PopData.SRbgTransMeanMinRest_mean_MinBg{a}(:,b),PopData.SRpulseMeanMinBg_mean_DivCon{a}(:,b),'*','Color',colorMat{5}(a,:)) ;
            hold on
        end
    end
end

figure % background response vs adaptation factor(averaged across all pulse concentration)
for a = 2:length(PopData.BgConcentrationRange); % for each background 
    if ~isempty(PopData.SRbgMeanMinRest_ApMean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_ApMean{a}) ; % for each DB at that background
            plot(PopData.SRbgMeanMinRest_ApMean{a}(b),nanmean(PopData.SRpulseMeanMinBg_mean_DivCon{a}(b,:)),'*','Color',colorMat{5}(a,:))
            hold on
        end
    end
end

figure % background rate (avearge for across all pulses) vs. pulse adaptation 
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        if nansum(PopData.NumTrials{a}(:,b))>0 ; % if there is data at this background and pulse concentration
            subplot(1,length(PopData.ConcentrationRange),b)
            plot(PopData.SRbgMeanMinRest_ApMean{a}(:),PopData.SRpulseMeanMinBg_mean_DivCon{a}(:,b),'*','Color',colorMat{5}(a,:)) ;
            hold on
        end
    end
end

figure % background rate vs. pulse adaptation 
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        if nansum(PopData.NumTrials{a}(:,b))>0 ; % if there is data at this background and pulse concentration
            subplot(1,length(PopData.ConcentrationRange),b)
            plot(PopData.SRbgMeanMinRest_mean{a}(:,b),PopData.SRpulseMeanMinBg_mean_DivCon{a}(:,b),'*','Color',colorMat{5}(a,:)) ;
            hold on
        end
    end
end

figure % average pulse response
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        if nansum(PopData.NumTrials{a}(:,b))>0 ; % if there is data at this background and pulse concentration
            subplot(1,length(PopData.ConcentrationRange),b)
            plot(PopData.time(1:1+RspPnts),RspEnsemble_mean_AdbMean{a}(b,:),'Color',colorMat{5}(a,:))
            hold on
            plot(PopData.time(1:1+RspPnts),RspEnsemble_mean_AdbMean{a}(b,:)+RspEnsemble_mean_AdbSem{a}(b,:),':','Color',colorMat{5}(a,:))
            plot(PopData.time(1:1+RspPnts),RspEnsemble_mean_AdbMean{a}(b,:)-RspEnsemble_mean_AdbSem{a}(b,:),':','Color',colorMat{5}(a,:))
            axis tight
        end 
    end
end

figure % fraction change (from average tuning curves)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    plot(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_AdbWMean_FracChange(a,:),'o-','Color',colorMat{5}(a,:))
    hold on
end
  
figure % fractional change (plotted by pulse)
for a = 1:length(PopData.BgConcentrationRange) ; % for each background 
    for b = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
        if nansum(PopData.NumTrials{a}(:,b))>0 ; % if there is data at this background and pulse concentration
            subplot(1,length(PopData.ConcentrationRange),b)
            plot(PopData.SRbgMeanMinRest_mean{a}(:,b),PopData.SRpulseMeanMinBg_FracChange{a}(:,b),'*','Color',colorMat{5}(a,:)) ;
            hold on
        end
    end
end

%% ForIgor

% dose-response curves (x axes)
ForIgor.logOdorPulsePlusBg = log10(PopData.ConcentrationRangePlusBg) ;
ForIgor.logOdorPulse = log10(PopData.ConcentrationRange) ;
ForIgor.logOdorPulseDivBg = log10(PopData.ConcentrationRangeDivBg) ;

% ForIgor.logPidPulsePlusBg = log10(PopData.PidConcentrationRangePlusBg) ;
% ForIgor.logPidPulse = log10(PopData.PidConcentrationRange) ;
% ForIgor.logPidPulseDivBg = log10(PopData.PidConcentrationRangeDivBg) ;

% dose-response curves (y axes) 
ForIgor.DrPulseMinRest = PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon ;
ForIgor.DrPulseMinBg = PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon ;
ForIgor.DrPulseMinBgNotNorm = PopData.SRpulseMeanMinBg_mean_AdbWMean ;

ForIgor.DrPulseMinRestSem = PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon ;
ForIgor.DrPulseMinBgSem = PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon ;
ForIgor.DrPulseMinBgNotNormSem = PopData.SRpulseMeanMinBg_mean_AdbSem ;

% dose-response curves shifted so that ORN and PN points can be plotted against each other
if strcmp(PopData.AnalysisFlag,'ORNdata') ;
    ForIgor.DrPulseMinRestShifted = [nan(length(PopData.BgConcentrationRange),1),PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon] ;
    ForIgor.DrPulseMinBgShifted = [nan(length(PopData.BgConcentrationRange),1),PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon] ;

    ForIgor.DrPulseMinRestSemShifted = [nan(length(PopData.BgConcentrationRange),1),PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon] ;
    ForIgor.DrPulseMinBgSemShifted = [nan(length(PopData.BgConcentrationRange),1),PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon] ;

elseif strcmp(PopData.AnalysisFlag,'PNdata') ;
    ForIgor.DrPulseMinRestShifted = [PopData.SRpulseMeanMinRest_mean_DivCon_AdbWMean_TimCon,nan(length(PopData.BgConcentrationRange),1)] ;
    ForIgor.DrPulseMinBgShifted = [PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean_TimCon,nan(length(PopData.BgConcentrationRange),1)] ;

    ForIgor.DrPulseMinRestSemShifted = [PopData.SRpulseMeanMinRest_mean_DivCon_AdbSem_TimCon,nan(length(PopData.BgConcentrationRange),1)] ;
    ForIgor.DrPulseMinBgSemShifted = [PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem_TimCon,nan(length(PopData.BgConcentrationRange),1)] ;
end

% fractional change
% ForIgor.FracChange = PopData.SRpulseMeanMinRest_AdbWMean_FracChange ;
% ForIgor.FracChangeSem = PopData.SRpulseMeanMinRest_AdbSem_FracChange ;

% psth (all cell means) 
ForIgor.PsthMinRest0 = PopData.PsthMinRest_mean_AdMean{1} ;
ForIgor.PsthMinBg0 = PopData.PsthMinBg_mean_AdMean{1} ;

ForIgor.PsthMinRest7 = PopData.PsthMinRest_mean_AdMean{2} ;
ForIgor.PsthMinBg7 = PopData.PsthMinBg_mean_AdMean{2} ;

ForIgor.PsthMinRest6 = PopData.PsthMinRest_mean_AdMean{3} ;
ForIgor.PsthMinBg6 = PopData.PsthMinBg_mean_AdMean{3} ;

ForIgor.PsthMinRest5 = PopData.PsthMinRest_mean_AdMean{4} ;
ForIgor.PsthMinBg5 = PopData.PsthMinBg_mean_AdMean{4} ;

ForIgor.PsthMinRest4 = PopData.PsthMinRest_mean_AdMean{5} ;
ForIgor.PsthMinBg4 = PopData.PsthMinBg_mean_AdMean{5} ;

% psth (all cell means) normalized by bg transient
ForIgor.PsthMinRestTransNorm0 = PopData.PsthMinRest_mean_AdMean_TransNorm{1} ;
ForIgor.PsthMinRestTransNorm7 = PopData.PsthMinRest_mean_AdMean_TransNorm{2} ;
ForIgor.PsthMinRestTransNorm6 = PopData.PsthMinRest_mean_AdMean_TransNorm{3} ;
ForIgor.PsthMinRestTransNorm5 = PopData.PsthMinRest_mean_AdMean_TransNorm{4} ;
ForIgor.PsthMinRestTransNorm4 = PopData.PsthMinRest_mean_AdMean_TransNorm{5} ;

% psth time
ForIgor.PsthTime = PopData.time ;

% bg rate 
ForIgor.bgMean = PopData.SRbgMeanMinRest_ApMean_AdbWMean ;
ForIgor.bgMeanSem = PopData.SRbgMeanMinRest_ApMean_AdbSem ;

% bg variance
ForIgor.bgStd =PopData.SRbgStd_ApMean_DivCon_AdbWMean_TimCon ;
ForIgor.bgStdSem =PopData.SRbgStd_ApMean_DivCon_AdbSem_TimCon ;

% bg SNR
ForIgor.bgSnr = PopData.SRbgSNR_mean_AdbWMean ;
ForIgor.bgSnrSem = PopData.SRbgSNR_mean_AdbSem ;

% peak comparison for methods
ForIgor.PeakMinBg0 = PopData.SRpulseMeanMinBg_mean{1} ;
ForIgor.PeakMinBg7 = PopData.SRpulseMeanMinBg_mean{2} ;
ForIgor.PeakMinBg6 = PopData.SRpulseMeanMinBg_mean{3} ;
ForIgor.PeakMinBg5 = PopData.SRpulseMeanMinBg_mean{4} ;
ForIgor.PeakMinBg4 = PopData.SRpulseMeanMinBg_mean{5} ;

% adaptation factors for methods
ForIgor.logOdorPulse = log10(PopData.ConcentrationRange) ;

ForIgor.PeakMinBgAdFact0 = PopData.SRpulseMeanMinBg_mean_DivCon{1} ;
ForIgor.PeakMinBgAdFact7 = PopData.SRpulseMeanMinBg_mean_DivCon{2} ;
ForIgor.PeakMinBgAdFact6 = PopData.SRpulseMeanMinBg_mean_DivCon{3} ;
ForIgor.PeakMinBgAdFact5 = PopData.SRpulseMeanMinBg_mean_DivCon{4} ;
ForIgor.PeakMinBgAdFact4 = PopData.SRpulseMeanMinBg_mean_DivCon{5} ;

ForIgor.PeakMinBgAdFact0sem = PopData.SRpulseMeanMinBg_sem_DivCon{1} ;
ForIgor.PeakMinBgAdFact7sem = PopData.SRpulseMeanMinBg_sem_DivCon{2} ;
ForIgor.PeakMinBgAdFact6sem = PopData.SRpulseMeanMinBg_sem_DivCon{3} ;
ForIgor.PeakMinBgAdFact5sem = PopData.SRpulseMeanMinBg_sem_DivCon{4} ;
ForIgor.PeakMinBgAdFact4sem = PopData.SRpulseMeanMinBg_sem_DivCon{5} ;

ForIgor.PeakMinBgAdFactWmeans = PopData.SRpulseMeanMinBg_mean_DivCon_AdbWMean ; 
ForIgor.PeakMinBgAdFactWmeansSem = PopData.SRpulseMeanMinBg_mean_DivCon_AdbSem ;


% stimulus schematic
StimSchematic = zeros(5,length(PopData.time)) ;
StimSchematic(:,iopb:iope) = .8 ;
for a=1:5 ;
    StimSchematic(a,ibpb:ibpe) = StimSchematic(a,ibpb:ibpe)+a-1 ;
end
ForIgor.StimSchematic = StimSchematic ;





