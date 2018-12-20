function ForIgor = PNadaptationAnalysisJ(Input,A)

% adaptation of "PNadaptationAnalysisI" 
% focuses on psth peak, different measures of bg and rest activity, determines discriminability 
% Input.id1 should be a cell array with trial numbers grouped according to vector specified in Input.id2.
% JC 6/27/13

for DB = 1:length(A) ; % for each cell

    % parameters
    sampleRate = 10000 ; % (hz) temp hard coded - should be saved in file 
    driftCheckTime = 0.25 ; %(sec) time at begining and end which current injected is inspected for changes
    bgTransTime = 5 ; %(sec) time after background odor starts to avoid transient response
    PsthBinTime = 0.1 ; % (sec) time width of psth bin
    OdorRspTime = .05 ; % (sec) time before and after peak of mean psth which odor response is assessed

    spikeDetectionParameters.absRefTime = 0.002 ; % sec
    spikeDetectionParameters.minRise = 0.3 ; % mV
    spikeDetectionParameters.minFall = 0.15 ; % mV
    spikeDetectionParameters.filterOrder = 1 ;
    spikeDetectionParameters.lpfCutOff = 4000 ; % hz
    spikeDetectionParameters.minRiseTime = 0.001 ; % sec
    spikeDetectionParameters.minFallTime = 0.001 ; % sec
    spikeDetectionParameters.maxPlatTime = 0.001 ; % sec 

    spikeDataPath = ['Z:/Cafaro Documents/Analysis/DetectedSpikes/'] ;

    id1 = 'OdorRsp' ;
    id2 = 'OdorConcentration' ;
    id3 = 'BgConcentration' ;

    % load data in matricies
    rootdir = ['Z:\Cafaro Data Backup\', Input(A(DB)).cellname(1:6),'Data'];

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

                temp = load([rootdir,'\',Input(A(DB)).cellname,'\','voltage_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                vData{a}{b}(loopNum,:) = temp.voltage ;

                temp = load([rootdir,'\',Input(A(DB)).cellname,'\','current_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                iData{a}{b}(loopNum,:) = temp.current ;

                temp = load([rootdir,'\',Input(A(DB)).cellname,'\','Ao0_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                ao0Data{a}{b}(loopNum,:) = temp.Ao0 ; % odor

                temp = load([rootdir,'\',Input(A(DB)).cellname,'\','Ao1_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                ao1Data{a}{b}(loopNum,:) = temp.Ao1 ; % background valve   

                temp = load([rootdir,'\',Input(A(DB)).cellname,'\','TrigTime_',Input(A(DB)).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
                tData{a}{b}(loopNum) = temp.Trigtime ;
            end
        end
    end

    % time vector
    time = [1:length(vData{1}{1})]/sampleRate ;

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

    % make odor valve pulse binary
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

    % make sure odor pulse was the same time
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

    % check that input current is not changing substantially during any of the trials
    driftCheckPnts = driftCheckTime*sampleRate ;
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

    % detect spikes in voltage data
    try PreviousDetect = load([spikeDataPath,'cell',num2str(A(DB))]) ; 
        if ~isequal(PreviousDetect.spikeDetectionParameters,spikeDetectionParameters) ; % if your previous spike parameters were not the same
            disp('spike detection parameters assumed are not as specified above')
        end
        spikePnt = PreviousDetect.spikePnt ;

    catch

        for a = 1:NumBackgrounds ; % for each background
            for b = 1:NumConcentrations ; % for concentration
                if NumTrials(a,b)>0 ;
                    [TempSpikePnt,SpikeData,NonSpikeData] = spikeFinder(vData{a}{b},sampleRate,spikeDetectionParameters) ;
                    spikePnt{a}{b}= TempSpikePnt ;
                end
            end
        end

        save([spikeDataPath,'cell',num2str(A(DB))], 'spikePnt','spikeDetectionParameters') ; % save spike data
    end

    % spike trains and psth
    PsthBinPnts = PsthBinTime*sampleRate ;
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SpikeTrain{a}{b}(c,:)= zeros(1,length(vData{a}{b}(c,:))) ;
                    SpikeTrain{a}{b}(c,spikePnt{a}{b}{c}) = 1 ;
                    SpikeTrainSmooth{a}{b}(c,:) = smooth(SpikeTrain{a}{b}(c,:),PsthBinPnts) ;
                    Psth{a}{b}(c,:) = SpikeTrainSmooth{a}{b}(c,:)*sampleRate ;
                end
                Psth_mean{a}{b} = mean(Psth{a}{b},1);
            end
        end
    end

    % index of background odor and odor pulse
    iopb = find(ao0DataB{2}{1}(1,:)~=0,1,'first')-1 ; % odor pulse begining
    iope = find(ao0DataB{2}{1}(1,:)~=0,1,'last')-1 ; % odor pulse ending

    ibpb = find(ao1DataB{2}{1}(1,:)~=0,1,'first')-1 ; % background pulse beginging
    ibpe = find(ao1DataB{2}{1}(1,:)~=0,1,'last')-1 ; % background pulse end
    
    % odor response time
    OdorRspPnts = OdorRspTime*sampleRate ;
    
    for a = 1:NumBackgrounds ;
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;         
                [m,mi] = max(Psth_mean{a}{b}(iopb:iopb+sampleRate)) ; % max point of mean psth within a second of odor pulse onset
                SRiorb(a,b) = mi-1+iopb - OdorRspPnts ; % point of odor pulse begining
                SRiore(a,b) = mi-1+iopb +  OdorRspPnts ; % point of odor pulse end
            end
        end
    end

    % spontaneous spike rate
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration   
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRrestMean{a}{b}(c) = mean(Psth{a}{b}(c,ceil(PsthBinPnts/2):ibpb)) ; % (avoid bining artifacts: background pulse begining)
                end
                SRrestMean_mean(a,b) = mean(SRrestMean{a}{b}) ;
                SRrestMean_sem(a,b) = std(SRrestMean{a}{b})/sqrt(NumTrials(a,b)) ;
            else
                SRrestMean_mean(a,b) = nan ;
                SRrestMean_sem(a,b) = nan ;
            end
        end
    end

    % assess background spike rate (during background odor after transient response)
    bgTransPnts = bgTransTime*sampleRate ;
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

    % odor pulse response spike rate
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

    % pulse minus rest
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMeanMinRest{a}{b} = SRpulseMean{a}{b} - SRrestMean{a}{b} ;
                SRpulseMeanMinRest_mean(a,b) = mean(SRpulseMeanMinRest{a}{b}) ;
                SRpulseMeanMinRest_std(a,b) = std(SRpulseMeanMinRest{a}{b}) ;
                SRpulseMeanMinRest_sem(a,b) = SRpulseMeanMinRest_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRpulseMeanMinRest_mean(a,b) = nan ;
                SRpulseMeanMinRest_std(a,b) = nan ;
                SRpulseMeanMinRest_sem(a,b) = nan ;
            end
        end
    end
    
    % backgound minus rest
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

    % pulse minus bg
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMeanMinBg{a}{b} = SRpulseMean{a}{b} - SRbgMean{a}{b} ;
                SRpulseMeanMinBg_mean(a,b) = mean(SRpulseMeanMinBg{a}{b}) ;
                SRpulseMeanMinBg_std(a,b) = std(SRpulseMeanMinBg{a}{b}) ;
                SRpulseMeanMinBg_sem(a,b) = SRpulseMeanMinBg_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRpulseMeanMinBg_mean(a,b) = nan ;
                SRpulseMeanMinBg_std(a,b) = nan ;
                SRpulseMeanMinBg_sem(a,b) = nan ;
            end
        end
    end
    
    % bg minus bg
    SRbgMeanMinBg_mean_DBwMean = zeros(1,NumBackgrounds) ; % x-x=0
    
    % rest and bg activity subtracted psth vectors
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    PsthMinRest{a}{b}(c,:) = Psth{a}{b}(c,:) - SRrestMean{a}{b}(c) ;
                    PsthMinBg{a}{b}(c,:) = Psth{a}{b}(c,:) - SRbgMean{a}{b}(c) ;
                end
            PsthMinRest_mean{a}{b} = mean(PsthMinRest{a}{b},1) ;
            PsthMinBg_mean{a}{b} = mean(PsthMinBg{a}{b},1) ;
            end
        end
    end

    % repository status
    temp = getGitInfo ; 
    RepVer = temp.hash ; % repository version 

    [tempS,tempR] = system('git status') ;
    if length(tempR)==62 ;
        RepStat = ['GitHub up to date ',RepVer] ; % if no unsynced files
    else
        RepStat = ['GitHub not up to date ',RepVer] ;
    end

    % figures
    Conc(1,:) = [.6,.5,.4,0,0,0,0] ;
    Conc(2,:) = [1,1,1,1,.6,.5,.4] ;
    Conc(3,:) = [.6,.5,.4,.3,.2,.1,0] ;
    for a=1:NumConcentrations ; % each concentration is a matrix and each background is a row within that matrix
        colorMat{a} = [Conc(3,a),Conc(3,a),Conc(3,a); Conc(2,a),Conc(1,a),Conc(1,a); Conc(1,a),Conc(1,a),Conc(2,a);... 
            Conc(1,a),Conc(2,a),Conc(1,a);Conc(2,a),Conc(2,a),Conc(1,a);...
            Conc(2,a),Conc(1,a),Conc(2,a);Conc(1,a),Conc(2,a),Conc(2,a)] ;
    end
        
    % spike detection
%     figure 
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             for c = 1:NumTrials(a,b) ;
%                 plot(time,vData{a}{b}(c,:)) 
%                 hold on
%                 plot(time(spikePnt{a}{b}{c}),vData{a}{b}(c,spikePnt{a}{b}{c}),'r*')
%                 title(num2str(odorRspTrials{a}{b}(c)))
%                 hold off
%                 pause
%             end
%         end
%     end
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
%     
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
%     % mean psth for different concentrations and backgrounds and time (rest subtracted data)
%     figure 
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(NumConcentrations+2,1,b)
%                 plot(time,PsthMinRest_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%                 hold on   
%                 axis tight
%             end
%         end
%     end
%     
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(NumConcentrations+2,1,NumConcentrations+2)
%                 plot(tDataN{a}{b},SRbgMeanMinRest{a}{b},'.','Color',colorMat{b}(a,:))
%                 hold on
%                 axis tight
%                 xlabel('trig time (sec)')
%                 ylabel('bg sr - rest sr (hz)')
%                 
%                 subplot(NumConcentrations+2,1,NumConcentrations+1)
%                 plot(tDataN{a}{b},SRpulseMeanMinRest{a}{b},'o','Color',colorMat{b}(a,:))
%                 hold on
%                 axis tight
%                 xlabel('trig time (sec)')
%                 ylabel('pulse sr - rest sr (hz)')
%             end
%         end
%     end
%     
%      % mean psth for different concentrations and backgrounds and time (bg subtracted data)
%     figure 
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(NumConcentrations+1,1,b)
%                 plot(time,PsthMinBg_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%                 hold on   
%                 axis tight
%             end
%         end
%     end
%     
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(NumConcentrations+1,1,NumConcentrations+1)
%                 plot(tDataN{a}{b},SRpulseMeanMinBg{a}{b},'o','Color',colorMat{b}(a,:))
%                 hold on
%                 xlabel('trig time (sec)')
%                 ylabel('pulse sr - bg sr (hz)')
%                 axis tight
%             end
%         end
%     end
% 
%     % comparing rest and background data across concentrations
%     figure 
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(1,2,1)
%                 plot(time,Psth_mean{a}{b},'Color',colorMat{b}(a,:))            
%                 hold on     
%     
%                 subplot(1,2,2)
%                 plot(time,PsthMinRest_mean{a}{b},'Color',colorMat{b}(a,:))            
%                 hold on                  
%             end
%         end
%     end
%    
%     % mean and std spike rate as a function of background
%     figure
%     for a = 1:NumBackgrounds ;
%         subplot(2,1,1)
%         errorbar(a, SRbgMeanMinRest_ApMean(a),SRbgMeanMinRest_ApSem(a),SRbgMeanMinRest_ApSem(a),'Color',colorMat{NumConcentrations}(a,:))
%         hold on
%         xlabel('background concentration')
%         ylabel('mean sr bg-rest')
% 
%         subplot(2,1,2)
%         errorbar(a, SRbgStd_ApMean(a),SRbgStd_ApSem(a),SRbgStd_ApSem(a),'Color',colorMat{NumConcentrations}(a,:))
%         hold on
%         xlabel('background concentration')
%         ylabel('std sr bg')
%     end
%     
%    % mean and std spike rate as a function of pulse and background
%     figure
%     for a = 1:NumBackgrounds ;
%         subplot(2,1,1)
%         errorbar(log10(Concentrations), SRpulseMeanMinBg_mean(a,:),SRpulseMeanMinBg_sem(a,:),SRpulseMeanMinBg_sem(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%         xlabel('pulse concentration')
%         ylabel('mean sr pulse-bg')
% 
%         subplot(2,1,2)
%         plot(log10(Concentrations), SRpulseMeanMinBg_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%         xlabel('background concentration')
%         ylabel('std sr bg')  
%     end
% 
%     %  spikes as a function of concentration
%     figure
%     subplot(2,1,1)
%     for a = 1:NumBackgrounds ;
%         errorbar([-10,log10(Concentrations)],[SRbgMeanMinRest_mean_DBwMean(a), SRpulseMeanMinRest_mean(a,:)],[SRbgStd_mean_DBwMean(a),SRpulseMean_std(a,:)],[SRbgStd_mean_DBwMean(a),SRpulseMean_std(a,:)],'*-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%     end
%     ylabel('spike rate')
%     xlabel('log concentration')
% 
%     subplot(2,1,1)
%     for a = 1:NumBackgrounds ;
%          errorbar([-10,log10(Concentrations)],[SRbgMeanMinBg_mean_DBwMean(a), SRpulseMeanMinBg_mean(a,:)],[SRbgStd_mean_DBwMean(a),SRpulseMeanMinBg_std(a,:)],[SRbgStd_mean_DBwMean(a),SRpulseMeanMinBg_std(a,:)],'*-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%     end
%     ylabel('pulse min bg spike rate')
%     xlabel('log concentration')
% 


    % temp results for population analysis
    
    % range of background and pulse concentrations assessed across
    % population and indicies to compartmentalize this DB correctly
    PopData.ConcentrationRange = 10.^[-8:-1] ;    
    PopData.BgConcentrationRange = {'0','Log8','Log7','Log6','Log5','Log4','Log3','wash'} ;
    PopData.BgConcentrationRangeNum = 10.^[-9:-3] ; % numerical version, control changed to 10^-9 and wash not included
    
    [c,PopPulsei] = intersect(PopData.ConcentrationRange,Concentrations) ;
    
    for a = 1:NumBackgrounds ;
        for b = 1:length(PopData.BgConcentrationRange) ;
            if strcmp(BgConcentrations{a},PopData.BgConcentrationRange{b}) ;
                PopBgi(a) = b ;
            end
        end
    end
    
    % {bg}(cell,pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            for b=1:length(PopData.ConcentrationRange) ;
                PopData.NumTrials{a} = nan(length(A),length(PopData.ConcentrationRange));
            end
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            PopData.NumTrials{PopBgi(a)}(DB,PopPulsei) = NumTrials(a,b) ;
        end
    end 
    
    % {bg}{pulse}(cell,time)  
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            for b=1:length(PopData.ConcentrationRange) ;
                PopData.Psth_mean{a}{b} = nan(length(A),length(time)) ;
                PopData.PsthMinRest_mean{a}{b} = nan(length(A),length(time)) ;
                PopData.PsthMinBg_mean{a}{b} = nan(length(A),length(time)) ;
            end
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                PopData.Psth_mean{PopBgi(a)}{PopPulsei(b)}(DB,:) = Psth_mean{a}{b} ;
                PopData.PsthMinRest_mean{PopBgi(a)}{PopPulsei(b)}(DB,:) = PsthMinRest_mean{a}{b} ;
                PopData.PsthMinBg_mean{PopBgi(a)}{PopPulsei(b)}(DB,:) = PsthMinBg_mean{a}{b} ;
            end
        end
    end 
    
    % {bg}(cell, pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRpulseMeanMinBg_mean{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinBg_sem{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMeanMinBg_std{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        PopData.SRpulseMeanMinBg_mean{PopBgi(a)}(DB,PopPulsei) = SRpulseMeanMinBg_mean(a,:) ;
        PopData.SRpulseMeanMinBg_sem{PopBgi(a)}(DB,PopPulsei) = SRpulseMeanMinBg_sem(a,:) ;
        PopData.SRpulseMeanMinBg_std{PopBgi(a)}(DB,PopPulsei) = SRpulseMeanMinBg_std(a,:) ;
    end
    
    % {bg} (cell)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRbgMeanMinRest_mean{a} = nan(1,length(A)) ;
            PopData.SRbgMeanMinRest_sem{a} = nan(1,length(A)) ;
            PopData.SRbgStd_mean{a} = nan(1,length(A)) ;
            PopData.SRbgStd_sem{a} = nan(1,length(A)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        PopData.SRbgMeanMinRest_mean{PopBgi(a)}(DB) = SRbgMeanMinRest_ApMean(a) ;
        PopData.SRbgMeanMinRest_sem{PopBgi(a)}(DB) = SRbgMeanMinRest_ApSem(a) ;
        PopData.SRbgStd_mean{PopBgi(a)}(DB) = SRbgStd_ApMean(a) ;
        PopData.SRbgStd_sem{PopBgi(a)}(DB) = SRbgStd_ApSem(a) ;
    end
    
    % {bg}{pulse}(cell,[contPulse,bgPulse])
%     if DB==1 ;
%         for a=2:length(PopData.BgConcentrationRange) ;
%             for b=1:length(PopData.ConcentrationRange) ; 
%                 PopData.SRpulseMinBg_mean_paired{a}{b} = nan(2,length(A)) ;
%                 PopData.SRpulseMinBg_sem_paired{a}{b} = nan(2,length(A))  ; 
%             end
%         end
%     end
%     for a = 2:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration 
%             if NumTrials(a,b)>0 ;
%                 PopData.SRpulseMinBg_mean_paired{PopBgi(a)}{PopPulsei(b)}(:,DB) = [SRpulseMinBg_mean(1,b);SRpulseMinBg_mean(a,b)] ;
%                 PopData.SRpulseMinBg_sem_paired{PopBgi(a)}{PopPulsei(b)}(:,DB) = [SRpulseMinBg_sem(1,b);SRpulseMinBg_sem(a,b)] ;  
%             end    
%         end
%     end
    disp(DB)
    % clear all variables except ForIgor and PopData structures
    clearvars -except Input A ForIgor PopData  
 end % Data block loop

% population analysis

% bg spike rate mean
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRbgMeanMinRest_mean{a}) ; % if there are DB at this background
        PopData.SRbgMeanMinRest_mean_AdbWMean(a) = nansum(PopData.SRbgMeanMinRest_mean{a}./PopData.SRbgMeanMinRest_sem{a})./nansum(1./PopData.SRbgMeanMinRest_sem{a}) ; % wieghted mean across data blocks
        PopData.SRbgMeanMinRest_mean_AdbSem(a) = nanstd(PopData.SRbgMeanMinRest_mean{a})./sqrt(sum(~isnan(PopData.SRbgMeanMinRest_mean{a}))) ; 
    end
end

% bg spike rate variability
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRbgMeanMinRest_mean{a}) ; % if there are DB at this background
        PopData.SRbgStd_mean_AdbWMean(a) = nansum(PopData.SRbgStd_mean{a}./PopData.SRbgStd_sem{a})./nansum(1./PopData.SRbgStd_sem{a}) ; % wieghted mean across data blocks
        PopData.SRbgStd_mean_AdbSem(a) = nanstd(PopData.SRbgStd_mean{a})./sqrt(sum(~isnan(PopData.SRbgStd_mean{a}))) ; 
    end
end

% pulse spike rate mean
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinBg_mean_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinBg_mean{a}./PopData.SRpulseMeanMinBg_sem{a},1)./nansum(1./PopData.SRpulseMeanMinBg_sem{a},1) ; % wiehgted mean across db
        PopData.SRpulseMeanMinBg_mean_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinBg_mean{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinBg_mean{a}),1)) ; 
    end
end

% pulse spike rate variability
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        PopData.SRpulseMeanMinBg_std_AdbWMean(a,:) = nansum(PopData.SRpulseMeanMinBg_std{a}.*sqrt(PopData.NumTrials{a}),1)./nansum(sqrt(PopData.NumTrials{a}),1) ; % wiehgted mean across db
        PopData.SRpulseMeanMinBg_std_AdbSem(a,:) = nanstd(PopData.SRpulseMeanMinBg_std{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMeanMinBg_std{a}),1)) ; 
    end
end

% pulse contrast axis
for a = 1:length(PopData.BgConcentrationRangeNum) ; % for each background
    PopData.ConcentrationRangeDivBg(a,:) = PopData.ConcentrationRange/PopData.BgConcentrationRangeNum(a) ;
end

% mean std of pulse and background responses
AssumedStd = 10 ; % to keep it simple I assume an std for all pulses and backgrounds is the same

% discriminability of background and pulse
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
     if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
         PopData.DfromBg(a,:) = PopData.SRpulseMeanMinBg_mean_AdbWMean(a,:)/AssumedStd ;
     end
end
    
% discriminability of pulses from nearby pulses
for a = 1:length(PopData.BgConcentrationRange) ; % for each background
     if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
         PopData.DfromPulse(a,:) = abs(diff([0,PopData.SRpulseMeanMinBg_mean_AdbWMean(a,:)]))/AssumedStd ;
     end
end

 
% figures
Conc(1,:) = [.7,.6,.5,.4,0,0,0,0] ;
Conc(2,:) = [1,1,1,1,.7,.6,.5,.4] ;
Conc(3,:) = [.7,.6,.5,.4,.3,.2,.1,0] ;
for a=1:length(PopData.ConcentrationRange) ; % each concentration is a matrix and each background is a row within that matrix
    colorMat{a} = [Conc(3,a),Conc(3,a),Conc(3,a); Conc(2,a),Conc(1,a),Conc(1,a); Conc(1,a),Conc(1,a),Conc(2,a);... 
        Conc(1,a),Conc(2,a),Conc(1,a);Conc(2,a),Conc(2,a),Conc(1,a);...
        Conc(2,a),Conc(1,a),Conc(2,a);Conc(1,a),Conc(2,a),Conc(2,a);...
        Conc(3,a),Conc(3,a),Conc(1,a)] ;
end

figure
for a = 1:length(PopData.BgConcentrationRange)-1 ; % for each background (but 'wash')
    if ~isempty(PopData.SRbgMeanMinRest_mean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_mean{a}) ; % for each DB at that background

            subplot(2,1,1)
            errorbar(a,PopData.SRbgMeanMinRest_mean{a}(b),PopData.SRbgMeanMinRest_sem{a}(b),PopData.SRbgMeanMinRest_sem{a}(b),'Color',colorMat{5}(a,:))
            hold on
            plot(a,PopData.SRbgMeanMinRest_mean_AdbWMean(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            ylabel('SRbgMeanMinRest')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)-1])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end-1))

            subplot(2,1,2)
            errorbar(a,PopData.SRbgMeanMinRest_mean_AdbWMean(a),PopData.SRbgMeanMinRest_mean_AdbSem(a),PopData.SRbgMeanMinRest_mean_AdbSem(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('SRbgMeanMinRest')
            xlabel('background')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)-1])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end-1))
        end
    end
end

figure
for a = 1:length(PopData.BgConcentrationRange)-1 ; % for each background (but 'wash')
    if ~isempty(PopData.SRbgMeanMinRest_mean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_mean{a}) ; % for each DB at that background

            subplot(2,1,1)
            errorbar(a,PopData.SRbgStd_mean{a}(b),PopData.SRbgStd_sem{a}(b),PopData.SRbgStd_sem{a}(b),'*','Color',colorMat{5}(a,:))
            hold on
            plot(a,PopData.SRbgStd_mean_AdbWMean(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            ylabel('SRbgStd')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)-1])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end-1))

            subplot(2,1,2)
            errorbar(a,PopData.SRbgStd_mean_AdbWMean(a),PopData.SRbgStd_mean_AdbSem(a),PopData.SRbgStd_mean_AdbSem(a),'o','LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('SRbgStd')
            xlabel('background')
            set(gca,'Xtick',[1:length(PopData.BgConcentrationRange)-1])
            set(gca,'XtickLabel',PopData.BgConcentrationRange(1:end-1))
        end
    end
end

figure
for a = 1:length(PopData.BgConcentrationRange)-1 ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_mean{a}) ; % for each DB at that background
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

figure
for a = 1:length(PopData.BgConcentrationRange)-1 ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_std{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_mean{a}) ; % for each DB at that background
            subplot(2,1,1)
            plot(log10(PopData.ConcentrationRange),PopData.SRpulseMeanMinBg_std{a}(b,:),'Color',colorMat{5}(a,:))
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

figure
for a = 1:length(PopData.BgConcentrationRange)-1 ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_std{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_mean{a}) ; % for each DB at that background
            subplot(2,1,1)
            plot(log10(PopData.ConcentrationRange),PopData.DfromBg(a,:),'Color',colorMat{5}(a,:))
            hold on
            plot(log10(PopData.ConcentrationRange),ones(1,length(PopData.ConcentrationRange)),'--k')
            ylabel('"D" from bg (std)')
            xlabel('log(pulse)')
            
            subplot(2,1,2)
            plot(log10(PopData.ConcentrationRange),PopData.DfromPulse(a,:),'Color',colorMat{5}(a,:))
            plot(log10(PopData.ConcentrationRange),ones(1,length(PopData.ConcentrationRange)),'--k')
            hold on
            ylabel('"D" from previous pulse (std)')
            xlabel('log(pulse)')
            
        end
    end
end
 
figure
for a = 1:length(PopData.BgConcentrationRange)-1 ; % for each background
    if ~isempty(PopData.SRpulseMeanMinBg_mean{a}) ; % if there are DB at this background
        for b = 1:length(PopData.SRbgMeanMinRest_mean{a}) ; % for each DB at that background
            errorbar(log10(PopData.ConcentrationRangeDivBg(a,:)),PopData.SRpulseMeanMinBg_mean_AdbWMean(a,:),PopData.SRpulseMeanMinBg_mean_AdbSem(a,:),PopData.SRpulseMeanMinBg_mean_AdbSem(a,:),'LineWidth',2,'Color',colorMat{5}(a,:))
            hold on
            ylabel('SRpulseMeanMinBg')
            xlabel('log(pulse/Bg)')    
        end
    end
end

