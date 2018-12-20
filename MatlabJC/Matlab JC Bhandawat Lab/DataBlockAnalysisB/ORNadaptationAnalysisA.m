function ForIgor = ORNadaptationAnalysisA(Input,A)

% this function will analyze ORN pulse responses +/- a background odor that
% is triggered during the trial
% Modified from PNadaptationAnalysisI.m
% JC 4/25/13

for DB = 1:length(A) ; % for each cell

    % parameters
    sampleRate = 10000 ; % (hz) temp hard coded - should be saved in file 
    driftCheckTime = 0.25 ; %(sec) time at begining and end which current injected is inspected for changes
    bgTransTime = 5 ; %(sec) time after background odor starts to avoid transient response
    PsthBinTime = 0.1 ; % (sec) time width of psth bin
    OdorRspTime = .1 ; % (sec) time before and after peak of mean psth which odor response is assessed

    spikeDetectionParameters.numClusters = 1 ;
    spikeDetectionParameters.numPCs = 2 ;
    spikeDetectionParameters.NegDiffThreshStd = 1.5 ; 
    SpikeThreshold = 10 ; % mv above base

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

    % detect spikes in voltage data
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
                        [spikePntGroup] = spikeSorter(vData{a}{b}(c,:), sampleRate, spikeDetectionParameters.NegDiffThreshStd, spikeDetectionParameters.numClusters, spikeDetectionParameters.numPCs, false) ;
                        spikePnt{a}{b}{c}= spikePntGroup{1} ;
                    end
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

    % index of current pulse and odor pulse
    iopb = find(ao0DataB{2}{1}(1,:)~=0,1,'first')-1 ; % odor pulse begining
    iope = find(ao0DataB{2}{1}(1,:)~=0,1,'last')-1 ; % odor pulse ending

    ibpb = find(ao1DataB{2}{1}(1,:)~=0,1,'first')-1 ; % background pulse beginging
    ibpe = find(ao1DataB{2}{1}(1,:)~=0,1,'last')-1 ; % background pulse end

    % odor response time
    for a = 1:NumBackgrounds ;
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;       
                [m,mi] = max(Psth_mean{a}{b}(iopb:iopb+sampleRate)) ; % max point of mean psth within a second of odor pulse onset
                SRiorb(a,b) = mi-1+iopb - OdorRspTime*sampleRate ; % point of odor pulse begining
                SRiore(a,b) = mi-1+iopb + OdorRspTime*sampleRate ; % point of odor pulse end

                [m,mi] = min(Psth_mean{a}{b}(iope:iope+sampleRate)) ; % min point of mean psth within a second of odor pulse offset
                SRiohb(a,b) = mi-1+iope - OdorRspTime*sampleRate ; % point of odor pulse begining
                SRiohe(a,b) = mi-1+iope + OdorRspTime*sampleRate ; % point of odor pulse end
            end
        end
    end

    % spontaneous spike rate
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration   
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRrest{a}{b}(c) = mean(Psth{a}{b}(c,1:ibpb)) ; %
                end
                SRrest_mean(a,b) = mean(SRrest{a}{b}) ;
                SRrest_std(a,b) = std(SRrest{a}{b}) ;
                SRrest_sem(a,b) = SRrest_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRrest_mean(a,b) = nan ;
                SRrest_std(a,b) = nan ; 
                SRrest_sem(a,b) = nan ;
            end
        end
    end

    % assess background transient psth
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRtbg{a}{b}(c) = max(Psth{a}{b}(c,ibpb:iopb)) ; % mV (start of background: begining of odor pulse)
                end
                SRtbg_mean(a,b) = mean(SRtbg{a}{b}) ;
                SRtbg_std(a,b) = std(SRtbg{a}{b}) ;
                SRtbg_sem(a,b) = SRtbg_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRtbg_mean(a,b) = nan ;
                SRtbg_std(a,b) = nan ;
                SRtbg_sem(a,b) = nan ;
            end
        end
    end

    % assess background spike rate (during background odor after transient response)
    bgTransPnts = bgTransTime*sampleRate ;
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRbg{a}{b}(c) = mean(Psth{a}{b}(c,ibpb+bgTransPnts:iopb)) ; % mV (start of background + transient time: begining of odor pulse)
                end
                SRbg_mean(a,b) = mean(SRbg{a}{b}) ;
                SRbg_std(a,b) = std(SRbg{a}{b}) ;
                SRbg_sem(a,b) = SRbg_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRbg_mean(a,b) = nan ;
                SRbg_std(a,b) = nan ;
                SRbg_sem(a,b) = nan ;
            end
        end
    end
    
    % assess background spike rate variability from psth
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRbgStd{a}{b}(c) = std(Psth{a}{b}(c,ibpb+bgTransPnts:iopb)) ; % mV (start of background + transient time: begining of odor pulse)
                end
                SRbgStd_mean(a,b) = mean(SRbgStd{a}{b}) ;
                SRbgStd_std(a,b) = std(SRbgStd{a}{b}) ;
            else
                SRbgStd_mean(a,b) = nan ;
                SRbgStd_std(a,b) = nan ;
            end
        end
    end    
    
    % assess background spike rate variability from ISI distribution
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                IsiBg{a}{b}=[] ;
                for c = 1:NumTrials(a,b) ;
                    IsiBg{a}{b} = [IsiBg{a}{b},diff(spikePnt{a}{b}{c}(spikePnt{a}{b}{c}>=ibpb+bgTransPnts & spikePnt{a}{b}{c}<=iopb))] ;
                end
                IsiBg_mean(a,b) = mean(IsiBg{a}{b});
                IsiBg_std(a,b) = std(IsiBg{a}{b}) ;
                IsiBg_cv(a,b) = IsiBg_std(a,b)/IsiBg_mean(a,b) ; 
            end
        end
    end

    % odor pulse response spike rate
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRpulse{a}{b}(c) = mean(Psth{a}{b}(c,SRiorb(a,b):SRiore(a,b))); %mV (odor response depol start: odor response depol end)
                    SRpulseHyp{a}{b}(c) = mean(Psth{a}{b}(c,SRiohb(a,b):SRiohe(a,b))); %mV (odor response hyperpol start: odor response hyperpol end)
                end

                SRpulse_mean(a,b) = mean(SRpulse{a}{b}) ;
                SRpulse_std(a,b) = std(SRpulse{a}{b}) ;
                SRpulse_sem(a,b) = SRpulse_std(a,b)/sqrt(NumTrials(a,b)) ;

                SRpulseHyp_mean(a,b) = mean(SRpulseHyp{a}{b}) ;
                SRpulseHyp_std(a,b) = std(SRpulseHyp{a}{b}) ; 
                SRpulseHyp_sem(a,b) = SRpulseHyp_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRpulse_mean(a,b) = nan ;
                SRpulse_std(a,b) = nan ;
                SRpulse_sem(a,b) = nan ;

                SRpulseHyp_mean(a,b) = nan ;
                SRpulseHyp_std(a,b) = nan ; 
                SRpulseHyp_sem(a,b) = nan ;
            end       
        end
    end

    % delta spike rate (pulse minus background)
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMinBg{a}{b} = SRpulse{a}{b} - SRbg{a}{b} ;
                SRpulseMinBg_mean(a,b) = mean(SRpulseMinBg{a}{b}) ;
                SRpulseMinBg_std(a,b) = std(SRpulseMinBg{a}{b}) ;
                SRpulseMinBg_sem(a,b) = SRpulseMinBg_std(a,b)/sqrt(NumTrials(a,b)) ;

                SRpulseHypMinBg{a}{b} = SRpulseHyp{a}{b} - SRbg{a}{b} ;
                SRpulseHypMinBg_mean(a,b) = mean(SRpulseHypMinBg{a}{b}) ;
                SRpulseHypMinBg_std(a,b) = std(SRpulseHypMinBg{a}{b}) ; 
                SRpulseHypMinBg_sem(a,b) = SRpulseHypMinBg_std(a,b)/sqrt(NumTrials(a,b)) ; 
            else
                SRpulseMinBg_mean(a,b) = nan ;
                SRpulseMinBg_std(a,b) = nan ;
                SRpulseMinBg_sem(a,b) = nan ;

                SRpulseHypMinBg_mean(a,b) = nan ;
                SRpulseHypMinBg_std(a,b) = nan ;
                SRpulseHypMinBg_sem(a,b) = nan ;
            end
        end
    end

    % delta spike rate (pulse minus background)NORMALIZED
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMinBg_norm{a}{b} = SRpulseMinBg{a}{b}/max(SRpulseMinBg_mean(:)) ;
                SRpulseMinBg_norm_mean(a,b) = mean(SRpulseMinBg_norm{a}{b}) ;
                SRpulseMinBg_norm_std(a,b) = std(SRpulseMinBg_norm{a}{b}) ;
                SRpulseMinBg_norm_sem(a,b) = SRpulseMinBg_norm_std(a,b)/sqrt(NumTrials(a,b)) ;

                SRpulseHypMinBg_norm{a}{b} = SRpulseHypMinBg{a}{b}/min(SRpulseHypMinBg_mean(:)) ;
                SRpulseHypMinBg_norm_mean(a,b) = mean(SRpulseHypMinBg_norm{a}{b}) ;
                SRpulseHypMinBg_norm_std(a,b) = std(SRpulseHypMinBg_norm{a}{b}) ;
                SRpulseHypMinBg_norm_sem(a,b) = SRpulseHypMinBg_norm_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRpulseMinBg_norm_mean(a,b) = nan ;
                SRpulseMinBg_norm_std(a,b) = nan ;
                SRpulseMinBg_norm_sem(a,b) = nan ;

                SRpulseHypMinBg_norm_mean(a,b) = nan ;
                SRpulseHypMinBg_norm_std(a,b) = nan ;
                SRpulseHypMinBg_norm_sem(a,b) = nan ;
            end
        end
    end

    % delta spike rate (pulse minus rest)
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMinRest{a}{b} = SRpulse{a}{b} - SRrest{a}{b} ;
                SRpulseMinRest_mean(a,b) = mean(SRpulseMinRest{a}{b}) ;
                SRpulseMinRest_std(a,b) = std(SRpulseMinRest{a}{b}) ;
                SRpulseMinRest_sem(a,b) = SRpulseMinRest_std(a,b)/sqrt(NumTrials(a,b)) ;

                SRpulseHypMinRest{a}{b} = SRpulseHyp{a}{b} - SRrest{a}{b} ;
                SRpulseHypMinRest_mean(a,b) = mean(SRpulseHypMinRest{a}{b}) ;
                SRpulseHypMinRest_std(a,b) = std(SRpulseHypMinRest{a}{b}) ;
                SRpulseHypMinRest_sem(a,b) = SRpulseHypMinRest_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRpulseMinRest_mean(a,b) = nan ;
                SRpulseMinRest_std(a,b) = nan ;
                SRpulseMinRest_sem(a,b) = nan ;

                SRpulseHypMinRest_mean(a,b) = nan ;
                SRpulseHypMinRest_std(a,b) = nan ;
                SRpulseHypMinRest_sem(a,b) = nan ;
            end
        end
    end

    % delta spike rate (pulse minus rest)NORMALIZED
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMinRest_norm{a}{b} = SRpulseMinRest{a}{b}/max(SRpulseMinRest_mean(:)) ;
                SRpulseMinRest_norm_mean(a,b) = mean(SRpulseMinRest_norm{a}{b}) ;
                SRpulseMinRest_norm_std(a,b) = std(SRpulseMinRest_norm{a}{b}) ;
                SRpulseMinRest_norm_sem(a,b) = SRpulseMinRest_norm_std(a,b)/sqrt(NumTrials(a,b)) ;

                SRpulseHypMinRest_norm{a}{b} = SRpulseHypMinRest{a}{b}/min(SRpulseHypMinRest_mean(:)) ;
                SRpulseHypMinRest_norm_mean(a,b) = mean(SRpulseHypMinRest_norm{a}{b}) ;
                SRpulseHypMinRest_norm_std(a,b) = std(SRpulseHypMinRest_norm{a}{b}) ;
                SRpulseHypMinRest_norm_sem(a,b) = SRpulseHypMinRest_norm_std(a,b)/sqrt(NumTrials(a,b)) ;
            else
                SRpulseMinRest_norm_mean(a,b) = nan ;
                SRpulseMinRest_norm_std(a,b) = nan ;
                SRpulseMinRest_norm_sem(a,b) = nan ;

                SRpulseHypMinRest_norm_mean(a,b) = nan ;
                SRpulseHypMinRest_norm_std(a,b) = nan ;
                SRpulseHypMinRest_norm_sem(a,b) = nan ;
            end
        end
    end

    % background activity subtracted voltage and psth vectors
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    PsthMinBg{a}{b}(c,:) = Psth{a}{b}(c,:) - SRbg{a}{b}(c) ;
                end
            PsthMinBg_mean{a}{b} = mean(PsthMinBg{a}{b},1) ;

            PsthMinBg_mean_norm{a}{b} = PsthMinBg_mean{a}{b}/max(PsthMinBg_mean{a}{b}(iopb:ibpe)) ;
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
          
    % spike raster
    figure
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            for c = 1:NumTrials(a,b) ;
                for d=1:length(spikePnt{a}{b}{c}) ;
                    plot([1,1]*spikePnt{a}{b}{c}(d),[odorRspTrials{a}{b}(c)-1,odorRspTrials{a}{b}(c)],'Color',colorMat{b}(a,:))
                    hold on
                end
            end
        end
    end
%     
% %    mean voltage and psth for different concentrations and backgrounds and time
%     figure 
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(NumConcentrations+2,2,b*2-1:b*2)
%                 plot(time,Psth_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%                 
%                 ylabel(num2str(Concentrations(b)))
%                 hold on
%                 
%             end
%         end
%     end
%     
%     subplot(NumConcentrations+2,2,2*NumConcentrations+1:2*NumConcentrations+2)
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 plot(tDataN{a}{b},SRpulse{a}{b},'o','Color',colorMat{b}(a,:))
%                 hold on
%             end
%         end
%     end
%     xlabel('trig time (sec)')
%     ylabel('spike rate (hz)')
%     
%     subplot(NumConcentrations+2,2,2*NumConcentrations+3:2*NumConcentrations+4)
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 plot(tDataN{a}{b},SRrest{a}{b},'*','Color',colorMat{b}(a,:))
%                 hold on
%                 plot(tDataN{a}{b},SRbg{a}{b},'.','Color',colorMat{b}(a,:))
%             end
%         end
%     end
%     xlabel('trig time (sec)')
%     ylabel('spike rate (hz)')
%     
%     axes('Position',[0 0 .02 1],'Visible','off');
%     text(0,.01,RepStat,'FontSize',5)
%     
%     %mean voltage and psth with background subtracted for different concentrations and backgrounds and time
%     figure 
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(NumConcentrations+2,2,b*2-1:b*2)
%                 plot(time,PsthMinBg_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%                 
%                 ylabel(num2str(Concentrations(b)))
%                 hold on
%                 %plot(time(1,[iorb(a,b),iore(a,b)]),Psth_mean{a}{b}(1,[iorb(a,b),iore(a,b)]),'o','Color',colorMat{NumConcentrations}(a,:))
%                 
%             end
%         end
%     end
%     
%     subplot(NumConcentrations+2,2,2*NumConcentrations+1:2*NumConcentrations+4)
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 plot(tDataN{a}{b},SRrest{a}{b},'*','Color',colorMat{b}(a,:))
%                 hold on
%                 plot(tDataN{a}{b},SRpulseMinBg{a}{b},'+','Color',colorMat{b}(a,:))
%             end
%         end
%     end
%     xlabel('trig time (sec)')
%     ylabel('spike rate (hz)')
%     
%     % comparing kinetics of pulse responses
%     figure 
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;  
%                 subplot(1,NumConcentrations,b)
%                 plot(time,PsthMinBg_mean_norm{a}{b},'Color',colorMat{b}(a,:))
%                 hold on
%                 xlim([time(iopb),time(ibpe)])
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
%                 plot(time,Psth_mean{a}{b}-SRrest_mean(a,b),'Color',colorMat{b}(a,:))            
%                 hold on                 
%                 
%             end
%         end
%     end
%     
%     % comparing background responses with pulse response
%     figure
%     for a = 1:NumBackgrounds ; % for each background
%         for b = 1:NumConcentrations ; % for concentration
%             if NumTrials(a,b)>0 ;
%                 subplot(3,1,1)
%                 plot(SRbg{a}{b}-SRrest{a}{b},SRpulse{a}{b}-SRbg{a}{b},'+','Color',colorMat{b}(a,:))
%                 xlabel('SRbg - SRrest (mV)')
%                 ylabel('SRpulse - SRbg (hz)')
%                 hold on        
%                 
%                 subplot(3,1,2)
%                 plot(SRtbg{a}{b}-SRbg{a}{b},SRpulse{a}{b}-SRbg{a}{b},'+','Color',colorMat{b}(a,:)) 
%                 xlabel('SRtbg - SRbg (mV)')
%                 ylabel('SRpulse - SRbg (hz)')
%                 hold on            
%                 
%                 subplot(3,1,3)
%                 plot(SRtbg{a}{b}-SRrest{a}{b},SRpulse{a}{b}-SRbg{a}{b},'+','Color',colorMat{b}(a,:)) 
%                 xlabel('SRtbg - SRrest (mV)')
%                 ylabel('SRpulse - SRbg (hz)')
%                 hold on              
%             end
%         end
%     end
%     
% 
%     %  spikes as a function of concentration
%     figure
%     subplot(3,2,1)
%     for a = 1:NumBackgrounds ;
%         errorbar(log10(Concentrations),SRrest_mean(a,:),SRrest_std(a,:),SRrest_std(a,:),'*:','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%         errorbar(log10(Concentrations),SRbg_mean(a,:),SRbg_std(a,:),SRbg_std(a,:),'*--','Color',colorMat{NumConcentrations}(a,:))
%         errorbar(log10(Concentrations),SRpulse_mean(a,:),SRpulse_std(a,:),SRpulse_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
%     end
%     ylabel('spike rate')
%     xlabel('log concentration')
% 
%     subplot(3,2,3)
%     for a = 1:NumBackgrounds ;
%         errorbar(log10(Concentrations),SRpulseMinBg_norm_mean(a,:),SRpulseMinBg_norm_std(a,:),SRpulseMinBg_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%     end
%     ylabel('pulse min bg spike rate')
%     xlabel('log concentration')
% 
%     subplot(3,2,5)
%     for a = 1:NumBackgrounds ;
%         errorbar(log10(Concentrations),SRpulseMinRest_norm_mean(a,:),SRpulseMinRest_norm_std(a,:),SRpulseMinRest_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%     end
%     ylabel('pulse min rest spike rate')
%     xlabel('log concentration')
% 
% 
%     % hyperpol
%     subplot(3,2,2)
%     for a = 1:NumBackgrounds ;
%         errorbar(log10(Concentrations),SRrest_mean(a,:),SRrest_std(a,:),SRrest_std(a,:),'*:','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%         errorbar(log10(Concentrations),SRbg_mean(a,:),SRbg_std(a,:),SRbg_std(a,:),'*--','Color',colorMat{NumConcentrations}(a,:))
%         errorbar(log10(Concentrations),SRpulseHyp_mean(a,:),SRpulseHyp_std(a,:),SRpulseHyp_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
%     end
%     ylabel('spike rate')
%     xlabel('log concentration')
% 
%     subplot(3,2,4)
%     for a = 1:NumBackgrounds ;
%         errorbar(log10(Concentrations),SRpulseHypMinBg_norm_mean(a,:),SRpulseHypMinBg_norm_std(a,:),SRpulseHypMinBg_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%     end
%     ylabel('pulse min bg spike rate')
%     xlabel('log concentration')
% 
%     subplot(3,2,6)
%     for a = 1:NumBackgrounds ;
%         errorbar(log10(Concentrations),SRpulseHypMinRest_norm_mean(a,:),SRpulseHypMinRest_norm_std(a,:),SRpulseHypMinRest_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
%         hold on
%     end
%     ylabel('pulse min rest spike rate')
%     xlabel('log concentration')
% 
% graphs for single page focus
    figure
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                subplot(NumConcentrations+2,2,b*2-1)
                plot(time,Psth_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
                axis tight
                ylabel(num2str(Concentrations(b)))
                hold on 
            end
        end
    end
    
    subplot(NumConcentrations+2,2,6:2:NumConcentrations*2)
    for a = 1:NumBackgrounds ;
        errorbar(log10(Concentrations),SRpulseMinBg_norm_mean(a,:),SRpulseMinBg_norm_std(a,:),SRpulseMinBg_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
        hold on
    end
    ylabel('pulse min bg rate')
    xlabel('log concentration')
    axis tight

    subplot(NumConcentrations+2,2,2:2:4)
    plot(time(1,iopb:SRiore(1,3)),vData{1}{3}(1,iopb:SRiore(1,3)))
    axis tight
    
    subplot(NumConcentrations+2,2,2*NumConcentrations+1:2*NumConcentrations+2)
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                plot(tDataN{a}{b},SRpulse{a}{b},'o','Color',colorMat{b}(a,:))
                hold on
            end
        end
    end
    xlabel('trig time (sec)')
    ylabel('spike rate (hz)')
    axis tight
    
    subplot(NumConcentrations+2,2,2*NumConcentrations+3:2*NumConcentrations+4)
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                plot(tDataN{a}{b},SRbg{a}{b},'.','Color',colorMat{b}(a,:))
                hold on
                plot(tDataN{a}{b},SRrest{a}{b},'*','Color',colorMat{b}(a,:))
            end
        end
    end
    xlabel('trig time (sec)')
    ylabel('spike rate (hz)')
    axis tight    



    
    % for igor

    identifier = ['LogConcentration','cell',num2str(A(DB))] ;
    ForIgor.(identifier) = log10(Concentrations) ;
    
    for a = 1:NumBackgrounds ; % for each background, 
        identifier = ['SRminBgMean','Bg',BgConcentrations{a},'ORNcell',num2str(A(DB))] ;
        ForIgor.(identifier) = SRpulseMinBg_norm_mean(a,:) ;
        
        identifier = ['SRminBgSem','Bg',BgConcentrations{a},'ORNcell',num2str(A(DB))] ;
        ForIgor.(identifier) = SRpulseMinBg_norm_sem(a,:) ;  
    end
    
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration 
            if NumTrials(a,b)>0 ;
                identifier = ['Psth','LogPulse',num2str(abs(log10(Concentrations(b)))),'Bg',BgConcentrations{a},'ORNcell',num2str(A(DB))] ;
                ForIgor.(identifier) = Psth_mean{a}{b} ;

                identifier = ['PsthMinBg','LogPulse',num2str(abs(log10(Concentrations(b)))),'Bg',BgConcentrations{a},'ORNcell',num2str(A(DB))] ;
                ForIgor.(identifier) = PsthMinBg_mean{a}{b} ;          
            end
        end
    end

    identifier = ['time','cell',num2str(A(DB))] ;
    ForIgor.(identifier) = time; 

    
    % temp results for population analysis
    
    % range of background and pulse concentrations assessed across
    % population and indicies to compartmentalize this DB correctly
    PopData.ConcentrationRange = 10.^[-8:-2] ;    
    PopData.BgConcentrationRange = {'0','Log8','Log7','Log6','Log5','Log4','Log3','wash'} ;
    
    [c,PopPulsei] = intersect(PopData.ConcentrationRange,Concentrations) ;
    
    for a = 1:NumBackgrounds ;
        for b = 1:length(PopData.BgConcentrationRange) ;
            if strcmp(BgConcentrations{a},PopData.BgConcentrationRange{b}) ;
                PopBgi(a) = b ;
            end
        end
    end
    
    % {bg}{pulse}(cell,time)  
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            for b=1:length(PopData.ConcentrationRange) ;
                PopData.Psth_mean{a}{b} = nan(length(A),length(time)) ;
                PopData.PsthMinBg_mean{a}{b} = nan(length(A),length(time)) ;
            end
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                PopData.Psth_mean{PopBgi(a)}{PopPulsei(b)}(DB,:) = Psth_mean{a}{b} ;
                PopData.PsthMinBg_mean{PopBgi(a)}{PopPulsei(b)}(DB,:) = PsthMinBg_mean{a}{b} ;
            end
        end
    end 
    
    % {bg}(cell, pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRpulseMinBg_norm_mean{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMinBg_norm_sem{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        PopData.SRpulseMinBg_norm_mean{PopBgi(a)}(DB,PopPulsei) = SRpulseMinBg_norm_mean(a,:) ;
        PopData.SRpulseMinBg_norm_sem{PopBgi(a)}(DB,PopPulsei) = SRpulseMinBg_norm_sem(a,:) ;
    end
    
        % {bg}(cell, pulse)
    if DB==1 ;
        for a=1:length(PopData.BgConcentrationRange) ;
            PopData.SRpulseMinBg_mean{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
            PopData.SRpulseMinBg_sem{a} = nan(length(A),length(PopData.ConcentrationRange)) ;
        end
    end
    for a = 1:NumBackgrounds ; % for each background
        PopData.SRpulseMinBg_mean{PopBgi(a)}(DB,PopPulsei) = SRpulseMinBg_mean(a,:) ;
        PopData.SRpulseMinBg_sem{PopBgi(a)}(DB,PopPulsei) = SRpulseMinBg_sem(a,:) ;
    end
    
    % {bg}{pulse}(cell,[contPulse,bgPulse])
    if DB==1 ;
        for a=2:length(PopData.BgConcentrationRange) ;
            for b=1:length(PopData.ConcentrationRange) ; 
                PopData.SRpulseMinBg_mean_paired{a}{b} = nan(2,length(A)) ;
                PopData.SRpulseMinBg_sem_paired{a}{b} = nan(2,length(A))  ; 
            end
        end
    end
    for a = 2:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration 
            if NumTrials(a,b)>0 ;
                PopData.SRpulseMinBg_mean_paired{PopBgi(a)}{PopPulsei(b)}(:,DB) = [SRpulseMinBg_mean(1,b);SRpulseMinBg_mean(a,b)] ;
                PopData.SRpulseMinBg_sem_paired{PopBgi(a)}{PopPulsei(b)}(:,DB) = [SRpulseMinBg_sem(1,b);SRpulseMinBg_sem(a,b)] ; 
            end    
        end
    end
    
    % clear all variables except ForIgor and PopData structures
    clearvars -except Input A ForIgor PopData  
end % Data block loop


% population calculations
for a = 1:length(PopData.SRpulseMinBg_norm_mean) ; % for each background assessed in the population 
    PopData.SRpulseMinBg_norm_mean_PopMean{a}= nanmean(PopData.SRpulseMinBg_norm_mean{a},1) ;
    PopData.SRpulseMinBg_norm_mean_PopSem{a} = nanstd(PopData.SRpulseMinBg_norm_mean{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMinBg_norm_mean{a}),1)) ;
    
    PopData.SRpulseMinBg_mean_PopMean{a}= nanmean(PopData.SRpulseMinBg_mean{a},1) ;
    PopData.SRpulseMinBg_mean_PopSem{a} = nanstd(PopData.SRpulseMinBg_mean{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMinBg_mean{a}),1)) ;
end



% population figures
Conc(1,:) = [.6,.5,.4,0,0,0,0] ;
Conc(2,:) = [1,1,1,1,.6,.5,.4] ;
Conc(3,:) = [.6,.5,.4,.3,.2,.1,0] ;
for a=1:length(PopData.ConcentrationRange) ; % each concentration is a matrix and each background is a row within that matrix
    colorMat{a} = [Conc(3,a),Conc(3,a),Conc(3,a); Conc(2,a),Conc(1,a),Conc(1,a); Conc(1,a),Conc(1,a),Conc(2,a);... 
        Conc(1,a),Conc(2,a),Conc(1,a);Conc(2,a),Conc(2,a),Conc(1,a);...
        Conc(2,a),Conc(1,a),Conc(2,a);Conc(1,a),Conc(2,a),Conc(2,a);Conc(2,a),Conc(2,a),Conc(2,a)] ;
end


figure
for a = 1:length(PopData.SRpulseMinBg_norm_mean) ; % for each background assessed in the population
    for b = 1:length(A) ; % for each data block in population
        subplot(2,2,1)
        errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMinBg_mean{a}(b,:),PopData.SRpulseMinBg_sem{a}(b,:),PopData.SRpulseMinBg_sem{a}(b,:),'+-','Color',colorMat{5}(a,:))
        hold on
        
        subplot(2,2,2)
        errorbar(log10(PopData.ConcentrationRange),PopData.SRpulseMinBg_norm_mean{a}(b,:),PopData.SRpulseMinBg_norm_sem{a}(b,:),PopData.SRpulseMinBg_norm_sem{a}(b,:),'+-','Color',colorMat{5}(a,:))
        hold on 
    end
        
    subplot(2,2,3)
    errorbar(log10(PopData.ConcentrationRange), PopData.SRpulseMinBg_mean_PopMean{a}, PopData.SRpulseMinBg_mean_PopSem{a}, PopData.SRpulseMinBg_mean_PopSem{a},'+-','Color',colorMat{5}(a,:))
    hold on

    subplot(2,2,4)
    errorbar(log10(PopData.ConcentrationRange), PopData.SRpulseMinBg_norm_mean_PopMean{a}, PopData.SRpulseMinBg_norm_mean_PopSem{a}, PopData.SRpulseMinBg_norm_mean_PopSem{a},'+-','Color',colorMat{5}(a,:))
    hold on         
end


figure
for a = 1:length(PopData.ConcentrationRange) ; % for each pulse concentration
    subplot(1,2,1)
    plot(nanmean(PopData.PsthMinBg_mean{1}{a},1),'color',colorMat{a}(1,:))
    hold on

    subplot(1,2,2)
    plot(nanmean(PopData.PsthMinBg_mean{1}{a},1)/max(nanmean(PopData.PsthMinBg_mean{1}{a},1)),'color',colorMat{a}(1,:))
    hold on
end



identifier = ['LogConcentration','AllMean'] ;
ForIgor.(identifier) = log10(PopData.ConcentrationRange) ;

for a = 1:length(PopData.SRpulseMinBg_norm_mean) ; % for each background assessed in the population 
    identifier = ['SRminBgMean','Bg',PopData.BgConcentrationRange{a},'ORNAllMean'] ;
    ForIgor.(identifier) = nanmean(PopData.SRpulseMinBg_norm_mean{a},1) ;

    identifier = ['SRminBgSem','Bg',PopData.BgConcentrationRange{a},'ORNAllMean'] ;
    ForIgor.(identifier) = nanmean(PopData.SRpulseMinBg_norm_sem{a},1) ;  
end

for a = 1:length(PopData.SRpulseMinBg_norm_mean) ; % for each background assessed in the population 
    identifier = ['SRminBgMean','Bg',PopData.BgConcentrationRange{a},'ORNAllSem'] ;
    ForIgor.(identifier) = nanstd(PopData.SRpulseMinBg_norm_mean{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMinBg_norm_mean{a}),1)) ;

    identifier = ['SRminBgSem','Bg',PopData.BgConcentrationRange{a},'ORNAllSem'] ;
    ForIgor.(identifier) = nanstd(PopData.SRpulseMinBg_norm_sem{a},[],1)./sqrt(sum(~isnan(PopData.SRpulseMinBg_norm_sem{a}),1)) ; 
end

for a = 1:length(PopData.Psth_mean) ; % for each background
    for b = 1:length(PopData.Psth_mean{a}) ; % for concentration 

            identifier = ['Psth','LogPulse',num2str(abs(log10(PopData.ConcentrationRange(b)))),'Bg',PopData.BgConcentrationRange{a},'ORNAllMean'] ;
            ForIgor.(identifier) = nanmean(PopData.Psth_mean{a}{b},1) ;

            identifier = ['PsthMinBg','LogPulse',num2str(abs(log10(PopData.ConcentrationRange(b)))),'Bg',PopData.BgConcentrationRange{a},'ORNAllMean'] ;
            ForIgor.(identifier) = nanmean(PopData.PsthMinBg_mean{a}{b},1) ;        
    end
end

for a = 1:length(PopData.Psth_mean) ; % for each background
    for b = 1:length(PopData.Psth_mean{a}) ; % for concentration 

            identifier = ['Psth','LogPulse',num2str(abs(log10(PopData.ConcentrationRange(b)))),'Bg',PopData.BgConcentrationRange{a},'ORNAllSem'] ;
            ForIgor.(identifier) = nanstd(PopData.Psth_mean{a}{b},[],1)./sqrt(sum(~isnan(PopData.Psth_mean{a}{b}),1)) ;

            identifier = ['PsthMinBg','LogPulse',num2str(abs(log10(PopData.ConcentrationRange(b)))),'Bg',PopData.BgConcentrationRange{a},'ORNAllSem'] ;
            ForIgor.(identifier) = nanstd(PopData.PsthMinBg_mean{a}{b},[],1)./sqrt(sum(~isnan(PopData.PsthMinBg_mean{a}{b}),1)) ;        
    end
end


% {bg}{pulse}(cell,[contPulse,bgPulse])
for a = 1:length(PopData.Psth_mean) ; % for each background
    for b = 1:length(PopData.Psth_mean{a}) ; % for concentration 
        if ~isempty(PopData.SRpulseMinBg_mean_paired{a}) ; % if that was assessed in the population  
            if ~isempty(PopData.SRpulseMinBg_mean_paired{a}{b}) ; % if that was assessed in the population  
                identifier = ['SRminBgMean','LogPulse',num2str(abs(log10(PopData.ConcentrationRange(b)))),'Bg',PopData.BgConcentrationRange{a},'ORNAllCnt'] ;
                ForIgor.(identifier) = PopData.SRpulseMinBg_mean_paired{a}{b}(1,:) ;

                identifier = ['SRminBgSem','LogPulse',num2str(abs(log10(PopData.ConcentrationRange(b)))),'Bg',PopData.BgConcentrationRange{a},'ORNAllCnt'] ;
                ForIgor.(identifier) = PopData.SRpulseMinBg_sem_paired{a}{b}(1,:) ;


                identifier = ['SRminBgMean','LogPulse',num2str(abs(log10(PopData.ConcentrationRange(b)))),'Bg',PopData.BgConcentrationRange{a},'ORNAllBg'] ;
                ForIgor.(identifier) = PopData.SRpulseMinBg_mean_paired{a}{b}(2,:) ;

                identifier = ['SRminBgSem','LogPulse',num2str(abs(log10(PopData.ConcentrationRange(b)))),'Bg',PopData.BgConcentrationRange{a},'ORNAllBg'] ;
                ForIgor.(identifier) = PopData.SRpulseMinBg_sem_paired{a}{b}(2,:) ;

            end
        end   
    end
end
    
   


