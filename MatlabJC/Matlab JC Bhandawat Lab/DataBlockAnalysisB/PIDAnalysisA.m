function ForIgor = PIDAnalysisA(Input,A,AnalysisFlag)

% to analyze PID measures of pulse +/- odor db 
% JC 3/31/14


for DB = 1:length(A) ; % for each cell

    %% parameters

    PopData.ConcentrationRange = [10^-9,10.^[-7:-2]] ; % 10^-9 is parafin oil
    PopData.BgConcentrationRange = {'0','Log7','Log6','Log5','Log4'} ;
    PopData.BgConcentrationRangeNum = [10^-9,10.^[-7:-4]] ; % numerical version, control changed to 10^-x and wash not included
    
    PopData.sampleRate = 10000 ; % (hz) temp hard coded - should be saved in file 
    bgTransTime = 5 ; %(sec) time after background odor starts to avoid transient response
    OdorRspTime = .05 ; % (sec) time before and after peak of mean psth which odor response is assessed
    PopData.RspTime = 1.5 ; % (sec) time after odor pulse onset over which odor respone ensemble is defined

    id1 = 'OdorRsp' ;
    id2 = 'OdorConcentration' ;
    id3 = 'BgConcentration' ;
    
    PopData.PidDataPath = ['Z:/Cafaro Documents/Analysis/PidData/'] ;

    %% load data in matricies
    rootdir = ['Z:\Cafaro Data Backup\', Input(A(DB)).cellname(1:6),'Data'];

    Concentrations = str2num(Input(A(DB)).(id2)) ;
    NumConcentrations = length(Concentrations) ;

    BgConcentrations = Input(A(DB)).(id3) ;
    NumBackgrounds = length(BgConcentrations) ;

    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration

            if length(Input(A(DB)).(id1){a}{b})>2 ;
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
            else
                NumTrials(a,b) = 0 ;
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

    %% index of background odor and odor pulse
    iopb = find(ao0DataB{2}{2}(1,:)~=0,1,'first')-1 ; % odor pulse begining
    iope = find(ao0DataB{2}{2}(1,:)~=0,1,'last')-1 ; % odor pulse ending

    ibpb = find(ao1DataB{2}{2}(1,:)~=0,1,'first')-1 ; % background pulse beginging
    ibpe = find(ao1DataB{2}{2}(1,:)~=0,1,'last')-1 ; % background pulse end
           
    %% mean voltage data
    for a = 1:NumBackgrounds ;
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ; 
                vData_mean{a}{b} = mean(vData{a}{b}) ;
                vData_sem{a}{b} = std(vData{a}{b})/sqrt(NumTrials(a,b)) ;
            end
        end
    end
       
    %% odor response time
    OdorRspPnts = OdorRspTime*PopData.sampleRate ;
    
    for a = 1:NumBackgrounds ;
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;                 
                [m,mi] = max(vData_mean{a}{b}(iopb:iopb+PopData.sampleRate)) ; % max point of mean psth within a second of odor pulse onset
               
                SRiorb(a,b) = mi-1+iopb - OdorRspPnts ; % point of odor pulse begining
                SRiore(a,b) = mi-1+iopb +  OdorRspPnts ; % point of odor pulse end
            end
        end
    end
       
    %% assess rest spike rate
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration   
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRrestMean{a}{b}(c) = mean(vData{a}{b}(c,1:ibpb)) ; % first point: background pulse start
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
                    SRbgMean{a}{b}(c) = mean(vData{a}{b}(c,ibpb+bgTransPnts:iopb)) ; % (start of background + transient time: begining of odor pulse)
                    SRbgStd{a}{b}(c) = std(vData{a}{b}(c,ibpb+bgTransPnts:iopb)) ; % variance during background
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
        
    %% odor pulse response spike rate
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    SRpulseMean{a}{b}(c) = mean(vData{a}{b}(c,SRiorb(a,b):SRiore(a,b))); %mV (odor response depol start: odor response depol end)
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
     
    %% pulse minus parfin (assuming parafin is listed as zero background lowest pulse strength in data block)
     for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                SRpulseMeanMinBgMinP_mean(a,b) = SRpulseMeanMinBg_mean(a,b) - SRpulseMeanMinBg_mean(1,1) ;
                SRpulseMeanMinBgMinP_sem(a,b) = sqrt(SRpulseMeanMinBg_std(a,b)^2+SRpulseMeanMinBg_std(1,1)^2) ;  
            else
                SRpulseMeanMinBgMinP_mean(a,b) = nan ;
                SRpulseMeanMinBgMinP_sem(a,b) = nan ;
            end
        end
     end
     
    %% fit pulse minus parafin with straight line begining at origin
     [MinPfit_R,MinPfit_m,MinPfit_b] = regression(log10(Concentrations(2:end)), log10(SRpulseMeanMinBgMinP_mean(1,2:end))) ;
    
    %% background steady state / pulse peak of same strength
    for a = 1:NumBackgrounds ; % for each background
         for b = 1:NumConcentrations ; % for concentration
             if NumTrials(a,b) && Concentrations(b)==PopData.BgConcentrationRangeNum(a) ; % CAUTION USING POPDATA vector which assumes all DB backgrounds are the same
                 SRbgDivPulse_mean(a) = SRbgMeanMinRest_mean(a,b)/SRpulseMeanMinBg_mean(1,b) ;
                 SRbgDivPulse_sem(a) = SRbgDivPulse_mean(a)*sqrt((SRbgMeanMinRest_sem(a,b)/SRbgMeanMinRest_mean(a,b))^2 + (SRpulseMeanMinBg_sem(1,b)/SRpulseMeanMinBg_mean(1,b))^2) ;
             end
         end
    end
    
    %% shifted pulse curve
    shiftedSRpulseMeanMinBg_mean = SRpulseMeanMinBg_mean(1,:)*mean(SRbgDivPulse_mean(2:end)) ; % dose-response * average shift factor
    
    %% rest and bg activity subtracted psth vectors
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    vDataMinRest{a}{b}(c,:) = vData{a}{b}(c,:) - SRrestMean{a}{b}(c) ;
                    vDataMinBg{a}{b}(c,:) = vData{a}{b}(c,:) - SRbgMean{a}{b}(c) ;
                end
            vDataMinRest_mean{a}{b} = mean(vDataMinRest{a}{b},1) ;
            vDataMinBg_mean{a}{b} = mean(vDataMinBg{a}{b},1) ;
            
            vDataMinRest_std{a}{b} = std(vDataMinRest{a}{b},[],1) ;
            vDataMinBg_std{a}{b} = std(vDataMinBg{a}{b},[],1) ;
            
            vDataMinRest_sem{a}{b} = vDataMinRest_std{a}{b}/sqrt(NumTrials(a,b)) ;
            vDataMinBg_sem{a}{b} = vDataMinBg_std{a}{b}/sqrt(NumTrials(a,b)) ;
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
       
    % mean voltage data for all data  
    figure
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                plot(time,vDataMinRest_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
                hold on              
                plot(time(SRiorb(a,b):SRiore(a,b)),vDataMinRest_mean{a}{b}(SRiorb(a,b):SRiore(a,b)),'y')
            end
        end
    end
    axis tight
    xlabel('time (sec)')
    ylabel('PID-offset (V)')
    
    % mean voltage data for each pulse concentration
    figure 
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                subplot(NumConcentrations,1,b)
                plot(time,vDataMinRest_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
                hold on              
            end
        end
    end
    axis tight
    xlabel('time (sec)')
    ylabel('PID-offset (V)')
    
    % mean voltage data for each pulse concentration
    figure 
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                subplot(NumConcentrations,1,b)
                plot(time,vDataMinBg_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
                hold on              
            end
        end
    end
    axis tight
    xlabel('time (sec')
    ylabel('PID-bg (V)')
    
    % stats as function of time
    figure
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                subplot(2,1,1)
                plot(tDataN{a}{b},SRbgMean{a}{b},'.','Color',colorMat{b}(a,:))
                hold on
                plot(tDataN{a}{b},SRpulseMean{a}{b},'o','Color',colorMat{b}(a,:))            
                axis tight
                set(gca,'yscale','log')
                
                subplot(2,1,2)
                plot(tDataN{a}{b},SRbgMean{a}{b},'.','Color',colorMat{b}(a,:))
                hold on
                plot(tDataN{a}{b},SRrestMean{a}{b},'*','Color',colorMat{b}(a,:))
                axis tight
                set(gca,'yscale','log')
            end
        end
    end
    xlabel('trig time (sec)')
    ylabel('PID offset,bg,pulse (V)')
    
    % stats as function of time
    figure
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                subplot(2,1,1)
                plot(tDataN{a}{b},SRpulseMeanMinBg{a}{b},'o','Color',colorMat{b}(a,:))   
                hold on
                axis tight
                set(gca,'yscale','log')
                
                subplot(2,1,2)
                plot(tDataN{a}{b},SRbgMeanMinRest{a}{b},'.','Color',colorMat{b}(a,:))
                hold on
                axis tight
                set(gca,'yscale','log')
            end
        end
    end
    xlabel('trig time (sec)')
    ylabel('PID pulse-bg, bg-offset(V)')
    
    % pulse peak dose-response curve
    figure
    for a = 1:NumBackgrounds ;
        errorbar(log10(Concentrations), SRpulseMeanMinBg_mean(a,:),SRpulseMeanMinBg_sem(a,:),SRpulseMeanMinBg_sem(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
        hold on
    end
    ylabel('PID pulse-bg')
    xlabel('pulse concentration')
    set(gca,'yscale','log')
    axis tight
    
    % pulse peak dose-response minus parafin trials
    figure
    errorbar(log10(Concentrations), SRpulseMeanMinBg_mean(1,:),SRpulseMeanMinBg_sem(1,:),SRpulseMeanMinBg_sem(1,:),'*-','Color',colorMat{NumConcentrations}(1,:))
    hold on
    errorbar(log10(Concentrations), SRpulseMeanMinBgMinP_mean(1,:),SRpulseMeanMinBgMinP_sem(1,:),SRpulseMeanMinBgMinP_sem(1,:),'o:','Color',colorMat{NumConcentrations}(1,:))
    ylabel('PID pulse-bg')
    xlabel('pulse concentration')
    set(gca,'yscale','log')
    axis tight
    
    % bg dose-response curve
    figure
    for a = 1:NumConcentrations ; % for concentration
        errorbar(log10(PopData.BgConcentrationRangeNum), SRbgMeanMinRest_mean(:,a),SRbgMeanMinRest_sem(:,a),SRbgMeanMinRest_sem(:,a),'*','Color',colorMat{a}(1,:))
        hold on
    end
    set(gca,'yscale','log')
    ylabel('PID bg-offset')
    xlabel('bg concentration')
    
    % bg and pulse dr curves together
    figure
    for a = 1:NumBackgrounds ;
        errorbar(log10(Concentrations), SRpulseMeanMinBg_mean(a,:),SRpulseMeanMinBg_sem(a,:),SRpulseMeanMinBg_sem(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
        hold on
    end
    for a = 1:NumConcentrations ; % for concentration
        errorbar(log10(PopData.BgConcentrationRangeNum), SRbgMeanMinRest_mean(:,a),SRbgMeanMinRest_sem(:,a),SRbgMeanMinRest_sem(:,a),'o','Color',colorMat{a}(1,:))
    end
    plot(log10(Concentrations),  shiftedSRpulseMeanMinBg_mean,':','Color',colorMat{NumConcentrations}(1,:))
    set(gca,'yscale','log')
    ylabel('PID pulse-bg,bg-offset')
    xlabel('odor concentration')
 
    % comparing background and pulse rising phases
    figure 
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                if a==1 ; % if pulse - bg
                    plot(vDataMinRest_mean{a}{b}(iopb:end),'Color',colorMat{NumConcentrations}(a,:))
                    hold on   
                else
                    plot(vDataMinRest_mean{a}{b}(ibpb:end),'Color',colorMat{NumConcentrations}(a,:))
                    hold on
                end
            end
        end
    end
    axis tight
    xlabel('time (sec)')
    ylabel('PID-offset (V)')
    
   
%% ForIgor

% pid pulse peak pnt curves
ForIgor.PidPk = SRpulseMeanMinBg_mean(:,2:end) ;
ForIgor.PidPkSem = SRpulseMeanMinBg_sem(:,2:end) ;
ForIgor.PidPkMinPara = SRpulseMeanMinBgMinP_mean(:,2:end) ;
ForIgor.PidPkMinParaSem = SRpulseMeanMinBgMinP_sem(:,2:end) ;
ForIgor.LogPulseConcentration = log10(Concentrations(2:end)) ; 

ForIgor.PidPkPara = SRpulseMeanMinBg_mean(1,1) ;
ForIgor.PidPkParaSem = SRpulseMeanMinBg_sem(1,1) ;

% pid steady state background curves
for a = 1:NumBackgrounds ; % for each background
    if NumTrials(a,a)>0 ;
        ForIgor.PidbgSs(a) = SRbgMeanMinRest_mean(a,a) ; % used the backgrounds when the same pulse was delivered
        ForIgor.PidbgSsSem(a) = SRbgMeanMinRest_sem(a,a) ;
    else
        ForIgor.PidbgSs(a) = nan ;
        ForIgor.PidbgSsSem(a) = nan ;
    end
end
ForIgor.LogBgConcentrations = [-7:1:-4] ;

% shifted peak dr curve
ForIgor.PidPkshift = shiftedSRpulseMeanMinBg_mean(1,2:end) ;

% comparing pulse and background
ForIgor.PidMinRest0shift = vDataMinRest_mean{1}{2}(iopb:end) ;
ForIgor.PidMinRest7shift =  vDataMinRest_mean{2}{2}(ibpb:end) ;

% pid waveforms minus
for a = 1:NumConcentrations ; % for each background
    if NumTrials(1,a)>0 ;
        ForIgor.PidMinBg0(a,:) = vDataMinBg_mean{1}{a} ;
        ForIgor.PidMinRest0(a,:) = vDataMinRest_mean{1}{a} ;
    else
        ForIgor.PidMinBg0(a,:) = nan(1,length(time)) ;
        ForIgor.PidMinRest0(a,:) = nan(1,length(time)) ;
    end
    
    if NumTrials(2,a)>0 ;
        ForIgor.PidMinBg7(a,:) = vDataMinBg_mean{2}{a} ;
        ForIgor.PidMinRest7(a,:) = vDataMinRest_mean{2}{a} ;
    else
        ForIgor.PidMinBg7(a,:) = nan(1,length(time)) ;
        ForIgor.PidMinRest7(a,:) = nan(1,length(time)) ;
    end
    
    if NumTrials(3,a)>0 ;
        ForIgor.PidMinBg6(a,:) = vDataMinBg_mean{3}{a} ;
        ForIgor.PidMinRest6(a,:) = vDataMinRest_mean{3}{a} ;
    else
        ForIgor.PidMinBg6(a,:) = nan(1,length(time)) ;
        ForIgor.PidMinRest6(a,:) = nan(1,length(time)) ;
    end
    
    if NumTrials(4,a)>0 ;
        ForIgor.PidMinBg5(a,:) = vDataMinBg_mean{4}{a} ;
        ForIgor.PidMinRest5(a,:) = vDataMinRest_mean{4}{a} ;
    else
        ForIgor.PidMinBg5(a,:) = nan(1,length(time)) ;
        ForIgor.PidMinRest5(a,:) = nan(1,length(time)) ;
    end
    
    if NumTrials(5,a)>0 ;
        ForIgor.PidMinBg4(a,:) = vDataMinBg_mean{5}{a} ;
        ForIgor.PidMinRest4(a,:) = vDataMinRest_mean{5}{a} ;
    else
        ForIgor.PidMinBg4(a,:) = nan(1,length(time)) ;
        ForIgor.PidMinRest4(a,:) = nan(1,length(time)) ;
    end
end
    
   

end 


