function ForIgor = CaDcMidgetAnalyzer(Input,Parameters,A,id,amp) ;

% this function will compare the spike output of cells stimulated by light
% recorded cell attached and cells stimulated by dynamic clamp recorded
% whole cell (for suplemental figure 5).

% 7/19/10 modified to get std of psth for each cell to compare variability
% of dc vs cell attached.

% JC 4/27/10

A={116,117,118,119,71,72,75,79,80,81,82} ;
id = {'CA','CA','CA','CA','dcSimShuffle','dcSimShuffle','dcSimShuffle','dcSimShuffle','dcSimShuffle','dcSimShuffle','dcSimShuffle'} ;

for a = 1:length(A) ; % for each A
    Smoothtime = 0.01 ; % secs  
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/060809Ec2']); % light stim seed 1
    [stm, error] = ITCReadEpochStm(246, 0, fp);
    stim = stm(100001:200000) ; % no prepoints
    
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A{a}).cellname]);
    
    epochs = str2num(Input(A{a}).(id{a})) ; %
    
    if strcmp(id{a},'CA') ; % if data is cell attached
        skpEpochs = 1 ;
    else 
        skpEpochs = 2 ;
    end
        
    round = 0;
    for b = 1:skpEpochs:length(epochs) ;
        round= round+1 ;
        [Data(round,:), error] = ITCReadEpoch(epochs(b), Input(A{a}).amp, fp) ; 
        [SI(round), error] = ITCGetSamplingInterval(epochs(b), fp); % get sampling interval    
        SI(round) = SI(round) * 1e-6; % Sampling interval in sec
    end
    

    if Input(A{a}).ITC18flag == 1 ;
        SI = SI*1.25 ;
    end
    time{a} = [SI(1):SI(1):SI(1)*length(Data)] ;

    % spike detection
    if strcmp(id{a},'CA') ; % if data is cell attached
        Data = highPassFilter(Data, 1/SI(1), 20) ; % remove drift below 20Hz
        DataShift = circshift(Data,[0,1]) ; % shift 
        
        peak = CCpeakFinder(Data(:)) ;
        threshold = 10*peak/abs(peak) ;
            
        SpikeTrain{a} = zeros(size(Data)) ; % prep spike trains

        if threshold>0 ;
            SpikeTrain{a}(Data>=threshold & DataShift<threshold) = 1 ; % find indicies where change in current is above threshold
        else
            SpikeTrain{a}(Data<=threshold & DataShift>threshold) = 1 ;
        end
        SpikeTrain{a}(:,1) = 0 ;
        
        [prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ;
        [stmPnts, error] = ITCGetStmPts(epochs(1),0,0,fp) ;
    else
        SpikePnts = SpikeDetection_WC(Data,-30,1/SI(1)) ;
        
        % spike trains
        SpikeTrain{a} = zeros(size(time{a})) ;
        for b=1:length(SpikePnts) ;
            SpikeTrain{a}(b,SpikePnts{b}) = 1 ;
        end
        prePnts = 10000 ;
        stmPnts = 100000 ;
    end
    
%     figure % spike detection check
%     for b=1:size(Data,1) ; % on each trial
%         plot(time{a},Data(b,:))
%         hold on
%         plot(time{a},SpikeTrain{a}(b,:)*10,'r')
%         hold off
% 
%         text(0,.9,num2str(epochs(b)),'units','norm')
%         pause
%     end
    
    % psth
    SmoothPnts = floor(Smoothtime/SI(1)) ;
    
    for b=1:size(SpikeTrain{a},1) ; % this was changed 7/19/10, before psth was calculated by smoothing mean of spike train
        prePSTH(b,:) = smooth([zeros(1,floor(SmoothPnts/2)),SpikeTrain{a}(b,prePnts+1:prePnts+stmPnts),zeros(1,floor(SmoothPnts/2))],SmoothPnts)' ;  
        
    end
    PSTH_short = prePSTH(:,(floor(SmoothPnts/2)+1):length(prePSTH)-floor(SmoothPnts/2))/SI(1) ;

    PSTH_mean = mean(PSTH_short) ;

    PSTH_norm = PSTH_short/mean(PSTH_mean) ;

    PSTH_norm_mean{a} = mean(PSTH_norm) ;
    PSTH_norm_std{a} = std(PSTH_norm) ;
    PSTH_norm_cv(a) = nanmean(PSTH_norm_std{a}./PSTH_norm_mean{a}) ;
    PSTH_norm_cv_SEM(a) = nanstd(PSTH_norm_std{a}./PSTH_norm_mean{a})/sqrt(length(PSTH_norm_std{a})) ;
    
%     figure
%     plot(time{a}(prePnts+1:prePnts+stmPnts),PSTH{a})
    
%     %linear filter
%    LinearFilter{a} = LinFilterFinder(repmat(stim(1:stmPnts),size(SpikeTrain{a},1),1),SpikeTrain{a}(:,prePnts+1:prePnts+stmPnts), 1/SI(1),40) ;
    
%     figure
%     plot(time{a}(prePnts+1:prePnts+stmPnts),LinearFilter{a})
  
%     % signal to noise ratio
%     [powerX,meanSpikeSpectrum,resSpikeSpectrum,meanSpikeSpectrum_smth{a},resSpikeSpectrum_smth,snrSpikeSpectrum_smth{a},VarSumResSpikeSpectrum] = snrSpikeSpectrum(SpikeTrain{a},1/SI(1),.001) ;
%     powerMin = find(powerX>=1,1,'first') ;
%     powerMax = find(powerX<=20,1,'last') ;
%     SumSnr(a) = sum(meanSpikeSpectrum(powerMin:powerMax))/sum(resSpikeSpectrum(powerMin:powerMax)) ; 
    
%     figure
%     plot(meanSpikeSpectrum_smth.Freq,snrSpikeSpectrum_smth{a})
    
    clearvars -except Input id A PSTH_norm_mean PSTH_norm_std time LinearFilter SumSnr meanSpikeSpectrum_smth snrSpikeSpectrum_smth SI  PSTH_norm_cv PSTH_norm_cv_SEM
end

KStest =  kstest2(PSTH_norm_cv(1:4),PSTH_norm_cv(5:end)) ; % are average cvs significantly different between CA and DC? (1=yes, 0=no) 

PSTH_norm_CA = [] ;
PSTH_norm_DC = [] ;

PSTH_norm_std_CA = [] ;
PSTH_norm_std_DC = [] ;

% LinearFilter_CA = [] ;
% LinearFilter_DC = [] ;
% 
% LinearFilter_norm_CA = [] ;
% LinearFilter_norm_DC = [] ;
% 
% snrSpikeSpectrum_smth_CA = [] ;
% snrSpikeSpectrum_smth_DC = [] ;
% 
% SumSnr_CA = SumSnr(1:4) ;
% SumSnr_DC = SumSnr(5:11) ;

% auttocorrelations of psth
for a=1:11 ;
    ac{a} = xcov(PSTH_norm_mean{a},'coef') ;
    actime{a} = ([1:length(ac{a})] - (length(ac{a})+1)/2)*SI(1) ;
end


figure
PSTHFig = gcf ;

figure
AcFig = gcf ;

% figure
% LinearFilterFig = gcf ;
% 
% figure
% SnrFig = gcf ;


for a=1:4;
    figure(PSTHFig)
    plot(time{a}(1:length(PSTH_norm_mean{a})),PSTH_norm_mean{a},'b') % Freds CA light stim has 0.1 sec worth prestim
    hold on
    plot(time{a}(1:length(PSTH_norm_mean{a})),PSTH_norm_mean{a}+PSTH_norm_std{a},'b--')
    PSTH_norm_CA = [PSTH_norm_CA;PSTH_norm_mean{a}] ;
    PSTH_norm_std_CA = [PSTH_norm_std_CA;PSTH_norm_std{a}] ;
    
    figure(AcFig)
    plot(actime{a},ac{a},'b')
    
%     figure(LinearFilterFig)
%     plot(time{a}(1:length(LinearFilter{a})),LinearFilter{a},'b') % Freds CA light stim has 0.1 sec worth prestim
%     hold on
%     LinearFilter_CA = [LinearFilter_CA;LinearFilter{a}] ;
%     LinearFilter_norm_CA = [LinearFilter_norm_CA;LinearFilter{a}/max(LinearFilter{a})] ;
%     
%     figure(SnrFig)
%     plot(meanSpikeSpectrum_smth{a}.Freq,snrSpikeSpectrum_smth{a},'b') % Freds CA light stim has 0.1 sec worth prestim
%     hold on
%     snrSpikeSpectrum_smth_CA = [snrSpikeSpectrum_smth_CA;snrSpikeSpectrum_smth{a}] ;
    
end

for a=5:11;
    figure(PSTHFig)
    plot(time{a}(1:length(PSTH_norm_mean{a})),PSTH_norm_mean{a},'r') % Freds CA light stim has 0.1 sec worth prestim
    hold on
    plot(time{a}(1:length(PSTH_norm_mean{a})),PSTH_norm_mean{a}+PSTH_norm_std{a},'r--')
    PSTH_norm_DC = [PSTH_norm_DC;PSTH_norm_mean{a}] ;
    PSTH_norm_std_DC = [PSTH_norm_std_DC;PSTH_norm_std{a}] ;
    
    figure(AcFig)
    plot(actime{a},ac{a},'r')
    
    
%     figure(LinearFilterFig)
%     plot(time{a}(1:length(LinearFilter{a})),LinearFilter{a},'r') % Freds CA light stim has 0.1 sec worth prestim
%     hold on
%     LinearFilter_DC = [LinearFilter_DC;LinearFilter{a}] ;
%     LinearFilter_norm_DC = [LinearFilter_norm_DC;LinearFilter{a}/max(LinearFilter{a})] ;
%     
%     figure(SnrFig)
%     plot(meanSpikeSpectrum_smth{a}.Freq,snrSpikeSpectrum_smth{a},'r') % Freds CA light stim has 0.1 sec worth prestim
%     hold on
%     snrSpikeSpectrum_smth_DC = [snrSpikeSpectrum_smth_DC;snrSpikeSpectrum_smth{a}] ;
    
end



% PSTH_CA_mean = mean(PSTH_CA) ;
% PSTH_DC_mean = mean(PSTH_DC) ;
% 
% PSTH_CA_sem = std(PSTH_CA)/sqrt(4) ;
% PSTH_DC_sem = std(PSTH_DC)/sqrt(11) ;
% 
% PSTH_norm_CA_mean = mean(PSTH_norm_CA) ;
% PSTH_norm_DC_mean = mean(PSTH_norm_DC) ;
% 
% PSTH_norm_CA_sem = std(PSTH_norm_CA)/sqrt(4) ;
% PSTH_norm_DC_sem = std(PSTH_norm_DC)/sqrt(11) ;
% 
% 
% LinearFilter_CA_mean = mean(LinearFilter_CA) ;
% LinearFilter_DC_mean = mean(LinearFilter_DC) ;
% 
% LinearFilter_CA_sem = std(LinearFilter_CA)/sqrt(4) ;
% LinearFilter_DC_sem = std(LinearFilter_DC)/sqrt(11) ;
% 
% LinearFilter_norm_CA_mean = mean(LinearFilter_norm_CA) ;
% LinearFilter_norm_DC_mean = mean(LinearFilter_norm_DC) ;
% 
% LinearFilter_norm_CA_sem = std(LinearFilter_norm_CA)/sqrt(4) ;
% LinearFilter_norm_DC_sem = std(LinearFilter_norm_DC)/sqrt(11) ;
% 
% 
% snrSpikeSpectrum_smth_CA_mean = mean(snrSpikeSpectrum_smth_CA) ;
% snrSpikeSpectrum_smth_DC_mean = mean(snrSpikeSpectrum_smth_DC) ;
% 
% snrSpikeSpectrum_smth_CA_sem = std(snrSpikeSpectrum_smth_CA)/sqrt(4) ;
% snrSpikeSpectrum_smth_DC_sem = std(snrSpikeSpectrum_smth_DC)/sqrt(11) ;


% figure
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_CA_mean,'b')
% hold on
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_CA_mean+PSTH_CA_sem,'b--')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_CA_mean-PSTH_CA_sem,'b--')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_DC_mean(1:length(PSTH_CA_mean)),'r')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_DC_mean(1:length(PSTH_CA_mean))+PSTH_DC_sem(1:length(PSTH_CA_mean)),'r--')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_DC_mean(1:length(PSTH_CA_mean))-PSTH_DC_sem(1:length(PSTH_CA_mean)),'r--')
% 
% figure
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_norm_CA_mean,'b')
% hold on
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_norm_CA_mean+PSTH_norm_CA_sem,'b--')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_norm_CA_mean-PSTH_norm_CA_sem,'b--')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_norm_DC_mean(1:length(PSTH_CA_mean)),'r')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_norm_DC_mean(1:length(PSTH_CA_mean))+PSTH_norm_DC_sem(1:length(PSTH_CA_mean)),'r--')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_norm_DC_mean(1:length(PSTH_CA_mean))-PSTH_norm_DC_sem(1:length(PSTH_CA_mean)),'r--')
% 
% 
% figure
% plot(time{a}(1:length(LinearFilter_CA_mean)),LinearFilter_CA_mean,'b')
% hold on
% plot(time{a}(1:length(LinearFilter_CA_mean)),LinearFilter_CA_mean+LinearFilter_CA_sem,'b--')
% plot(time{a}(1:length(LinearFilter_CA_mean)),LinearFilter_CA_mean-LinearFilter_CA_sem,'b--')
% plot(time{a}(1:length(LinearFilter_CA_mean)),LinearFilter_DC_mean(1:length(LinearFilter_CA_mean)),'r')
% plot(time{a}(1:length(LinearFilter_CA_mean)),LinearFilter_DC_mean(1:length(LinearFilter_CA_mean))+LinearFilter_DC_sem(1:length(LinearFilter_CA_mean)),'r--')
% plot(time{a}(1:length(LinearFilter_CA_mean)),LinearFilter_DC_mean(1:length(LinearFilter_CA_mean))-LinearFilter_DC_sem(1:length(LinearFilter_CA_mean)),'r--')
% 
% figure
% plot(time{a}(1:length(LinearFilter_norm_CA_mean)),LinearFilter_norm_CA_mean,'b')
% hold on
% plot(time{a}(1:length(LinearFilter_norm_CA_mean)),LinearFilter_norm_CA_mean+LinearFilter_norm_CA_sem,'b--')
% plot(time{a}(1:length(LinearFilter_norm_CA_mean)),LinearFilter_norm_CA_mean-LinearFilter_norm_CA_sem,'b--')
% plot(time{a}(1:length(LinearFilter_norm_CA_mean)),LinearFilter_norm_DC_mean(1:length(LinearFilter_norm_CA_mean)),'r')
% plot(time{a}(1:length(LinearFilter_norm_CA_mean)),LinearFilter_norm_DC_mean(1:length(LinearFilter_norm_CA_mean))+LinearFilter_norm_DC_sem(1:length(LinearFilter_norm_CA_mean)),'r--')
% plot(time{a}(1:length(LinearFilter_norm_CA_mean)),LinearFilter_norm_DC_mean(1:length(LinearFilter_norm_CA_mean))-LinearFilter_norm_DC_sem(1:length(LinearFilter_norm_CA_mean)),'r--')
% 
% 
% figure
% plot(meanSpikeSpectrum_smth{1}.Freq,snrSpikeSpectrum_smth_CA_mean,'b')
% hold on
% plot(meanSpikeSpectrum_smth{1}.Freq,snrSpikeSpectrum_smth_CA_mean+snrSpikeSpectrum_smth_CA_sem,'b--')
% plot(meanSpikeSpectrum_smth{1}.Freq,snrSpikeSpectrum_smth_CA_mean-snrSpikeSpectrum_smth_CA_sem,'b--')
% plot(meanSpikeSpectrum_smth{11}.Freq,snrSpikeSpectrum_smth_DC_mean,'r')
% plot(meanSpikeSpectrum_smth{11}.Freq,snrSpikeSpectrum_smth_DC_mean+snrSpikeSpectrum_smth_DC_sem,'r--')
% plot(meanSpikeSpectrum_smth{11}.Freq,snrSpikeSpectrum_smth_DC_mean-snrSpikeSpectrum_smth_DC_sem,'r--')

% for a=1:4;
%     for b=1:7 ;
%         cc(a,b) = mean(PSTH_norm_CA(a,:).*PSTH_norm_DC(b,1:length(PSTH_CA_mean)))/sqrt(mean(PSTH_norm_CA(a,:).^2)*mean(PSTH_norm_DC(b,1:length(PSTH_CA_mean)).^2));
%     end
% end
% [rw,clmn] = find(cc==max(cc(:))) ;
% 
% for a=1:4;
%     cc2(a) = mean(PSTH_norm_CA(rw,:).*PSTH_norm_CA(a,:))/sqrt(mean(PSTH_norm_CA(rw,:).^2)*mean(PSTH_norm_CA(a,:).^2));
% end
% diffCAi = find(cc2==min(cc2)) ;
% 
% figure
% plot(time{1}(1:length(PSTH_CA_mean)),PSTH_norm_CA(rw,:),'b')
% hold on
% plot(time{1}(1:length(PSTH_CA_mean)),PSTH_norm_DC(clmn,1:length(PSTH_CA_mean)),'r')
% plot(time{1}(1:length(PSTH_CA_mean)),PSTH_norm_CA(diffCAi,:),'c')

% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_norm_DC_mean(1:length(PSTH_CA_mean))+PSTH_norm_DC_sem(1:length(PSTH_CA_mean)),'r--')
% plot(time{a}(1:length(PSTH_CA_mean)),PSTH_norm_DC_mean(1:length(PSTH_CA_mean))-PSTH_norm_DC_sem(1:length(PSTH_CA_mean)),'r--')

cai = 3 ;
cai2 = 4 ;
cai3 = 2 ;
dci = 5 ;


% ForIgor

identifier = ['PsthNormMidgetCA',num2str(A{cai})] ;
ForIgor.(identifier) = PSTH_norm_mean{cai}(1:length(PSTH_norm_mean{cai}));

identifier = ['PsthNormMidgetDC',num2str(A{dci})] ;
ForIgor.(identifier) = PSTH_norm_mean{dci}(1:length(PSTH_norm_mean{cai}));


identifier = ['PsthNormStdMidgetCA',num2str(A{cai})] ;
ForIgor.(identifier) = PSTH_norm_std{cai}(1:length(PSTH_norm_mean{cai}));

identifier = ['PsthNormPlusStdMidgetDC',num2str(A{dci})] ;
ForIgor.(identifier) = PSTH_norm_mean{dci}(1:length(PSTH_norm_mean{cai})) + PSTH_norm_std{dci}(1:length(PSTH_norm_mean{cai}));

identifier = ['PsthNormMinusStdMidgetDC',num2str(A{dci})] ;
ForIgor.(identifier) = PSTH_norm_mean{dci}(1:length(PSTH_norm_mean{cai})) - PSTH_norm_std{dci}(1:length(PSTH_norm_mean{cai}));


identifier = ['AcMidgetCA',num2str(A{cai})] ;
ForIgor.(identifier) = ac{cai} ;

identifier = ['AcMidgetCA',num2str(A{cai2})] ;
ForIgor.(identifier) = ac{cai2} ;

identifier = ['AcMidgetCA',num2str(A{cai3})] ;
ForIgor.(identifier) = ac{cai3} ;

identifier = ['AcMidgetDC',num2str(A{cai})] ;
ForIgor.(identifier) = ac{dci} ;

identifier = ['AcTimeMidgetCA',num2str(A{cai2})] ;
ForIgor.(identifier) = actime{cai} ;

identifier = ['AcTimeMidgetDC',num2str(A{cai})] ;
ForIgor.(identifier) = actime{dci} ;


identifier = 'timePsthNormMidget' ;
ForIgor.(identifier) = time{1}(1:length(PSTH_norm_mean{cai})) ;




    