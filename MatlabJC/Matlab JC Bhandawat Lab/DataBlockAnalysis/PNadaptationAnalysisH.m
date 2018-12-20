function ForIgor = PNadaptationAnalysisH(Input,A)

% this function will analyze PN pulse responses +/- a background odor that
% is triggered during the trial
% similar to G but uses different points for voltage and psth max and min response values
% Input.id1 should be a cell array with trial numbers grouped according to vector specified in Input.id2.
% JC 1/8/13


% parameters
sampleRate = 10000 ; % (hz) temp hard coded - should be saved in file 
driftCheckTime = 0.25 ; %(sec) time at begining and end which current injected is inspected for changes
bgTransTime = 5 ; %(sec) time after background odor starts to avoid transient response
PsthBinTime = 0.1 ; % (sec) time width of psth bin
OdorRspTime = .1 ; % (sec) time before and after peak of mean psth which odor response is assessed

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
rootdir = ['Z:\Cafaro Data Backup\', Input(A).cellname(1:6),'Data'];

Concentrations = str2num(Input(A).(id2)) ;
NumConcentrations = length(Concentrations) ;

BgConcentrations = Input(A).(id3) ;
NumBackgrounds = length(BgConcentrations) ;

for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
    
        odorRspTrials{a}{b} = str2num(Input(A).(id1){a}{b}) ;
        NumTrials(a,b) = length(odorRspTrials{a}{b}) ;
        loopNum = 0 ;
        for c = 1:NumTrials(a,b) ; % for each trial
            loopNum = loopNum+1 ;

            temp = load([rootdir,'\',Input(A).cellname,'\','voltage_',Input(A).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
            vData{a}{b}(loopNum,:) = temp.voltage ;

            temp = load([rootdir,'\',Input(A).cellname,'\','current_',Input(A).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
            iData{a}{b}(loopNum,:) = temp.current ;

            temp = load([rootdir,'\',Input(A).cellname,'\','Ao0_',Input(A).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
            ao0Data{a}{b}(loopNum,:) = temp.Ao0 ; % odor

            temp = load([rootdir,'\',Input(A).cellname,'\','Ao1_',Input(A).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
            ao1Data{a}{b}(loopNum,:) = temp.Ao1 ; % background valve   

            temp = load([rootdir,'\',Input(A).cellname,'\','TrigTime_',Input(A).cellname,'_',num2str(odorRspTrials{a}{b}(c))]) ;
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
try PreviousDetect = load([spikeDataPath,'cell',num2str(A)]) ; 
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

    save([spikeDataPath,'cell',num2str(A)], 'spikePnt','spikeDetectionParameters') ; % save spike data
end


% voltage data
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            vData_mean{a}{b} = mean(vData{a}{b}) ; % mean voltage trace
        end
    end
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
            [m,mi] = max(vData_mean{a}{b}(iopb:iopb+sampleRate)) ; % max point of mean vData within a second of odor pulse onset
            Viorb(a,b) = mi-1+iopb - OdorRspTime*sampleRate ; % point of odor pulse begining
            Viore(a,b) = mi-1+iopb + OdorRspTime*sampleRate ; % point of odor pulse end
            
            [m,mi] = min(vData_mean{a}{b}(iope:iope+sampleRate)) ; % min point of mean vData within a second of odor pulse offset
            Viohb(a,b) = mi-1+iope - OdorRspTime*sampleRate ; % point of odor pulse begining
            Viohe(a,b) = mi-1+iope + OdorRspTime*sampleRate ; % point of odor pulse end            

            [m,mi] = max(Psth_mean{a}{b}(iopb:iopb+sampleRate)) ; % max point of mean psth within a second of odor pulse onset
            SRiorb(a,b) = mi-1+iopb - OdorRspTime*sampleRate ; % point of odor pulse begining
            SRiore(a,b) = mi-1+iopb + OdorRspTime*sampleRate ; % point of odor pulse end
            
            [m,mi] = min(Psth_mean{a}{b}(iope:iope+sampleRate)) ; % min point of mean psth within a second of odor pulse offset
            SRiohb(a,b) = mi-1+iope - OdorRspTime*sampleRate ; % point of odor pulse begining
            SRiohe(a,b) = mi-1+iope + OdorRspTime*sampleRate ; % point of odor pulse end
        end
    end
end

% assess resting potential 
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                Vrest{a}{b}(c) = mean(vData{a}{b}(c,1:ibpb)) ; % mV (start of trial: begining of background odor pulse
            end
            Vrest_mean(a,b) = mean(Vrest{a}{b}) ;
            Vrest_std(a,b) = std(Vrest{a}{b}) ;
            Vrest_sem(a,b) = Vrest_std(a,b)/sqrt(NumTrials(a,b)) ;
        else
            Vrest_mean(a,b) = nan ;
            Vrest_std(a,b) = nan ;
            Vrest_sem(a,b) = nan ;
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
    
% assess background transient potential
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                Vtbg{a}{b}(c) = max(vData{a}{b}(c,ibpb:iopb)) ; % mV (start of background: begining of odor pulse)
            end
            Vtbg_mean(a,b) = mean(Vtbg{a}{b}) ;
            Vtbg_std(a,b) = std(Vtbg{a}{b}) ;
            Vtbg_sem(a,b) = Vtbg_std(a,b)/sqrt(NumTrials(a,b)) ;
        else
            Vtbg_mean(a,b) = nan ;
            Vtbg_std(a,b) = nan ;
            Vtbg_sem(a,b) = nan ;
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

% assess background potential (during background odor after transient response)
bgTransPnts = bgTransTime*sampleRate ;
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                Vbg{a}{b}(c) = mean(vData{a}{b}(c,ibpb+bgTransPnts:iopb)) ; % mV (start of background + transient time: begining of odor pulse)
            end
            Vbg_mean(a,b) = mean(Vbg{a}{b}) ;
            Vbg_std(a,b) = std(Vbg{a}{b}) ;
            Vbg_sem(a,b) = Vbg_std(a,b)/sqrt(NumTrials(a,b)) ;
        else
            Vbg_mean(a,b) = nan ;
            Vbg_std(a,b) = nan ;
            Vbg_sem(a,b) = nan ;
        end
    end
end

% assess background spike rate (during background odor after transient response)
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

% odor pulse response potential
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                Vpulse{a}{b}(c) = mean(vData{a}{b}(c,Viorb(a,b):Viore(a,b))) ; %mV (odor response depol start: odor response depol end)
                VpulseHyp{a}{b}(c) = mean(vData{a}{b}(c,Viohb(a,b):Viohe(a,b))) ; %mV (odor response hyp start: odor response hyp end)
            end
            
            Vpulse_mean(a,b) = mean(Vpulse{a}{b}) ; %mV
            Vpulse_std(a,b) = std(Vpulse{a}{b}) ; %mV
            Vpulse_sem(a,b) = Vpulse_std(a,b)/sqrt(NumTrials(a,b)) ;
            
            VpulseHyp_mean(a,b) = mean(VpulseHyp{a}{b}) ; %mV
            VpulseHyp_std(a,b) = std(VpulseHyp{a}{b}) ; %mV
            VpulseHyp_sem(a,b) = VpulseHyp_std(a,b)/sqrt(NumTrials(a,b)) ;
        else
            Vpulse_mean(a,b) = nan ; %mV
            Vpulse_std(a,b) = nan ; %mV
            Vpulse_sem(a,b) = nan ;
            
            VpulseHyp_mean(a,b) = nan ; %mV
            VpulseHyp_std(a,b) = nan ; %mV
            VpulseHyp_sem(a,b) = nan ;
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

% delta potential (pulse minus background)
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            VpulseMinBg{a}{b} = Vpulse{a}{b}- Vbg{a}{b} ; % mV
            VpulseMinBg_mean(a,b) = mean(VpulseMinBg{a}{b}) ; % mV
            VpulseMinBg_std(a,b) = std(VpulseMinBg{a}{b}) ; % mV
            VpulseMinBg_sem(a,b) = VpulseMinBg_std(a,b)/sqrt(NumTrials(a,b)) ;
            
            VpulseHypMinBg{a}{b} = VpulseHyp{a}{b}- Vbg{a}{b} ; % mV
            VpulseHypMinBg_mean(a,b) = mean(VpulseHypMinBg{a}{b}) ; % mV
            VpulseHypMinBg_std(a,b) = std(VpulseHypMinBg{a}{b}) ; % mV
            VpulseHypMinBg_sem(a,b) = VpulseHypMinBg_std(a,b)/sqrt(NumTrials(a,b)) ;
        else
            VpulseMinBg_mean(a,b) = nan ;
            VpulseMinBg_std(a,b) = nan ;
            VpulseMinBg_sem(a,b) = nan ;
            
            VpulseHypMinBg_mean(a,b) = nan ;
            VpulseHypMinBg_std(a,b) = nan ; 
            VpulseHypMinBg_sem(a,b) = nan ;
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

% delta potential (pulse minus background) NORMALIZED
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            VpulseMinBg_norm{a}{b} = VpulseMinBg{a}{b}/max(VpulseMinBg_mean(:)) ; % mV
            VpulseMinBg_norm_mean(a,b) = mean(VpulseMinBg_norm{a}{b}) ; % mV
            VpulseMinBg_norm_std(a,b) = std(VpulseMinBg_norm{a}{b}) ; % mV
             VpulseMinBg_norm_sem(a,b) =  VpulseMinBg_norm_std(a,b)/sqrt(NumTrials(a,b)) ;
            
            VpulseHypMinBg_norm{a}{b} = VpulseHypMinBg{a}{b}/min(VpulseHypMinBg_mean(:)) ; % mV
            VpulseHypMinBg_norm_mean(a,b) = mean(VpulseHypMinBg_norm{a}{b}) ; % mV
            VpulseHypMinBg_norm_std(a,b) = std(VpulseHypMinBg_norm{a}{b}) ; % mV
            VpulseHypMinBg_norm_sem(a,b) = VpulseHypMinBg_norm_std(a,b)/sqrt(NumTrials(a,b)) ;
        else
            VpulseMinBg_norm_mean(a,b) = nan ;
            VpulseMinBg_norm_std(a,b) = nan ;
            VpulseMinBg_norm_sem(a,b) = nan ;
            
            VpulseHypMinBg_norm_mean(a,b) = nan ;
            VpulseHypMinBg_norm_std(a,b) = nan ; 
            VpulseHypMinBg_norm_sem(a,b) = nan ;
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

% delta potential (pulse minus rest)
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            VpulseMinRest{a}{b} = Vpulse{a}{b}- Vrest{a}{b} ; % mV
            VpulseMinRest_mean(a,b) = mean(VpulseMinRest{a}{b}) ; % mV
            VpulseMinRest_std(a,b) = std(VpulseMinRest{a}{b}) ; % mV
            VpulseMinRest_sem(a,b) = VpulseMinRest_std(a,b)/sqrt(NumTrials(a,b)) ;
            
            VpulseHypMinRest{a}{b} = VpulseHyp{a}{b}- Vrest{a}{b} ; % mV
            VpulseHypMinRest_mean(a,b) = mean(VpulseHypMinRest{a}{b}) ; % mV
            VpulseHypMinRest_std(a,b) = std(VpulseHypMinRest{a}{b}) ; % mV
            VpulseHypMinRest_sem(a,b) = VpulseHypMinRest_std(a,b)/sqrt(NumTrials(a,b)) ;
        else
            VpulseMinRest_mean(a,b) = nan ;
            VpulseMinRest_std(a,b) = nan ;
            VpulseMinRest_sem(a,b) = nan ;
            
            VpulseHypMinRest_mean(a,b) = nan ;
            VpulseHypMinRest_std(a,b) = nan ;
            VpulseHypMinRest_sem(a,b) = nan ;
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

% delta potential (pulse minus rest) NORMALIZED
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            VpulseMinRest_norm{a}{b} = VpulseMinRest{a}{b}/max(VpulseMinRest_mean(:)) ; % mV
            VpulseMinRest_norm_mean(a,b) = mean(VpulseMinRest_norm{a}{b}) ; % mV
            VpulseMinRest_norm_std(a,b) = std(VpulseMinRest_norm{a}{b}) ; % mV
            VpulseMinRest_norm_sem(a,b) =  VpulseMinRest_norm_std(a,b)/sqrt(NumTrials(a,b)) ;
            
            VpulseHypMinRest_norm{a}{b} = VpulseHypMinRest{a}{b}/min(VpulseHypMinRest_mean(:)) ; % mV
            VpulseHypMinRest_norm_mean(a,b) = mean(VpulseHypMinRest_norm{a}{b}) ; % mV
            VpulseHypMinRest_norm_std(a,b) = std(VpulseHypMinRest_norm{a}{b}) ; % mV
            VpulseHypMinRest_norm_sem(a,b) = VpulseHypMinRest_norm_std(a,b)/sqrt(NumTrials(a,b)) ;
        else
            VpulseMinRest_norm_mean(a,b) = nan ;
            VpulseMinRest_norm_std(a,b) = nan ;
            VpulseMinRest_norm_sem(a,b) = nan ;
            
            VpulseHypMinRest_norm_mean(a,b) = nan ;
            VpulseHypMinRest_norm_std(a,b) = nan ;
            VpulseHypMinRest_norm_sem(a,b) = nan ;
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
                vDataMinBg{a}{b}(c,:) = vData{a}{b}(c,:) - Vbg{a}{b}(c) ; 
                PsthMinBg{a}{b}(c,:) = Psth{a}{b}(c,:) - SRbg{a}{b}(c) ;
            end
        vDataMinBg_mean{a}{b} = mean(vDataMinBg{a}{b},1) ;
        PsthMinBg_mean{a}{b} = mean(PsthMinBg{a}{b},1) ;
        
        vDataMinBg_mean_norm{a}{b} = vDataMinBg_mean{a}{b}/max(vDataMinBg_mean{a}{b}(iopb:ibpe)) ; % data normalized to pulse depol amp
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
%     
% % spike detection
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         for c = 1:NumTrials(a,b) ;
%             plot(time,vData{a}{b}(c,:)) 
%             hold on
%             plot(time(spikePnt{a}{b}{c}),vData{a}{b}(c,spikePnt{a}{b}{c}),'r*')
%             title(num2str(odorRspTrials{a}{b}(c)))
%             hold off
%             pause
%         end
%     end
% end
%       
% % spike raster
% figure
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         for c = 1:NumTrials(a,b) ;
%             for d=1:length(spikePnt{a}{b}{c}) ;
%                 plot([1,1]*spikePnt{a}{b}{c}(d),[odorRspTrials{a}{b}(c)-1,odorRspTrials{a}{b}(c)],'Color',colorMat{b}(a,:))
%                 hold on
%             end
%         end
%     end
% end
% 
% %mean voltage and psth for different concentrations and backgrounds and time
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             subplot(NumConcentrations+2,2,b*2-1)
%             plot(time,vData_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             
%             ylabel(num2str(Concentrations(b)))
%             hold on
%             
%             subplot(NumConcentrations+2,2,b*2)
%             plot(time,Psth_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             
%             hold on
%             
%         end
%     end
% end
% 
% subplot(NumConcentrations+2,2,2*NumConcentrations+1:2*NumConcentrations+2)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(tDataN{a}{b},Vrest{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%             plot(tDataN{a}{b},Vbg{a}{b},'.','Color',colorMat{b}(a,:))
%             plot(tDataN{a}{b},Vpulse{a}{b},'o','Color',colorMat{b}(a,:))
%         end
%     end
% end
% ylabel('potential (mV)')
% 
% subplot(NumConcentrations+2,2,2*NumConcentrations+3:2*NumConcentrations+4)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(tDataN{a}{b},SRrest{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%             plot(tDataN{a}{b},SRbg{a}{b},'.','Color',colorMat{b}(a,:))
%             plot(tDataN{a}{b},SRpulse{a}{b},'o','Color',colorMat{b}(a,:))
%         end
%     end
% end
% xlabel('trig time (sec)')
% ylabel('spike rate (hz)')
% 
% axes('Position',[0 0 .02 1],'Visible','off');
% text(0,.01,RepStat,'FontSize',5)
% 
% %mean voltage and psth with background subtracted for different concentrations and backgrounds and time
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             subplot(NumConcentrations+2,2,b*2-1)
%             plot(time,vDataMinBg_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             
%             ylabel(num2str(Concentrations(b)))
%             hold on
%             %plot(time(1,[iorb(a,b),iore(a,b)]),vData_mean{a}{b}(1,[iorb(a,b),iore(a,b)]),'o','Color',colorMat{NumConcentrations}(a,:))
%             
%             subplot(NumConcentrations+2,2,b*2)
%             plot(time,PsthMinBg_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             
%             hold on
%             %plot(time(1,[iorb(a,b),iore(a,b)]),Psth_mean{a}{b}(1,[iorb(a,b),iore(a,b)]),'o','Color',colorMat{NumConcentrations}(a,:))
%             
%         end
%     end
% end
% 
% subplot(NumConcentrations+2,2,2*NumConcentrations+1:2*NumConcentrations+2)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(tDataN{a}{b},Vrest{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%             plot(tDataN{a}{b},VpulseMinBg{a}{b},'+','Color',colorMat{b}(a,:))
%         end
%     end
% end
% ylabel('potential (mV)')
% 
% subplot(NumConcentrations+2,2,2*NumConcentrations+3:2*NumConcentrations+4)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(tDataN{a}{b},SRrest{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%             plot(tDataN{a}{b},SRpulseMinBg{a}{b},'+','Color',colorMat{b}(a,:))
%         end
%     end
% end
% xlabel('trig time (sec)')
% ylabel('spike rate (hz)')
% 
% % comparing kinetics of pulse responses
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             subplot(2,NumConcentrations,b)
%             plot(time,vDataMinBg_mean_norm{a}{b},'Color',colorMat{b}(a,:))
%             hold on
%             xlim([time(iopb),time(ibpe)])
% 
%             subplot(2,NumConcentrations,b+NumConcentrations)
%             plot(time,PsthMinBg_mean_norm{a}{b},'Color',colorMat{b}(a,:))
%             hold on
%             xlim([time(iopb),time(ibpe)])
%         end
%     end
% end
%  
% 
% % comparing rest and background data across concentrations
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             subplot(2,2,1)
%             plot(time,vData_mean{a}{b},'Color',colorMat{b}(a,:)) 
%             hold on
% 
%             subplot(2,2,2)
%             plot(time,Psth_mean{a}{b},'Color',colorMat{b}(a,:))            
%             hold on     
%             
%             subplot(2,2,3)
%             plot(time,vData_mean{a}{b}-Vrest_mean(a,b),'Color',colorMat{b}(a,:)) 
%             hold on
% 
%             subplot(2,2,4)
%             plot(time,Psth_mean{a}{b}-SRrest_mean(a,b),'Color',colorMat{b}(a,:))            
%             hold on                 
%             
%         end
%     end
% end
% 
% % comparing background responses with pulse response
% figure
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             subplot(3,2,1)
%             plot(Vbg{a}{b}-Vrest{a}{b},Vpulse{a}{b}-Vbg{a}{b},'+','Color',colorMat{b}(a,:)) 
%             xlabel('Vbg - Vrest (mV)')
%             ylabel('Vpulse - Vbg (hz)')
%             hold on
% 
%             subplot(3,2,2)
%             plot(SRbg{a}{b}-SRrest{a}{b},SRpulse{a}{b}-SRbg{a}{b},'+','Color',colorMat{b}(a,:))
%             xlabel('SRbg - SRrest (mV)')
%             ylabel('SRpulse - SRbg (hz)')
%             hold on        
% 
%             subplot(3,2,3)
%             plot(Vtbg{a}{b}-Vbg{a}{b},Vpulse{a}{b}-Vbg{a}{b},'+','Color',colorMat{b}(a,:)) 
%             xlabel('Vtbg - Vbg (mV)')
%             ylabel('Vpulse - Vbg (hz)')
%             hold on
%             
%             subplot(3,2,4)
%             plot(SRtbg{a}{b}-SRbg{a}{b},SRpulse{a}{b}-SRbg{a}{b},'+','Color',colorMat{b}(a,:)) 
%             xlabel('SRtbg - SRbg (mV)')
%             ylabel('SRpulse - SRbg (hz)')
%             hold on            
%  
%             subplot(3,2,5)
%             plot(Vtbg{a}{b}-Vrest{a}{b},Vpulse{a}{b}-Vbg{a}{b},'+','Color',colorMat{b}(a,:)) 
%             xlabel('Vtbg - Vrest (mV)')
%             ylabel('Vpulse - Vbg (hz)')
%             hold on
%             
%             subplot(3,2,6)
%             plot(SRtbg{a}{b}-SRrest{a}{b},SRpulse{a}{b}-SRbg{a}{b},'+','Color',colorMat{b}(a,:)) 
%             xlabel('SRtbg - SRrest (mV)')
%             ylabel('SRpulse - SRbg (hz)')
%             hold on              
%         end
%     end
% end
% 
% % comparing spikerate and potential (spont and evoked)
% figure 
% subplot(3,1,1)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(Vrest{a}{b},SRrest{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%             plot(Vbg{a}{b},SRbg{a}{b},'.','Color',colorMat{b}(a,:))
%             plot(Vpulse{a}{b},SRpulse{a}{b},'o','Color',colorMat{b}(a,:))
%         end
%     end
% end
% xlabel('potential (mV)')
% ylabel('spike rate (hz)')
% 
% % comparing spikerate and potential (change in spont to evoked)
% subplot(3,1,2) 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(VpulseMinBg{a}{b},SRpulseMinBg{a}{b},'+','Color',colorMat{b}(a,:))
%             hold on
%         end
%     end
% end
% xlabel('delta potential (mV)')
% ylabel('delta spike rate (hz)')
% 
% % comparing odor evoked spikerate and vrest
% subplot(3,1,3)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(Vrest{a}{b},SRpulseMinBg{a}{b},'+','Color',colorMat{b}(a,:))
%             hold on
%         end
%     end
% end
% xlabel('rest potential (mV)')
% ylabel('delta spike rate (hz)')
 
% potential as a function of concentration
figure
subplot(3,2,1)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),Vrest_mean(a,:),Vrest_std(a,:),Vrest_std(a,:),'*:','Color',colorMat{NumConcentrations}(a,:))
    hold on
    errorbar(log10(Concentrations),Vbg_mean(a,:),Vbg_std(a,:),Vbg_std(a,:),'*--','Color',colorMat{NumConcentrations}(a,:))
    errorbar(log10(Concentrations),Vpulse_mean(a,:),Vpulse_std(a,:),Vpulse_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
end
ylabel('potential')
xlabel('log concentration')

subplot(3,2,3)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),VpulseMinBg_norm_mean(a,:),VpulseMinBg_norm_std(a,:),VpulseMinBg_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
    hold on
end
ylabel('pulse min bg potential')
xlabel('log concentration')

subplot(3,2,5)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),VpulseMinRest_norm_mean(a,:),VpulseMinRest_norm_std(a,:),VpulseMinRest_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
    hold on
end
ylabel('pulse min rest potential')
xlabel('log concentration')

% hyperpolerization
subplot(3,2,2)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),Vrest_mean(a,:),Vrest_std(a,:),Vrest_std(a,:),'*:','Color',colorMat{NumConcentrations}(a,:))
    hold on
    errorbar(log10(Concentrations),Vbg_mean(a,:),Vbg_std(a,:),Vbg_std(a,:),'*--','Color',colorMat{NumConcentrations}(a,:))
    errorbar(log10(Concentrations),VpulseHyp_mean(a,:),VpulseHyp_std(a,:),VpulseHyp_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
end
ylabel('potential')
xlabel('log concentration')

subplot(3,2,4)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),VpulseHypMinBg_norm_mean(a,:),VpulseHypMinBg_norm_std(a,:),VpulseHypMinBg_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
    hold on
end
ylabel('pulse min bg potential')
xlabel('log concentration')

subplot(3,2,6)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),VpulseHypMinRest_norm_mean(a,:),VpulseHypMinRest_norm_std(a,:),VpulseHypMinRest_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
    hold on
end
ylabel('pulse min rest potential')
xlabel('log concentration')

%  spikes as a function of concentration
figure
subplot(3,2,1)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),SRrest_mean(a,:),SRrest_std(a,:),SRrest_std(a,:),'*:','Color',colorMat{NumConcentrations}(a,:))
    hold on
    errorbar(log10(Concentrations),SRbg_mean(a,:),SRbg_std(a,:),SRbg_std(a,:),'*--','Color',colorMat{NumConcentrations}(a,:))
    errorbar(log10(Concentrations),SRpulse_mean(a,:),SRpulse_std(a,:),SRpulse_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
end
ylabel('spike rate')
xlabel('log concentration')

subplot(3,2,3)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),SRpulseMinBg_norm_mean(a,:),SRpulseMinBg_norm_std(a,:),SRpulseMinBg_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
    hold on
end
ylabel('pulse min bg spike rate')
xlabel('log concentration')

subplot(3,2,5)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),SRpulseMinRest_norm_mean(a,:),SRpulseMinRest_norm_std(a,:),SRpulseMinRest_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
    hold on
end
ylabel('pulse min rest spike rate')
xlabel('log concentration')


% hyperpol
subplot(3,2,2)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),SRrest_mean(a,:),SRrest_std(a,:),SRrest_std(a,:),'*:','Color',colorMat{NumConcentrations}(a,:))
    hold on
    errorbar(log10(Concentrations),SRbg_mean(a,:),SRbg_std(a,:),SRbg_std(a,:),'*--','Color',colorMat{NumConcentrations}(a,:))
    errorbar(log10(Concentrations),SRpulseHyp_mean(a,:),SRpulseHyp_std(a,:),SRpulseHyp_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
end
ylabel('spike rate')
xlabel('log concentration')

subplot(3,2,4)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),SRpulseHypMinBg_norm_mean(a,:),SRpulseHypMinBg_norm_std(a,:),SRpulseHypMinBg_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
    hold on
end
ylabel('pulse min bg spike rate')
xlabel('log concentration')

subplot(3,2,6)
for a = 1:NumBackgrounds ;
    errorbar(log10(Concentrations),SRpulseHypMinRest_norm_mean(a,:),SRpulseHypMinRest_norm_std(a,:),SRpulseHypMinRest_norm_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
    hold on
end
ylabel('pulse min rest spike rate')
xlabel('log concentration')


% for igor

% % this is standard
% identifier = ['LogConcentration','cell',num2str(A)] ;
% ForIgor.(identifier) = log10(Concentrations) ;
% 
% for a = 1:NumBackgrounds ; % for each background, 
%     identifier = ['SRminBgMean','Bg',BgConcentrations{a},'cell',num2str(A)] ;
%     ForIgor.(identifier) = SRpulseMinBg_norm_mean(a,:) ;
%     
%     identifier = ['SRminBgSem','Bg',BgConcentrations{a},'cell',num2str(A)] ;
%     ForIgor.(identifier) = SRpulseMinBg_norm_sem(a,:) ;  
%     
%     identifier = ['VminBgMean','Bg',BgConcentrations{a},'cell',num2str(A)] ;
%     ForIgor.(identifier) = VpulseMinBg_norm_mean(a,:) ;
%     
%     identifier = ['VminBgSem','Bg',BgConcentrations{a},'cell',num2str(A)] ;
%     ForIgor.(identifier) = VpulseMinBg_norm_sem(a,:) ; 
% end

% this is not standard and expands data onto a common concentration axis to make it easier to average across cells
ConcentrationRange = [-8:-3] ;
identifier = ['LogConcentration','cell',num2str(A)] ;
ForIgor.(identifier) = ConcentrationRange ;

[c,ri] = intersect(ConcentrationRange,log10(Concentrations)) ;

for a = 1:NumBackgrounds ; % for each background, 
    identifier = ['SRminBgMean','Bg',BgConcentrations{a},'cell',num2str(A)] ;
    ForIgor.(identifier) = nan(1,length(ConcentrationRange)) ;
    ForIgor.(identifier)(ri) = SRpulseMinBg_norm_mean(a,:) ;
    
    identifier = ['SRminBgSem','Bg',BgConcentrations{a},'cell',num2str(A)] ;
    ForIgor.(identifier) = nan(1,length(ConcentrationRange)) ;
    ForIgor.(identifier)(ri) = SRpulseMinBg_norm_sem(a,:) ;  
    
    identifier = ['VminBgMean','Bg',BgConcentrations{a},'cell',num2str(A)] ;
    ForIgor.(identifier) = nan(1,length(ConcentrationRange)) ;
    ForIgor.(identifier)(ri) = VpulseMinBg_norm_mean(a,:) ;
    
    identifier = ['VminBgSem','Bg',BgConcentrations{a},'cell',num2str(A)] ;
    ForIgor.(identifier) = nan(1,length(ConcentrationRange)) ;
    ForIgor.(identifier)(ri) = VpulseMinBg_norm_sem(a,:) ; 
end

% unnormalized for pairwise comparison (NOT FINISHED- CAUTION!)
for a = 1:NumBackgrounds ; 
    for b=1:NumConcentrations ;
        identifier = ['SRminBgMean','Bg',BgConcentrations{a},'Pulse',num2str(abs(log10(Concentrations(b)))),'cell',num2str(A)] ;
        ForIgor.(identifier) = nan(1,length(ConcentrationRange)) ;
        ForIgor.(identifier)(ri) = SRpulseMinBg_mean(a,b) ;

        identifier = ['SRminBgSem','Bg',BgConcentrations{a},'Pulse',num2str(abs(log10(Concentrations(b)))),'cell',num2str(A)] ;
        ForIgor.(identifier) = nan(1,length(ConcentrationRange)) ;
        ForIgor.(identifier)(ri) = SRpulseMinBg_sem(a,:) ;  

        identifier = ['VminBgMean','Bg',BgConcentrations{a},'Pulse',num2str(abs(log10(Concentrations(b)))),'cell',num2str(A)] ;
        ForIgor.(identifier) = nan(1,length(ConcentrationRange)) ;
        ForIgor.(identifier)(ri) = VpulseMinBg_mean(a,:) ;

        identifier = ['VminBgSem','Bg',BgConcentrations{a},'Pulse',num2str(abs(log10(Concentrations(b)))),'cell',num2str(A)] ;
        ForIgor.(identifier) = nan(1,length(ConcentrationRange)) ;
        ForIgor.(identifier)(ri) = VpulseMinBg_sem(a,:) ; 
    end
end

for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration 
        if NumTrials(a,b)>0 ;
            identifier = ['Psth','LogPulse',num2str(abs(log10(Concentrations(b)))),'Bg',BgConcentrations{a},'cell',num2str(A)] ;
            ForIgor.(identifier) = Psth_mean{a}{b} ;

            identifier = ['PsthMinBg','LogPulse',num2str(abs(log10(Concentrations(b)))),'Bg',BgConcentrations{a},'cell',num2str(A)] ;
            ForIgor.(identifier) = PsthMinBg_mean{a}{b} ;        

            identifier = ['vData','LogPulse',num2str(abs(log10(Concentrations(b)))),'Bg',BgConcentrations{a},'cell',num2str(A)] ;
            ForIgor.(identifier) = vData_mean{a}{b} ;

            identifier = ['vDataMinBg','LogPulse',num2str(abs(log10(Concentrations(b)))),'Bg',BgConcentrations{a},'cell',num2str(A)] ;
            ForIgor.(identifier) = vDataMinBg_mean{a}{b} ;   
        end
    end
end

identifier = ['time','cell',num2str(A)] ;
ForIgor.(identifier) = time;  




