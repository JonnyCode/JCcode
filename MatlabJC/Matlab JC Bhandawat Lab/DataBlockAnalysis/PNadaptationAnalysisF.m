function ForIgor = PNadaptationAnalysisF(Input,A)

% this function is adapted from "PNadaptationAnalysisE", but uses an odor response time based on peak psth response.
% Input.id1 should be a cell array with trial numbers grouped according to vector specified in Input.id2.
% JC 10/28/12 


% parameters
sampleRate = 10000 ; % (hz) temp hard coded - should be saved in file 
driftCheckTime = 0.25 ; %(sec) time at begining and end which current injected is inspected for changes
absRefTime = 0.002 ; % (sec) 
minRise = .3 ; % (mV)
minFall = .15 ; % (mV)
PsthBinTime = 0.1 ; % (sec) time width of psth bin
OdorRspTime = .1 ; % (sec) time before and after peak of mean psth which odor response is assessed


id1 = 'OdorRsp' ;
id2 = 'OdorConcentration' ;
id3 = 'OdorBgTimes' ;

% load data in matricies
rootdir = ['Z:\Cafaro Data Backup\', Input(A).cellname(1:6),'Data'];

NumBackgrounds = length(Input(A).(id1)) ;

Concentrations = str2num(Input(A).(id2)) ;
NumConcentrations = length(Concentrations) ;

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
            ao1Data{a}{b}(loopNum,:) = temp.Ao1 ; % vext   

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

% round Vext pulse to nearest 10 mV and get rid of single sample pulses
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        for c = 1:NumTrials(a,b) ; % for each trial
            ao1DataR{a}{b}(c,:) = round(ao1Data{a}{b}(c,:)*100)/100 ; 
            
            for d=1:length(ao1DataR{a}{b}(c,:))-2 ;
                if ao1DataR{a}{b}(c,d)~=ao1DataR{a}{b}(c,d+1) && ao1DataR{a}{b}(c,d)==ao1DataR{a}{b}(c,d+2);
                    ao1DataR{a}{b}(c,d+1)= ao1DataR{a}{b}(c,d) ;
                end
            end
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

% make sure the R input check was at the same time
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            ao1DataRdiff = ao1DataR{a}{b} - repmat(ao1DataR{1}{1}(1,:),NumTrials(a,b),1) ;
            if sum(abs(ao1DataRdiff(:)))~=0 ;
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
                disp(['significant I drift in trial',num2str(c)]) ;
            end
        end
    end
end

% detect spikes in voltage data
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            [TempSpikePnt,SpikeData,NonSpikeData] = spikeFinder(vData{a}{b},sampleRate,absRefTime,minRise,minFall) ;
            spikePnt{a}{b}= TempSpikePnt ;
        end
    end
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
iopb = find(ao0DataB{1}{1}(1,:)~=0,1,'first')-1 ; % odor pulse begining
iope = find(ao0DataB{1}{1}(1,:)~=0,1,'last') ; % odor pulse ending

iipb = find(ao1DataR{1}{1}(1,:)~=0,1,'first') ; % current pulse beginging
iipe = find(ao1DataR{1}{1}(1,1:iopb)~=0,1,'last') ; % current pulse end

% odor response time
for a = 1:NumBackgrounds ;
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            [m,mi] = max(Psth_mean{a}{b}(iopb:iopb+sampleRate)) ; % max point of mean psth within a second of odor pulse onset
            iorb(a,b) = mi-1+iopb - OdorRspTime*sampleRate ; % point of odor pulse begining
            iore(a,b) = mi-1+iopb + OdorRspTime*sampleRate ; % point of odor pulse end
        end
    end
end

% assess resting potential
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                Vrest{a}{b}(c) = mean(vData{a}{b}(c,iipe:iopb)) ; % mV (end of I pulse: begining of odor pulse
            end
            Vrest_mean(a,b) = mean(Vrest{a}{b}) ;
            Vrest_std(a,b) = std(Vrest{a}{b}) ;
        else
            Vrest_mean(a,b) = nan ;
            Vrest_std(a,b) = nan ;
        end
    end
end

% spontaneous spike rate

for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration   
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                SpontSpikeRate{a}{b}(c) = mean(Psth{a}{b}(c,iipe:iopb)) ; %
            end

            SpontSpikeRate_mean(a,b) = mean(SpontSpikeRate{a}{b}) ;
            SpontSpikeRate_std(a,b) = std(SpontSpikeRate{a}{b}) ;
        else
            SpontSpikeRate_mean(a,b) = nan ;
            SpontSpikeRate_std(a,b) = nan ;                
        end
    end
end
    

% odor response potential
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                Pot{a}{b}(c) = mean(vData{a}{b}(c,iorb(a,b):iore(a,b))) ; %mV (odor response depol start: odor response depol end)
                deltaPot{a}{b}(c) = Pot{a}{b}(c)- Vrest{a}{b}(c) ; % mV
            end
            Pot_mean(a,b) = mean(Pot{a}{b}) ; %mV
            deltaPot_mean(a,b) = mean(deltaPot{a}{b}) ; % mV

            Pot_std(a,b) = std(Pot{a}{b}) ; %mV
            deltaPot_std(a,b) = std(deltaPot{a}{b}) ; % mV   
        else
            Pot_mean(a,b) = nan ; %mV
            deltaPot_mean(a,b) = nan ; % mV

            Pot_std(a,b) = nan ; %mV
            deltaPot_std(a,b) = nan ; % mV 
        end
    end
end


% odor response spike rate

for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                RespSpikeRate{a}{b}(c) = mean(Psth{a}{b}(c,iorb(a,b):iore(a,b))); %mV (odor response depol start: odor response depol end)
            end
              
            RespSpikeRate_mean(a,b) = mean(RespSpikeRate{a}{b}) ;
            RespSpikeRate_std(a,b) = std(RespSpikeRate{a}{b}) ;
        else
            RespSpikeRate_mean(a,b) = nan ;
            RespSpikeRate_std(a,b) = nan ;
        end       
    end
end

% delta spike rate
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            deltaSpikeRate{a}{b} = RespSpikeRate{a}{b} - SpontSpikeRate{a}{b} ;
            deltaSpikeRate_mean(a,b) = mean(deltaSpikeRate{a}{b}) ;
            deltaSpikeRate_std(a,b) = std(deltaSpikeRate{a}{b}) ;
        else
            deltaSpikeRate_mean(a,b) = nan ;
            deltaSpikeRate_std(a,b) = nan ;
        end
    end
end

% spont activity subtracted voltage and psth
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            for c = 1:NumTrials(a,b) ;
                vDataDelta{a}{b}(c,:) = vData{a}{b}(c,:) - Vrest{a}{b}(c) ;
                PsthDelta{a}{b}(c,:) = Psth{a}{b}(c,:) - SpontSpikeRate{a}{b}(c) ;
            end
        vDataDelta_mean{a}{b} = mean(vDataDelta{a}{b},1) ;
        PsthDelta_mean{a}{b} = mean(PsthDelta{a}{b},1) ;
        end
    end
end

% make background times second time stamps
if ~isempty(Input(A).(id3)) ;
    for a=1:length(Input(A).(id3)) ;
        bgTime(a) = (datenum(Input(A).(id3){a})- FirstTime)*24*60^2 ;
    end
else
    bgTime = [] ;
end

% assess input resisitance
if (iipe-driftCheckPnts)<iipb ; % if the points you want to check are before the begining of the pulse
    disp('I pulse is not that long') ;
else
    for a = 1:NumBackgrounds ; % for each background
        for b = 1:NumConcentrations ; % for concentration
            if NumTrials(a,b)>0 ;
                for c = 1:NumTrials(a,b) ;
                    Ipulse = mean(iData{a}{b}(c,iipe-driftCheckPnts:iipe)) - mean(iData{a}{b}(c,1:iipb)) ;
                    Vresp = mean(vData{a}{b}(c,iipe-driftCheckPnts:iipe)) - mean(vData{a}{b}(c,1:iipb)) ;
                    Rin{a}{b}(c) = Vresp/Ipulse ; % g ohms
                end
            end
        end
    end
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

% % input resistance
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(tDataN{a}{b},Rin{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%         end
%     end
% end
% ylabel('Rinput G Ohms')
% xlabel('trig time (sec)')
% 
% for a=1:length(bgTime) ;
%     plot([bgTime(a),bgTime(a)],get(gca,'ylim'),'k')
% end
    

% % potential as a function of concentration
% figure
% subplot(2,1,1)
% for a = 1:NumBackgrounds ;
%     errorbar(log10(Concentrations),Vrest_mean(a,:),Vrest_std(a,:),Vrest_std(a,:),'*:','Color',colorMat{NumConcentrations}(a,:))
%     hold on
%     errorbar(log10(Concentrations),Pot_mean(a,:),Pot_std(a,:),Pot_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
% end
% ylabel('potential')
% xlabel('log concentration')
% 
% subplot(2,1,2)
% for a = 1:NumBackgrounds ;
%     errorbar(log10(Concentrations),deltaPot_mean(a,:),deltaPot_std(a,:),deltaPot_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
%     hold on
% end
% ylabel('delta potential')
% xlabel('log concentration')
% 
% 
% %  spikes as a function of concentration
% figure
% subplot(2,1,1)
% for a = 1:NumBackgrounds ;
%     errorbar(log10(Concentrations),SpontSpikeRate_mean(a,:),SpontSpikeRate_std(a,:),SpontSpikeRate_std(a,:),'*:','Color',colorMat{NumConcentrations}(a,:))
%     hold on
%     errorbar(log10(Concentrations),RespSpikeRate_mean(a,:),RespSpikeRate_std(a,:),RespSpikeRate_std(a,:),'*-','Color',colorMat{NumConcentrations}(a,:))
% end
% ylabel('spike rate')
% xlabel('log concentration')
% 
% subplot(2,1,2)
% for a = 1:NumBackgrounds ;
%     errorbar(log10(Concentrations),deltaSpikeRate_mean(a,:),deltaSpikeRate_std(a,:),deltaSpikeRate_std(a,:),'+-','Color',colorMat{NumConcentrations}(a,:))
%     hold on
% end
% ylabel('delta spike rate')
% xlabel('log concentration')

% % mean voltage for different concentrations and backgrounds and time
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             subplot(NumConcentrations+1,2,b*2-1)
%             plot(time,vData_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             title(num2str(Concentrations(b)))
%             hold on
%             plot(time(1,[iorb(a,b),iore(a,b)]),vData_mean{a}{b}(1,[iorb(a,b),iore(a,b)]),'o','Color',colorMat{NumConcentrations}(a,:))
%             
%             subplot(NumConcentrations+1,2,b*2)
%             plot(time,vData_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             xlim([3,4])
%             hold on
%             plot(time(1,[iorb(a,b),iore(a,b)]),vData_mean{a}{b}(1,[iorb(a,b),iore(a,b)]),'o','Color',colorMat{NumConcentrations}(a,:))
%         end
%     end
% end
% 
% subplot(NumConcentrations+1,2,2*NumConcentrations+1:2*NumConcentrations+2)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(tDataN{a}{b},Vrest{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%             plot(tDataN{a}{b},Pot{a}{b},'o','Color',colorMat{b}(a,:))
%         end
%     end
% end
% ylabel('potential (mV)')
% xlabel('trig time (sec)')
% %legend('rest','evoked')
% 
% for a=1:length(bgTime) ;
%     plot([bgTime(a),bgTime(a)],get(gca,'ylim'),'k')
% end
% % 
% % 
% % 
% % psth and time stuff 
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             subplot(NumConcentrations+1,2,b*2-1)
%             plot(time,Psth_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             title(num2str(Concentrations(b)))
%             hold on
%             plot(time(1,[iorb(a,b),iore(a,b)]),Psth_mean{a}{b}(1,[iorb(a,b),iore(a,b)]),'o','Color',colorMat{NumConcentrations}(a,:))
%             
%             subplot(NumConcentrations+1,2,b*2)
%             plot(time,Psth_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             xlim([3,4])
%             hold on
%             plot(time(1,[iorb(a,b),iore(a,b)]),Psth_mean{a}{b}(1,[iorb(a,b),iore(a,b)]),'o','Color',colorMat{NumConcentrations}(a,:))
%         end
%     end
% end
% 
% subplot(NumConcentrations+1,2,2*NumConcentrations+1:2*NumConcentrations+2)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(tDataN{a}{b},SpontSpikeRate{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%             plot(tDataN{a}{b},RespSpikeRate{a}{b},'o','Color',colorMat{b}(a,:))
%         end
%     end
% end
% xlabel('trig time (sec)')
% ylabel('spike rate (hz)')
% %legend('spont','evoked')
%         
% for a=1:length(bgTime) ;
%     plot([bgTime(a),bgTime(a)],get(gca,'ylim'),'k')
% end

% mean voltage and psth for different concentrations and backgrounds and time
% figure 
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             subplot(NumConcentrations+2,2,b*2-1)
%             plot(time,vData_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             xlim([2.5,5])
%             ylabel(num2str(Concentrations(b)))
%             hold on
%             %plot(time(1,[iorb(a,b),iore(a,b)]),vData_mean{a}{b}(1,[iorb(a,b),iore(a,b)]),'o','Color',colorMat{NumConcentrations}(a,:))
%             
%             subplot(NumConcentrations+2,2,b*2)
%             plot(time,Psth_mean{a}{b},'Color',colorMat{NumConcentrations}(a,:))
%             xlim([2.5,5])
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
%             plot(tDataN{a}{b},Pot{a}{b},'o','Color',colorMat{b}(a,:))
%         end
%     end
% end
% ylabel('potential (mV)')
% 
% 
% for a=1:length(bgTime) ;
%     plot([bgTime(a),bgTime(a)],get(gca,'ylim'),'k')
% end
% 
% subplot(NumConcentrations+2,2,2*NumConcentrations+3:2*NumConcentrations+4)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(tDataN{a}{b},SpontSpikeRate{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%             plot(tDataN{a}{b},RespSpikeRate{a}{b},'o','Color',colorMat{b}(a,:))
%         end
%     end
% end
% xlabel('trig time (sec)')
% ylabel('spike rate (hz)')
%         
% for a=1:length(bgTime) ;
%     plot([bgTime(a),bgTime(a)],get(gca,'ylim'),'k')
% end


% comparing spikerate and potential (spont and evoked)
figure 
subplot(3,1,1)
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            plot(Vrest{a}{b},SpontSpikeRate{a}{b},'*','Color',colorMat{b}(a,:))
            hold on
            plot(Pot{a}{b},RespSpikeRate{a}{b},'o','Color',colorMat{b}(a,:))
        end
    end
end
xlabel('potential (mV)')
ylabel('spike rate (hz)')

% comparing spikerate and potential (change in spont to evoked)
subplot(3,1,2) 
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            plot(deltaPot{a}{b},deltaSpikeRate{a}{b},'+','Color',colorMat{b}(a,:))
            hold on
        end
    end
end
xlabel('delta potential (mV)')
ylabel('delta spike rate (hz)')

% comparing odor evoked spikerate and vrest
subplot(3,1,3)
for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration
        if NumTrials(a,b)>0 ;
            plot(Vrest{a}{b},deltaSpikeRate{a}{b},'+','Color',colorMat{b}(a,:))
            hold on
        end
    end
end
xlabel('rest potential (mV)')
ylabel('delta spike rate (hz)')


% % comparing Vrest to Vpot on each trial
% figure(3)
% subplot(2,1,1)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(Vrest{a}{b},deltaPot{a}{b},'*','Color',colorMat{b}(a,:))
%             hold on
%         end
%     end
% end
% xlabel('Vrest (mV)')
% ylabel('mean odor evoked delta pot (mV)')
% 
% % comparing spont spikes to odor evoked on each trial
% subplot(2,1,2)
% for a = 1:NumBackgrounds ; % for each background
%     for b = 1:NumConcentrations ; % for concentration
%         if NumTrials(a,b)>0 ;
%             plot(SpontSpikeRate{a}{b},deltaSpikeRate{a}{b},'o','Color',colorMat{b}(a,:))
%             hold on
%         end
%     end
% end
% xlabel('spontaneous (hz)')
% ylabel('odor evoked delta spike rate (hz)')
%  
 
% for igor

% identifier = ['LogConcentration','cell',num2str(A)] ;
% ForIgor.(identifier) = log10(Concentrations) ;
% 
% for a = 1:NumBackgrounds ; % for each background
%     identifier = ['DeltaSRmean',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = deltaSpikeRate_mean(a,:) ;
%     
%     identifier = ['DeltaSRstd',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = deltaSpikeRate_std(a,:) ;   
% end

for a = 1:NumBackgrounds ; % for each background
    for b = 1:NumConcentrations ; % for concentration 
        if NumTrials(a,b)>0 ;
            identifier = ['Psth','LogPulse',num2str(abs(log10(Concentrations(b)))),'Back',num2str(a),'cell',num2str(A)] ;
            ForIgor.(identifier) = Psth_mean{a}{b} ;

            identifier = ['PsthDelta','LogPulse',num2str(abs(log10(Concentrations(b)))),'Back',num2str(a),'cell',num2str(A)] ;
            ForIgor.(identifier) = PsthDelta_mean{a}{b} ;        

            identifier = ['vData','LogPulse',num2str(abs(log10(Concentrations(b)))),'Back',num2str(a),'cell',num2str(A)] ;
            ForIgor.(identifier) = vData_mean{a}{b} ;

            identifier = ['vDataDelta','LogPulse',num2str(abs(log10(Concentrations(b)))),'Back',num2str(a),'cell',num2str(A)] ;
            ForIgor.(identifier) = vDataDelta_mean{a}{b} ;   
        end
    end
end

identifier = ['time','cell',num2str(A)] ;
ForIgor.(identifier) = time;  




