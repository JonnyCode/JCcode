% faseb 08 fig 4

lowroddata{1} = {'052208Bc1', '[446:8:518]','[445:8:518]','[537:8:579]','[536:8:579]'} ;
lowroddata{2} = {'052208Bc3','[692:8:764]','[691:8:764]','[783:8:855]','[782:8:855]'} ;
lowroddata{3} = {'052208Bc2','[394:8:426]','[393:8:426]','[451:2:469]','[450:2:469]'} ;
lowroddata{4} = {'051308Bc2','[431:8:503]','[430:8:503]','[523:8:555]','[522:8:555]'} ;
lowroddata{5} = {'051308Bc1','[541:8:613]','[540:8:613]','[637:8:709]','[636:8:709]'} ;

roddata{1} = {'052208Bc1', '[634:8:706]','[633:8:706]','[764:8:836]','[763:8:836]'} ;
roddata{2} = {'060408Bc1','[664:8:696]','[663:8:696]','[802:8:834]','[801:8:834]'} ;
roddata{3} = {'050808Bc1','[289:8:361]','[288:8:361]','[401:8:473]','[400:8:473]'} ;
roddata{4} = {'052908Bc3','[457:8:513]','[456:8:513]','[533:8:605]','[532:8:605]'} ;
roddata{5} = {'042308Bc2','[363:372]','[413:422]','[515:524]','[505:514]'} ;
roddata{6} = {'052208Bc3','[873:8:945]','[872:8:945]','[960:8:1032]','[959:8:1032]'} ;
roddata{7} = {'052208Bc2','[484:8:556]','[483:8:556]','[575:8:647]','[574:8:647]'} ;
roddata{8} = {'042308Bc1','[261:270]','[251:260]','[378:387]','[368:377]'}


lowrodspikedata{1} = {'052208Bc1','[46:8:199]','[47:8:199]'} ;
lowrodspikedata{2} = {'052208Bc2','[11:8:164]','[12:8:164]'} ;
lowrodspikedata{3} = {'052208Bc3','[411:8:484]','[412:8:484]'} ;
lowrodspikedata{4} = {'051308Bc4','[66:8:219]','[67:8:219]'} ;
lowrodspikedata{5} = {'051308Bc3','[227:8:300]','[228:8:300]'} ;
lowrodspikedata{6} = {'051308Bc2','[48:8:161]','[49:8:161]'} ;
lowrodspikedata{7} = {'051308Bc1','[84:8:237]','[85:8:237]'} ;
lowrodspikedata{8} = {'052908Bc3','[26:8:179]','[27:8:179]'} ;


rodspikedata{1} = {'052208Bc1','[316:8:389,206:8:279]','[317:8:389,207:8:279]'} ;
rodspikedata{2} = {'052208Bc2','[181:8:334]','[182:8:334]'} ;
rodspikedata{3} = {'052208Bc3','[501:8:654]','[502:8:654]'} ;
rodspikedata{4} = {'051308Bc4','[250:8:379]','[251:8:379]'} ;
rodspikedata{5} = {'051308Bc3','[387:8:460]','[388:8:460]'} ;
rodspikedata{6} = {'051308Bc2','[168:8:241]','[169:8:241]'} ;
rodspikedata{7} = {'051308Bc1','[244:8:317]','[245:8:317]'} ;
rodspikedata{8} = {'052908Bc3','[186:8:339]','[187:8:339]'} ;
rodspikedata{9} = {'052908Bc1','[59:8:132]','[60:8:132]'} ;
rodspikedata{10} = {'060408Bc1','[493:8:606]','[494:8:606]'} ;

for a = 1:length(lowroddata) ;
   CellName_str = lowroddata{a}{1} ;
   epochs_str = lowroddata{a}(2:5) ;
   [fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',CellName_str]);


for set = 1:length(epochs_str) ;
    epochNums{set} = str2num(epochs_str{set}) ;
    for trace = 1:length(epochNums{set}) ;
        [Data{set}(trace,:), error] = ITCReadEpoch(epochNums{set}(trace), 0, fp);
    end
    Trace_mean{set} = mean(Data{set},1) ;
    Trace_mean{set} = Trace_mean{set} - mean(Trace_mean{set}(1:5000)) ;
    if set<3
        Trace_mean{set} = Trace_mean{set}/-61.6 ;
    else
        Trace_mean{set} = Trace_mean{set}/61.6 ;
    end
    
    figure(a)
    plot(Trace_mean{set},'b')
    hold on

    Gmean(a,set) = mean(Trace_mean{set}(5000:10000)) ;
    
    if a ==1
        example{set} = Trace_mean{set} ;
    end
end

% clear Data Trace_mean
end
    
ForIgor.examplelowrodexcdec = example{1} ;
ForIgor.examplelowrodexcinc = example{2} ;
ForIgor.examplelowrodinhdec = example{3} ;
ForIgor.examplelowrodinhinc = example{4} ;

ForIgor.lowrodIncG = Gmean(:,4) ;
ForIgor.lowrodDecG = Gmean(:,3) ;

ForIgor.lowrodIncGMean = mean(ForIgor.lowrodIncG) ;
ForIgor.lowrodDecGMean = mean(ForIgor.lowrodDecG) ;

ForIgor.examplerodexcdec = example{1} ;
ForIgor.examplerodexcinc = example{2} ;
ForIgor.examplerodinhdec = example{3} ;
ForIgor.examplerodinhinc = example{4} ;

ForIgor.rodIncG = Gmean(:,4) ;
ForIgor.rodDecG = Gmean(:,3) ;

ForIgor.rodIncGMean = mean(ForIgor.rodIncG) ;
ForIgor.rodDecGMean = mean(ForIgor.rodDecG) ;


for a = 1:length(lowrodspikedata)
    CellName_str = lowrodspikedata{a}{1} ;
    epochs_str = lowrodspikedata{a}(2:3) ;
    [fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',CellName_str]);


for set = 1:length(epochs_str) ;
    epochNums{set} = str2num(epochs_str{set}) ;
    for trace = 1:length(epochNums{set}) ;
        [Data(trace,:), error] = ITCReadEpoch(epochNums{set}(trace), 0, fp);
        Data(trace,:) = Data(trace,:) - mean(Data(trace,1:5000)) ;
    end
        SpikePnts = SpikeDetection(Data,10,10000) ;
        
        for b= 1:length(SpikePnts) ;
        PSTH(b,:) = zeros(1,length(Data)) ;                     % make a vector of zeros
        PSTH(b,SpikePnts{b}) = 1 ;                              % indicate a one when there was a spike
        PSTHlong(b,:) = smooth([zeros(1,250),PSTH(b,:)],500) ;   % smooth PSTH by 100 points (10ms) - the zero padding shifts the time of the average to the later times as appropriate
        PSTHsmooth(b,:) = PSTHlong(b,1:end-250) ;                % make the PSTH the correct length
        end
    PSTHMean{a}(set,:) = mean(PSTHsmooth) ;
    clear PSTH PSTHlong PSTHsmooth Data
end
    
end

for a = 1:length(PSTHMean)
    PSTHInc(a,:) = PSTHMean{a}(1,:) ;
    PSTHDec(a,:) = PSTHMean{a}(2,:) ;
end

ForIgor.lowrodIncPSTH1 = PSTHInc(1,:)*10000 ;
ForIgor.lowrodIncPSTHMean = mean(PSTHInc)*10000;

ForIgor.lowrodDeccPSTH1 = PSTHDec(1,:)*10000 ;
ForIgor.lowrodDecPSTHMean = mean(PSTHDec)*10000 ;

ForIgor.rodIncPSTH1 = PSTHInc(1,:)*10000 ;
ForIgor.rodIncPSTHMean = mean(PSTHInc)*10000 ;

ForIgor.rodDeccPSTH1 = PSTHDec(1,:)*10000 ;
ForIgor.rodDecPSTHMean = mean(PSTHDec)*10000 ;


