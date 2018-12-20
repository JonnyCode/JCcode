dp = '/Volumes/lab/Experiments/Array/Analysis/2017-01-16-0/data009-map_KR/data009-map_KR'
dataRun = load_data(dp) ;
dataRun = load_neurons(dataRun) ;

PopDist_all = PopDistFinder(dataRun.spikes, dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_all_bin20 = PopDistFinder(dataRun.spikes, dataRun.triggers([1:3:600]),'PsthBinTime', 0.02,'BinSearchNumber', 2) ;

PopDist_all_bin100 = PopDistFinder(dataRun.spikes, dataRun.triggers([1:3:600]),'PsthBinTime', 0.1,'BinSearchNumber', 2) ;

PopDist_all_bin1k = PopDistFinder(dataRun.spikes, dataRun.triggers([1:3:600]),'PsthBinTime', 1,'BinSearchNumber', 2) ;

PopDist_all_Search0 = PopDistFinder(dataRun.spikes, dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 0) ;

PopDist_all_Search4 = PopDistFinder(dataRun.spikes, dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 4) ;

figure
plot(2,mean(PopDist_all.AcrossTrials(:)),'k*')
hold on
plot(0,mean(PopDist_all_Search0.AcrossTrials(:)),'k*')
plot(4,mean(PopDist_all_Search4.AcrossTrials(:)),'k*')
xlabel('number bins search')
ylabel('pop distance noise')

figure
plot(10,mean(PopDist_all.AcrossStimMinusAcrossTrials(:)),'k*')
hold on
plot(20,mean(PopDist_all_bin20.AcrossStimMinusAcrossTrials(:)),'k*')
plot(100,mean(PopDist_all_bin100.AcrossStimMinusAcrossTrials(:)),'k*')
plot(1000,mean(PopDist_all_bin1k.AcrossStimMinusAcrossTrials(:)),'k*')
xlabel('size bin')
ylabel('pop distance sig-noise')

dp = '/Volumes/lab/Experiments/Array/Analysis/2017-01-16-0/data006_KR/data006_KR'
dataRunBw = load_data(dp) ;
dataRunBw = load_neurons(dataRunBw) ;
dataRunBw = load_params(dataRunBw,'cell_type_depth', 5) ;

OffBiphasic_ids = get_cell_ids(dataRunBw,'off biphasic') ;
OffLessBiphasic_ids = get_cell_ids(dataRunBw,'off less biphasic') ;
OnBiphasic_ids = get_cell_ids(dataRunBw,'on biphasic') ;
OnLargeBiphasic_ids = get_cell_ids(dataRunBw,'on large biphasic') ;
DsOnOffUp_ids = get_cell_ids(dataRunBw,'ds onoff Up') ;
DsOnOffDown_ids = get_cell_ids(dataRunBw,'ds onoff Down') ;
DsOnOffLeft_ids = get_cell_ids(dataRunBw,'ds onoff Left') ;
DsOnOffRight_ids = get_cell_ids(dataRunBw,'ds onoff Right') ;
DsOnUp_ids = get_cell_ids(dataRunBw,'ds on Up') ;
DsOnRight_ids = get_cell_ids(dataRunBw,'ds on Right') ;
DsOnLeft_ids = get_cell_ids(dataRunBw,'ds on Left') ;

OffBiphasic_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,OffBiphasic_ids)) ;
OffLessBiphasic_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,OffLessBiphasic_ids)) ;
OnBiphasic_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,OnBiphasic_ids)) ;
OnLargeBiphasic_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,OnLargeBiphasic_ids)) ;
DsOnOffUp_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,DsOnOffUp_ids)) ;
DsOnOffDown_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,DsOnOffDown_ids)) ;
DsOnOffLeft_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,DsOnOffLeft_ids)) ;
DsOnOffRight_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,DsOnOffRight_ids)) ;
DsOnUp_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,DsOnUp_ids)) ;
DsOnRight_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,DsOnRight_ids)) ;
DsOnLeft_i = get_cell_indices(dataRun,intersect(dataRun.cell_ids,DsOnLeft_ids)) ;

PopDist_OffBiphasic = PopDistFinder(dataRun.spikes(OffBiphasic_i), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_OnBiphasic = PopDistFinder(dataRun.spikes(OnBiphasic_i), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_OffLessBiphasic = PopDistFinder(dataRun.spikes(OffLessBiphasic_i), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_OnLargeBiphasic = PopDistFinder(dataRun.spikes(OnLargeBiphasic_i), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_WithoutOffBiphasic = PopDistFinder(dataRun.spikes(setdiff([1:length(dataRun.spikes)],OffBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_WithoutOnBiphasic = PopDistFinder(dataRun.spikes(setdiff([1:length(dataRun.spikes)],OnBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_WithoutOffBiphasic_WithoutOnBiphasic = PopDistFinder(dataRun.spikes(setdiff([1:length(dataRun.spikes)],[OffBiphasic_i,OnBiphasic_i])), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_AllDs = PopDistFinder(dataRun.spikes([DsOnOffUp_i,DsOnOffDown_i,DsOnOffLeft_i,DsOnOffRight_i]), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

figure
plot(mean(PopDist_all.AcrossStim,1)-mean(PopDist_all.AcrossTrials,1),'k')
hold on
% plot(mean(PopDist_OffBiphasic.AcrossStim,1)-mean(PopDist_OffBiphasic.AcrossTrials,1),'g')
%plot(mean(PopDist_OnBiphasic.AcrossStim,1)-mean(PopDist_OnBiphasic.AcrossTrials,1),'g')
%plot(mean(PopDist_WithoutOnBiphasic.AcrossStim,1)-mean(PopDist_WithoutOnBiphasic.AcrossTrials,1),'r')
% plot(mean(PopDist_OffLessBiphasic.AcrossStim,1)-mean(PopDist_OffLessBiphasic.AcrossTrials,1),'g')
% plot(mean(PopDist_OnLargeBiphasic.AcrossStim,1)-mean(PopDist_OnLargeBiphasic.AcrossTrials,1),'c')
% plot(mean(PopDist_SmallSet1.AcrossStim,1)-mean(PopDist_SmallSet1.AcrossTrials,1),'y')
% plot(mean(PopDist_WithoutOffBiphasic.AcrossStim,1)-mean(PopDist_WithoutOffBiphasic.AcrossTrials,1),'r')
% plot(mean(PopDist_WithoutOffBiphasic_WithoutOnBiphasic.AcrossStim,1)-mean(PopDist_WithoutOffBiphasic_WithoutOnBiphasic.AcrossTrials,1),'b--')
plot(mean(PopDist_AllDs.AcrossStim,1)-mean(PopDist_AllDs.AcrossTrials,1),'g')


SmallSet1_ids =[1158,1502,1531,4009,1577,1637,1338,1697,1397,1756,1606,2626,...
    1579,1699,1471,1863,2866,1263,1636,1684,2628,1578,1698,2509] ;

SmallSet1_i = get_cell_indices(dataRun,intersect(SmallSet1_ids,dataRun.cell_ids)) ;

PopDist_SmallSet1 = PopDistFinder(dataRun.spikes(SmallSet1_i), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_WithoutOffBiphasic = PopDistFinder(dataRun.spikes(setdiff(SmallSet1_i,OffBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_WithoutOnBiphasic = PopDistFinder(dataRun.spikes(setdiff(SmallSet1_i,OnBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_WithoutOffBiphasic_WithoutOnBiphasic = PopDistFinder(dataRun.spikes(setdiff(SmallSet1_i,[OffBiphasic_i,OnBiphasic_i])), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_WithoutOffLessBiphasic = PopDistFinder(dataRun.spikes(setdiff(SmallSet1_i,OffLessBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_WithoutOnLargeBiphasic = PopDistFinder(dataRun.spikes(setdiff(SmallSet1_i,OnLargeBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_OffBiphasic = PopDistFinder(dataRun.spikes(intersect(SmallSet1_i,OffBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_OnBiphasic = PopDistFinder(dataRun.spikes(intersect(SmallSet1_i,OnBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_OffBiphasic_OnBiphasic = PopDistFinder(dataRun.spikes(intersect(SmallSet1_i,[OffBiphasic_i,OnBiphasic_i])), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_OffLessBiphasic = PopDistFinder(dataRun.spikes(intersect(SmallSet1_i,OffLessBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

PopDist_SmallSet1_OnLargeBiphasic = PopDistFinder(dataRun.spikes(intersect(SmallSet1_i,OnLargeBiphasic_i)), dataRun.triggers([1:3:600]),'PsthBinTime', 0.01,'BinSearchNumber', 2) ;

figure
plot(mean(PopDist_SmallSet1.AcrossStimMinusAcrossTrials,1),'k')
hold on
%plot(mean(PopDist_SmallSet1_WithoutOffBiphasic.AcrossStimMinusAcrossTrials,1),'r')
%plot(mean(PopDist_SmallSet1_WithoutOnBiphasic.AcrossStimMinusAcrossTrials,1),'r')
%plot(mean(PopDist_SmallSet1_WithoutOffBiphasic_WithoutOnBiphasic.AcrossStimMinusAcrossTrials,1),'c--')
%plot(mean(PopDist_SmallSet1_WithoutOnLargeBiphasic.AcrossStimMinusAcrossTrials,1),'b--')
plot(mean(PopDist_SmallSet1_WithoutOffLessBiphasic.AcrossStimMinusAcrossTrials,1),'r')
%plot(mean(PopDist_SmallSet1_OffBiphasic.AcrossStimMinusAcrossTrials,1),'g')
%plot(mean(PopDist_SmallSet1_OnBiphasic.AcrossStimMinusAcrossTrials,1),'g')
%plot(mean(PopDist_SmallSet1_OffBiphasic_OnBiphasic.AcrossStimMinusAcrossTrials,1),'c')
%plot(mean(PopDist_SmallSet1_OnLargeBiphasic.AcrossStimMinusAcrossTrials,1),'b')
plot(mean(PopDist_SmallSet1_OffLessBiphasic.AcrossStimMinusAcrossTrials,1),'g')

figure
plot(mean(PopDist_SmallSet1.r_mean(1:end-1),1),mean(PopDist_SmallSet1.AcrossStimMinusAcrossTrials,1),'k')
hold on
plot(mean(PopDist_SmallSet1_OffBiphasic.r_mean(1:end-1),1),mean(PopDist_SmallSet1_OffBiphasic.AcrossStimMinusAcrossTrials,1),'r')
plot(mean(PopDist_SmallSet1_OnBiphasic.r_mean(1:end-1),1),mean(PopDist_SmallSet1_OnBiphasic.AcrossStimMinusAcrossTrials,1),'g')
plot(mean(PopDist_SmallSet1_OnLargeBiphasic.r_mean(1:end-1),1),mean(PopDist_SmallSet1_OnLargeBiphasic.AcrossStimMinusAcrossTrials,1),'b')
plot(mean(PopDist_SmallSet1_OffLessBiphasic.r_mean(1:end-1),1),mean(PopDist_SmallSet1_OffLessBiphasic.AcrossStimMinusAcrossTrials,1),'y')

figure
plot(mean(PopDist_all.AcrossStim,1)-mean(PopDist_all.AcrossTrials,1),'k')
hold on
plot(mean(PopDist_OffBiphasic.AcrossStim,1)-mean(PopDist_OffBiphasic.AcrossTrials,1),'b')
plot(mean(PopDist_OnBiphasic.AcrossStim,1)-mean(PopDist_OnBiphasic.AcrossTrials,1),'r')
plot(mean(PopDist_OffLessBiphasic.AcrossStim,1)-mean(PopDist_OffLessBiphasic.AcrossTrials,1),'g')
plot(mean(PopDist_OnLargeBiphasic.AcrossStim,1)-mean(PopDist_OnLargeBiphasic.AcrossTrials,1),'c')
plot(mean(PopDist_SmallSet1.AcrossStim,1)-mean(PopDist_SmallSet1.AcrossTrials,1),'y')
plot(mean(PopDist_WithoutOffBiphasic.AcrossStim,1)-mean(PopDist_WithoutOffBiphasic.AcrossTrials,1),'b--')
plot(mean(PopDist_WithoutOffBiphasic_WithoutOnBiphasic.AcrossStim,1)-mean(PopDist_WithoutOffBiphasic_WithoutOnBiphasic.AcrossTrials,1),'b--')

figure
plot(mean(PopDist_all.AcrossStimOverAcrossTrials,1),'k')
hold on
plot(mean(PopDist_OffBiphasic.AcrossStimOverAcrossTrials,1),'b')
plot(mean(PopDist_OnBiphasic.AcrossStimOverAcrossTrials,1),'r')
plot(mean(PopDist_OffLessBiphasic.AcrossStimOverAcrossTrials,1),'g')
plot(mean(PopDist_OnLargeBiphasic.AcrossStimOverAcrossTrials,1),'c')
plot(mean(PopDist_SmallSet1.AcrossStimOverAcrossTrials,1),'y')

figure
plot(mean(PopDist_all.r_mean,1),'k')
hold on
plot(mean(PopDist_OffBiphasic.r_mean,1),'b')
plot(mean(PopDist_OnBiphasic.r_mean,1),'r')
plot(mean(PopDist_OffLessBiphasic.r_mean,1),'g')
plot(mean(PopDist_OnLargeBiphasic.r_mean,1),'c')
plot(mean(PopDist_SmallSet1.r_mean,1),'y')
