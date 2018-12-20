function ForIgor = DfAnalyzer(DataBlock, DB, Params)

% this function will analyze full field light flashe response from array data


% JC 4/22/15
% DB = 3 ;

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
interFlashIntVar = .005 ; % (sec) expected precision of trigger intervals
minNumFlashes = 5 ; % min number of flashes in a block 

% compare spatial and temporal receptive fields across data defined sets 

% load data
dataRun = load_data(DataBlock(DB).DfPath) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_params(dataRun) ;

% params
NumCells = length(dataRun.spikes) ;

% find sets of flash stimuli
PutativeSet_i = diff(dataRun.triggers)-DataBlock(DB).DfParams.interFlashInt < interFlashIntVar &...
    diff(dataRun.triggers)-DataBlock(DB).DfParams.interFlashInt > -interFlashIntVar ; % if time interval before next trigger is within range
st=1 ; 
trigger_set_i{1} = [] ;
for a=1:length(dataRun.triggers) ; % for each trigger
    trigger_set_i{st} = [trigger_set_i{st},a] ; % put it in a set   
    if a<length(dataRun.triggers) ; % if its not the last trigger
        if PutativeSet_i(a) ; % next trigger is has the right interval
            st = st ; % keep it in the set
        else
            st = st+1 ; % put it in a new set
            trigger_set_i{st} = [] ;
        end
    end
end

% get psth in each cell for each set of stimuli 
ds = 1 ; % data set
for ts=1:length(trigger_set_i) ; % for each trigger set  
    if length(trigger_set_i{ts})>minNumFlashes ; % if more than X flashes
        for c = 1:NumCells ; % for each cell
            [psthTemp,psthTime] = get_smooth_psth(dataRun.spikes{c}, dataRun.triggers(trigger_set_i{ts}), 'start', -2, 'stop', DataBlock(DB).DfParams.interFlashInt) ; 
            psth{ds}(c,:) = psthTemp ;
            
            preStimPnts = [1:find(psthTime==0)] ; % points before flash
            preStimMean = mean(psthTemp(preStimPnts)) ; % mean before flash
            preStimStd = std(psthTemp(preStimPnts)) ; % std before flash
            
            psth_delta{ds}(c,:) = psthTemp -psthTemp(find(psthTime==0)) ; % response change
            [v,i] = max(abs(psth_delta{ds}(c,:))) ; % largest response change
            psth_delta_peak(ds,c) = psth_delta{ds}(c,i) ;
            
            RspPnts = find(abs(psth_delta{ds}(c,preStimPnts(end)+1:end))> preStimStd*NumStdaboveNoise) ; % "response points"
            psth_absRspSum(ds,c) = sum(abs(psth_delta{ds}(c,preStimPnts(end)+RspPnts))) ; % sum of response points
        end
        psth_mean(ds,:) = mean(psth{ds}) ; % average across all cells
        psth_std(ds,:) = std(psth{ds}) ;
        
        ds = ds+1 ;
    end   
end

if size(psth_mean,1) ~= length(DataBlock(DB).DfParams.Ftime)  ~= length(DataBlock(DB).DfParams.NDF) ; % if the number of identified flash blocks is not as expected
    error('Params or data length mismatch')
end

% estimate the relative stim intensity assuming flashes are within linear response range
Irel = DataBlock(DB).DfParams.Ftime./10.^DataBlock(DB).DfParams.NDF ;
Irel_unique = unique(Irel) ;

% Organise by and average across cells with the same intensity levels 
psth_delta_peak_Irel = zeros(length(Irel_unique), NumCells) ;

for ds = 1:size(psth_mean,1) ;
    [v,i] = intersect(Irel_unique,Irel(ds)) ; % index of unique intensity
    NumSame = sum(Irel_unique(i)==Irel) ; % number of same stim
    psth_delta_peak_Irel(i,:) = psth_delta_peak_Irel(i,:)+ psth_delta_peak(ds,:)/NumSame ;
    
    psth_absRspSum_Irel(i,:) = psth__absRspSum_Irel(i,:)+ psth_absRspSum(ds,:)/NumSame ;
end
psth_delta_peak_Irel_mean = mean(psth_delta_peak_Irel,2) ;
psth_delta_peak_Irel_meanOn = mean(psth_delta_peak_Irel(:,sum(psth_delta_peak_Irel,1)>0),2) ;
psth_delta_peak_Irel_meanOff = mean(psth_delta_peak_Irel(:,sum(psth_delta_peak_Irel,1)<0),2) ;

psth_absRspSum_Irel_mean = mean(psth_absRspSum_Irel,2) ;

% histograms
HistBins = [-200:5:200] ;
for ds = 1:size(psth_mean,1) ;
    HistValues= hist(psth_delta_peak(ds,:),HistBins) ; 
    pdf(ds,:) = HistValues/sum(HistValues) ;
end
    
% figures
Ftimes = unique(DataBlock(DB).DfParams.Ftime) ;
NDFs = unique(DataBlock(DB).DfParams.NDF) ;

for ds = 1:size(psth_mean,1) ;
    plot_i(ds) = (find(NDFs==DataBlock(DB).DfParams.NDF(ds))-1)*length(Ftimes) + find(Ftimes==DataBlock(DB).DfParams.Ftime(ds)) ; % subplot i
end

figure
for ds = 1:size(psth_mean,1) ;
    subplot(length(NDFs),length(Ftimes),plot_i(ds))
    plot(psthTime,psth_mean(ds,:),'k','linewidth',3) ;
    hold on
    %plot(psthTime,psth{ds},'c') ;
    plot(psthTime,psth_mean(ds,:)+psth_std(ds,:),'b:') ;
    plot(psthTime,psth_mean(ds,:)-psth_std(ds,:),'b:') ;
    text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
end

figure
for ds = 1:size(psth_mean,1) ;
    subplot(length(NDFs),length(Ftimes),plot_i(ds))
    plot(HistBins,pdf(ds,:)) ;
    text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
end

figure
plot(log(Irel_unique),psth_delta_peak_Irel)
hold on
plot(log(Irel_unique),psth_delta_peak_Irel_mean,'k','linewidth',2)
plot(log(Irel_unique),psth_delta_peak_Irel_meanOn,'k:','linewidth',2)
plot(log(Irel_unique),psth_delta_peak_Irel_meanOff,'k:','linewidth',2)

figure
plot(log(Irel_unique),psth_absRspSum_Irel)
hold on
plot(log(Irel_unique),psth_absRspSum_Irel_mean,'k','linewidth',2)

% for igor
ForIgor = struct() ;

VecName = ['psthTime','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psthTime) ; % psth time vector

VecName = ['histBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,HistBins) ; % pdf bins

for ds = 1:size(psth_mean,1) ;
    VecName = ['psth','NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),'Fdur',num2str(DataBlock(DB).DfParams.Ftime(ds)),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor, VecName, psth_mean(ds,:)) ; % psth vector
    
    VecName = ['pdf','NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),'Fdur',num2str(DataBlock(DB).DfParams.Ftime(ds)),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor, VecName, pdf(ds,:)) ; % pdf vector
end

% avearge intensity respones curves
VecName = ['Irel','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Irel_unique) ; % relative intensity vector

VecName = ['peak_mean','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_mean) ; % relative intensity vector

VecName = ['peak_mean_On','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_meanOn) ; % relative intensity vector

VecName = ['peak_mean_Off','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_meanOff) ; % relative intensity vector

    




