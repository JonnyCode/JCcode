function ForIgor = KoStepConcatAnalyzer(DataBlock, DB, Params)

% this function will analyze full field light step response from array data
% before and after cells are knocked out (KO).  

% The data will not be mapped via binary white noise but will be concatinated.  The times of drug wash in and out must be specified. 

% JC 12/3/15


% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
bl = 20 ; % number of gwgb pulses averaged in each block

DrugT = DataBlock(DB).FfPulseConcatDrugTimes ;

% load concatiniated data
dataRun = load_data(DataBlock(DB).FfPulseConcat) ;
dataRun = load_neurons(dataRun) ;

% load first set of concatinated data (assumes only first set of data may have been interupted)  
dataRunTemp = load_data(DataBlock(DB).FfPulse{1}) ;
dataRunTemp = load_neurons(dataRunTemp) ;

% find appropriate set of triggers
FirstSetExcess = rem(length(dataRunTemp.triggers),4) ;
trigs = [dataRun.triggers(4:length(dataRunTemp.triggers)-FirstSetExcess)',dataRun.triggers(length(dataRunTemp.triggers)+5:end)'] ;

blg = 4*bl ; % block group (total number of triggers per block)
num_bl = floor(length(trigs)/(blg)) ; % number of blocks 

for cl=1:length(dataRun.spikes) ; % for each cell
    for b = 1:num_bl ; % for each block of trials
        b_time(b) = trigs(b*blg-blg+1) ; % time of first data block trigger
        [psthTemp, binsTemp] = get_psth(dataRun.spikes{cl}, trigs(b*blg-blg+1:4:b*blg),'stop',12) ;
        psth{cl}(b,:) = psthTemp ;
    end 
    psth_mean(cl,:) = mean(psth{cl}) ; % mean psth
end
psth_time = binsTemp ;
 
BrkPnts = floor([0:size(psth_mean,2)/4:size(psth_mean,2)]) ; % points of stim change

for cl=1:length(dataRun.spikes) ; % for each cell
    [mx,mi] = max(psth_mean(cl,:)) ; % max of average psth 
    si = find(BrkPnts<mi,1,'last') ; % find stimulus that caused the max response

    [psth_peak(cl,:),mi] = max(psth{cl}(:,BrkPnts(si)+1:BrkPnts(si+1)),[],2) ;
    psth_rsp_mean(cl,:) = mean(psth{cl}(:,BrkPnts(si)+1:BrkPnts(si+1)),2) ;
    psth_rsp_trans(cl,:) = psth_peak(cl,:)./psth_rsp_mean(cl,:) ; % response transience (peak/mean)
    hp = psth_peak(cl,:)/2 ; % half peak
    for b = 1:num_bl ; % for each block of trials
        psth_duty(cl,b) = sum(psth{cl}(b,BrkPnts(si)+1:BrkPnts(si+1))>hp(b))/BrkPnts(2) ; % fraction of points above half peak
        psth_peak_time(cl,b) = psth_time(mi(b)) ; % (s) time post stim change of peak  
    end
end

% averages across cells
psth_peak_mean = mean(psth_peak,1) ;
psth_rsp_mean_mean = mean(psth_rsp_mean,1) ;
psth_rsp_trans_mean = mean(psth_rsp_trans,1) ;
psth_duty_mean = mean(psth_duty,1) ;
psth_peak_time_mean = mean(psth_peak_time,1) ;

psth_peak_std = std(psth_peak,0,1) ;
psth_rsp_mean_std = std(psth_rsp_mean,0,1) ;
psth_rsp_trans_std = std(psth_rsp_trans,0,1) ;
psth_duty_std = std(psth_duty,0,1) ;
psth_peak_time_std = std(psth_peak_time,0,1) ;


% figures
figure
for cl=1:length(dataRun.spikes) ; % for each cell
    clf
    
    subplot(4,1,1)
    plot(psth_time,psth_mean(cl,:),'c','linewidth',4)
    hold on
    %plot(psth_time,psth{cl})
    for b = 1:num_bl ; % for each block of trials
        if DrugT(1)<b_time(b) && b_time(b)<DrugT(2) ; 
            plot(psth_time,psth{cl}(b,:),'r','linewidth',2*(1.1-b/num_bl))
        elseif DrugT(1)>b_time(b) ;
            plot(psth_time,psth{cl}(b,:),'k','linewidth',2*(1.1-b/num_bl))
        elseif DrugT(2)<b_time(b)
            plot(psth_time,psth{cl}(b,:),'b','linewidth',2*(1.1-b/num_bl))
        end    
        hold on
    end
    ylabel('firing rate (hz)')
    xlabel('time (s)')

        
    subplot(4,1,2)
    [ax,h1,h2] = plotyy(b_time,psth_peak(cl,:),b_time,psth_rsp_mean(cl,:)) ;
    h1.Marker = '*' ; h1.LineStyle = 'none' ; h2.Marker = '*' ; h2.LineStyle = 'none' ;
    hold on
    plot(DrugT,[psth_peak(cl,1),psth_peak(cl,1)],'r')
    ylabel('peak rate (hz)')
    xlabel('block start time (s)')
    
    subplot(4,1,3)
    plot(b_time,psth_peak_time(cl,:),'*')
    hold on
    plot(DrugT,[psth_peak_time(cl,1),psth_peak_time(cl,1)],'r')
    ylabel('peak time (s)')
    xlabel('block start time (s)')
    
    subplot(4,1,4)
    [ax,h1,h2] = plotyy(b_time,psth_duty(cl,:),b_time,1./psth_rsp_trans(cl,:)) ;
    h1.Marker = '*' ; h1.LineStyle = 'none' ; h2.Marker = '*' ; h2.LineStyle = 'none' ;
    hold on
    plot(DrugT,[psth_duty(cl,1),psth_duty(cl,1)],'r')
    ylabel('duty (%)')
    xlabel('block start time (s)')
    
    pause
end   


figure
subplot(3,1,1)
[ax,h1,h2] = plotyy(b_time,psth_peak_mean,b_time,psth_rsp_mean_mean) ;
h1.Marker = '*' ; h1.LineStyle = 'none' ; h2.Marker = '*' ; h2.LineStyle = 'none' ;
hold on
plot(DrugT,[psth_peak_mean(1),psth_peak_mean(1)],'r')
ylabel('peak rate (hz)')
xlabel('block start time (s)')

subplot(3,1,2)
plot(b_time,psth_peak_time_mean,'*')
hold on
plot(DrugT,[psth_peak_time_mean(1),psth_peak_time_mean(1)],'r')
ylabel('peak time (s)')
xlabel('block start time (s)')

subplot(3,1,3)
[ax,h1,h2] = plotyy(b_time,psth_duty_mean,b_time,1./psth_rsp_trans_mean) ;
h1.Marker = '*' ; h1.LineStyle = 'none' ; h2.Marker = '*' ; h2.LineStyle = 'none' ;
hold on
plot(DrugT,[psth_duty_mean(1),psth_duty_mean(1)],'r')
ylabel('duty (%)')
xlabel('block start time (s)')




check = [] ;

