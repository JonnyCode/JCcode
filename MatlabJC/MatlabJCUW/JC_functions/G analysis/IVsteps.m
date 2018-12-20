function ForIgor = IVsteps(Input,Parameters,id,A) ;
 
% function will get IV plots for time points requested
 
epochs_str = Input(A).(id) ; 
 
 % get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/mouse/',Input(A).cellname]);

for set = 1:length(epochs_str) ; % for each potential
    epochNums{set} = str2num(epochs_str{set}) ;
    for trace = 1:length(epochNums{set}) ;
        [Data{set}(trace,:), error] = ITCReadEpoch(epochNums{set}(trace), 0, fp);
        [SI, error] = ITCGetSamplingInterval(epochNums{set}(trace), fp); % get sampling interval
        SI = SI * 1e-6; % Sampling interval in sec
        [prePnts, error] = ITCGetStmPrePts(epochNums{set}(trace), 0, 0, fp) ;
    end
    DataMean(set,:) = mean(Data{set},1) ; % mean of epochs at    
    DataMean(set,:) = DataMean(set,:) - mean(DataMean(set,1:prePnts)) ; % subtract of offset

    identifier = ['Idata',num2str(set),id,num2str(A)] ;
    ForIgor.(identifier) = DataMean(set,:) ;
end


time = [SI:SI:(SI*length(DataMean))] ;

figure
subplot(2,1,1)
plot(time,DataMean)
hold on
xlabel('time (sec)')
ylabel('pA')
  
% interactive feature to pick points for a IV plot
q1 = input('select point for IV plot (y or n)?') ; % check if they want the interactive part
while strcmp(q1,'y') % while they still want another IV plot
    display('click leaf trace you want to add key word tag and then hit keyboard')
    pause
    ax = gca ; % get current axis
    cp = get(ax,'CurrentPoint') ; % this gives last point clicked by mouse matrix
    cp = cp(1,1:2) ; % this gives point x,y of click in appropriate units
    
    Parameters.IVtimes = [Parameters.IVtimes,cp(1)] ; % times you want IV plot
    q1 = input('select another time (y or n)?') ; % check if they want to select another trace   
end

color = ['b','r','g','k','c','y','b','r','g','k','c','y','b','r','g','k','c','y','b','r','g','k','c','y'] ;

for a = 1:length(Parameters.IVtimes) ; % for each time point
    IVpnts(a) = find(time<Parameters.IVtimes(a),1,'last')
    subplot(2,1,1)
    plot(time(IVpnts(a)),DataMean(:,IVpnts(a)),[color(a),'o']) % show selected time points
    subplot(2,1,2)
    plot(str2num(Input(A).([id,'Pots'])),DataMean(:,IVpnts(a)),[color(a),'*-']) %IV plot
    hold on

    identifier = ['IatTime',num2str(a),id,num2str(A)] ;
    ForIgor.(identifier) = DataMean(:,IVpnts(a)) ;
end

% identifier = ['IVtimes',id,num2str(A)] ;
% ForIgor.(identifier) = Parameters.IVtimes ; 

identifier = ['Pots',id,num2str(A)] ;
ForIgor.(identifier) = str2num(Input(A).([id,'Pots'])) ;

identifier = ['time',id,num2str(A)] ;
ForIgor.(identifier) = time ;
