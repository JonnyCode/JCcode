function ForIgor = Icomparison(Input,A) ; 

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);

epochsControl = str2num(Input(A).StepControl) ;
epochsExp = str2num(Input(A).StepExp) ;

for a = 1:length(epochsControl) ; % for each spike epoch
    [dataControl(a,:), error] = ITCReadEpoch(epochsControl(a), 0, fp) ;    % get data
    [prePnts, error] = ITCGetStmPrePts(epochsControl(a), 0, 0, fp) ; 
    dataControl(a,:) = dataControl(a,:) - mean(dataControl(a,1:prePnts)) ;
end

for a = 1:length(epochsExp) ; % for each spike epoch
    [dataExp(a,:), error] = ITCReadEpoch(epochsExp(a), 0, fp) ;    % get data
    [prePnts, error] = ITCGetStmPrePts(epochsExp(a), 0, 0, fp) ;
    dataExp(a,:) = dataExp(a,:) - mean(dataExp(a,1:prePnts)) ;
    
    [SI, error] = ITCGetSamplingInterval(epochsExp(a), fp); % get sampling interval
    SI = SI * 1e-6; % Sampling interval in sec
end

if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

time = [1:length(dataControl)]*SI ;

meandataControl = mean(dataControl) ;
meandataExp = mean(dataExp) ;

meandataControl_sum = sum(meandataControl) ;
meandataExp_sum = sum(meandataExp) ;

[maxd,iControl] = max(abs(meandataControl)) ;
meandataControl_peak = meandataControl(iControl) ;

[maxd,iExp] = max(abs(meandataExp)) ;
meandataExp_peak = meandataExp(iExp) ;

%figures
figure
plot(time,meandataExp,'r')
hold on
plot(time,meandataControl)
plot(time(iControl),meandataControl_peak,'ob')
plot(time(iExp),meandataExp_peak,'or')

% for igor
identifier = ['IcontrolPeak',num2str(A)] ;
ForIgor.(identifier) = meandataControl_peak ;

identifier = ['IexpPeak',num2str(A)] ;
ForIgor.(identifier) = meandataExp_peak ;

% individual example
identifier = ['IcontrolPeakTime',num2str(A)] ;
ForIgor.(identifier) = time(iControl) ;

identifier = ['IexpPeakTime',num2str(A)] ;
ForIgor.(identifier) = time(iExp) ;

identifier = ['IcontrolMean',num2str(A)] ;
ForIgor.(identifier) = meandataControl ;

identifier = ['IexpMean',num2str(A)] ;
ForIgor.(identifier) = meandataExp ;

identifier = ['Itime',num2str(A)] ;
ForIgor.(identifier) = time ;



end
