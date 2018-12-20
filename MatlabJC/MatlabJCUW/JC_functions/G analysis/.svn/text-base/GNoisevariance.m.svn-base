function ForIgor = GNoisevariance(Input,Parameters,id,A) ;

% figure in driving force for variance exc and inh are likely to contribute to total signal
if strcmp(id,'Exc') ;
    DF = (-60-0) ;
    DFfactor = (-50-0)^2 ; % assuming V=-50 and Exc Rev = 0, squarred because a signal times a factor will have a varaince =  var original signal times the factor^2 
    else if strcmp(id,'Inh') ;
            DF = (0 - -60) ;
            DFfactor = (-50--80)^2 ; % Inh Rev = -80 ;
        end
end


% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

% get spike data and sampling interval
epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each epoch  
    [predata(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get cell attached data
    data(a,:) = predata(a,Parameters.PrePnts:Parameters.PrePnts+Parameters.StmPnts)/DF - mean(predata(a,1:Parameters.PrePnts))/DF ; % subtract off pre pnt mean
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp) ;  % get the light stimulus (this will be ignored if not cell attached)
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp) ;
    SI(a) = SI(a) * 1e-6 ; % Sampling interval in sec
end

time = [SI(1):SI(a):SI(a)*length(data)] ; 

for a = 1:10 ; % for each individual stim
    dataMean(a,:) = mean(data(a:10:end,:)) ; % mean of epochs caused by the same stim
    dataVar(a,:) = var(data(a:10:end,:)) ; % var of epochs caused by the same stim 
    ccMeanVar(a,:) = xcov(dataMean(a,:),dataVar(a,:)) ; % cross correlation of mean and variance 
end

dataVarMean = mean(mean(dataVar,2)) ; % the mean of the mean variance
dataVarStd = std(mean(dataVar,2)) ; % the std of the mean variance

dataCvMean = mean(mean(sqrt(dataVar)./abs(dataMean),2)) ; % mean of the CV (std/mean)
dataCvStd = std(mean(sqrt(dataVar)./abs(dataMean),2)) ; % mean of the CV (std/mean)

ccMVmean = mean(ccMeanVar) ; % mean cross correlation
ccX = [SI(a):SI(a):SI(a)*length(ccMVmean)] - (SI(a)*length(ccMVmean)+SI(a))/2 ; % cross corr x values
        
EstIVarMean = mean(mean(dataVar*DFfactor,2)) ; % the mean of the mean variance that g would contribute to current
EstIVarStd =  std(mean(dataVar*DFfactor,2)) ; % the std of the mean variance

figure
subplot(2,2,1:2)
plot(time,dataMean(3,:))
hold on
plot(time,dataVar(3,:),'r')
plot(time,data(3:10:end,:),'k')
title([id,A])
legend('mean','var')

subplot(2,2,3)
plot(ccX,ccMVmean)

subplot(2,2,4)
text(.1,.9,['gVariance = ',num2str(dataVarMean),' +/- ',num2str(dataVarStd)])
text(.1,.8,['CV = ',num2str(dataCvMean),' +/- ',num2str(dataCvStd)])
text(.1,.7,['Estmated I Variance = ',num2str(EstIVarMean),' +/- ',num2str(EstIVarStd)])

% prep for export

ForIgor.Nada = 'Not recorded' ;

end




