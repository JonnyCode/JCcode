% chapter3 thesis midget step response data 100% contrast step

cellname = '060809Ec2' ;
epochsExcStr = '[149:153]' ; 
epochsInhStr = '[184:188]' ;

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',cellname]);

epochsExc = str2num(epochsExcStr) ;
epochsInh = str2num(epochsInhStr) ;

for a=1:length(epochs)
    [dataExc(a,:), error] = ITCReadEpoch(epochsExc(a), 0, fp);
    [dataInh(a,:), error] = ITCReadEpoch(epochsInh(a), 0, fp);
end
   
[SI, error] = ITCGetSamplingInterval(epochsExc(1), fp); 
SI = SI*10^-6 ;
time = [1:length(dataExc)]*SI ;

Exc = dataExc/-61 ;
Inh = dataInh/61 ; 

Exc = Exc-min(Exc(:)) ;
Inh = Inh-min(Inh(:)) ;

figure
plot(time,mean(Exc),'b')
hold on
plot(time,mean(Inh),'r')

% for igor

ForIgor.midgetExcStep = mean(Exc) ;
ForIgor.midgetInhStep = mean(Inh) ;
ForIgor.time = time ;

cd ~/data_analysis/TempToIgor/
exportStructToHDF5(ForIgor,'MidgetGsteps.h5','/')