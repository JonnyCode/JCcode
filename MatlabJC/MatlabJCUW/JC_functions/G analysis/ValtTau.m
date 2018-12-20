function ForIgor = ValtTau(Input,id,A) ;

% this function will analyze alternating voltage exp current trace during
% inh blockade and assess the time it takes to achieve the exc rev
% potential. JC 11/2/09

% get data

try
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);
catch 
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
end

epochs = str2num(Input(A).(id)) ;
round = 0 ;
for a = 1:3:length(epochs) ; % for each spike epoch
    round = round +1 ;

    [dataExc(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data

    [dataInh(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    % get data
    
    [dataAltV(round,:), error] = ITCReadEpoch(epochs(a+2), 0, fp) ;    % get data
    [SI(round), error] = ITCGetSamplingInterval(epochs(a+2), fp); % get sampling interval
    SI(round) = SI(round) * 1e-6; % Sampling interval in sec
    
    [lightCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 0,fp); 
    [voltageCommand(round,:), error] = ITCReadEpochStm(epochs(a+2), 1,fp); % get voltage command
end

lightCommand = lightCommand(:,1:length(dataExc)) ;
lightCommand(lightCommand<0)=0 ;

voltageCommand = voltageCommand(:,1:length(dataExc)) ;

if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

cyclepnts = 100 ; % number of sample points in a cycle, pnts between leaving hold1 and returning (Also gets rid of first cycle)    
FirstAltPnt = (cyclepnts/2)+1 ; % first sample point you want to plot after begining of step from alternation 
%FirstAltPnt = 30 ;
LastAltPnt = FirstAltPnt ;  % last "                                                                    "

[prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

if strcmp('AltDS',id)
    SpatialStimParams = hdf5load(['~/Data/mouse/',Input(A).cellname,'_spatial.h5']) ;
    frameRate = 60 ;
    StrucString = ['params_epoch',num2str(epochs(1))] ; 
    Struct = SpatialStimParams.(StrucString) ;
    prePnts = floor(Struct.spatial_prepts/(frameRate*SI(1))) ;
    postPnts = floor(Struct.spatial_postpts/(frameRate*SI(1))) ;
end

samplerate = 1/SI(1) ; % Hz at which data was collected
time = [SI(1):SI(1):SI(1)*length(dataExc(1,:))] ; % time vector in seconds

% reshape data for each voltage cycle
for a=1:size(dataAltV,1) ;
    cutfactor = rem(length(dataAltV(a,500*cyclepnts+1:prePnts)),cyclepnts) ;
    IdataPre_epoch{a} = reshape(dataAltV(a,500*cyclepnts+1:prePnts-cutfactor)',cyclepnts,numel(dataAltV(a,500*cyclepnts+1:prePnts-cutfactor))/cyclepnts)' ;
    meanIdataPre(a,:) = mean(IdataPre_epoch{a}) ;
    varIdataPre(a,:) = var(IdataPre_epoch{a}) ;

    cutfactor = rem(length(dataAltV(a,prePnts+1:end-postPnts)),cyclepnts) ;
    Idata_epoch{a} = reshape(dataAltV(a,prePnts+1:end-postPnts-cutfactor)',cyclepnts,numel(dataAltV(a,prePnts+1:end-postPnts-cutfactor))/cyclepnts)' ;
    meanIdata(a,:) = mean(Idata_epoch{a}) ;
    varIdata(a,:) = var(Idata_epoch{a}) ;
end

diffvarIdata = mean(varIdata) - mean(varIdataPre) ;

figure
subplot(3,2,1)
plot(mean(meanIdata),'k')
hold on
plot(mean(meanIdataPre),'r')
title('mean I')

subplot(3,2,2)
plot(mean(varIdata),'k')
hold on
plot(mean(varIdataPre),'r')
title('var I')

subplot(3,2,3:6)
plot([SI(1):SI(1):cyclepnts*SI(1)],diffvarIdata,'b')
title('diff variance I')
hold on

CapacitancePntsExc = 1 ; %input('number of points to avoid in first half') ;
CapacitancePntsInh = 1 ; %input('number of points to avoid in second half') ;

% find time it takes to move achieve a plateau where the average rate of 
% change is at least 90% the rate of change during the defined plateau time.

PlateauPnts = 5 ; % number of points that define the plateau

IPlateauExc = mean(diffvarIdata(cyclepnts/2-PlateauPnts:cyclepnts/2)) ; 
IPlateauInh = mean(diffvarIdata(cyclepnts-PlateauPnts:end)) ;

% fit exponetial
fitExc = nlinfit([CapacitancePntsExc*SI(1):SI(1):(cyclepnts/2)*SI(1)],diffvarIdata(CapacitancePntsExc:cyclepnts/2),'exponential',[IPlateauExc IPlateauInh-IPlateauExc .004]) ;
fitInh = nlinfit([CapacitancePntsInh*SI(1):SI(1):(cyclepnts/2)*SI(1)],diffvarIdata(cyclepnts/2+CapacitancePntsInh:end),'exponential',[IPlateauInh IPlateauExc-IPlateauInh .004]) ;

expExc = fitExc(1) + fitExc(2) .* exp(-[SI(1):SI(1):(cyclepnts/2)*SI(1)] ./ fitExc(3));
expInh = fitInh(1) + fitInh(2) .* exp(-[SI(1):SI(1):(cyclepnts/2)*SI(1)] ./ fitInh(3));

% figures

subplot(3,2,3:6)
plot([cyclepnts/2-10:cyclepnts/2]*SI(1),IPlateauExc,'g')
plot([cyclepnts-10:cyclepnts]*SI(1),IPlateauInh,'r')
plot(CapacitancePntsExc*SI(1),diffvarIdata(CapacitancePntsExc),'g*')
plot((CapacitancePntsInh+cyclepnts/2)*SI(1),diffvarIdata(CapacitancePntsInh+cyclepnts/2),'r*')
plot([SI(1):SI(1):(cyclepnts/2)*SI(1)],expExc,'g')
plot([(cyclepnts/2)*SI(1)+SI(1):SI(1):cyclepnts*SI(1)],expInh,'r')
text(.003,IPlateauExc,['ExcTau= ',num2str(fitExc(3))]) ;
text(.008,IPlateauInh,['InhTau= ',num2str(fitInh(3))]) ;

% 
% % for Igor

% identifier = ['lightCommand',id,num2str(A)] ;
% ForIgor.(identifier) = lightCommand(1,:) ;
% 
% identifier = ['voltCommand',id,num2str(A)] ;
% ForIgor.(identifier) = voltageCommand(1,:) ;
% 
% identifier = ['Ialt',id,num2str(A)] ;
% ForIgor.(identifier) = dataAltV(1,:) ;
% 
% identifier = ['time',id,num2str(A)] ;
% ForIgor.(identifier) = time ;
% 
% for a=1:size(Idata_epoch{1},1) ;
%     identifier = ['Iindcycle',num2str(a),id,num2str(A)] ;
%     ForIgor.(identifier) = Idata_epoch{1}(a,:) ;
% end
% 
% 
% identifier = ['varStim',id,num2str(A)] ;
% ForIgor.(identifier) = mean(varIdata) ;
% 
% identifier = ['varPre',id,num2str(A)] ;
% ForIgor.(identifier) = mean(varIdataPre) ;
% 
% identifier = ['vardiff',id,num2str(A)] ;
% ForIgor.(identifier) = diffvarIdata ;
% 
% identifier = ['cycletime',id,num2str(A)] ;
% ForIgor.(identifier) = [SI(1):SI(1):cyclepnts*SI(1)] ;
% 
% 
% identifier = ['ExponInh',id,num2str(A)] ;
% ForIgor.(identifier) = expInh ;
% 
% identifier = ['ExponTime',id,num2str(A)] ;
% ForIgor.(identifier) = [(cyclepnts/2)*SI(1)+SI(1):SI(1):cyclepnts*SI(1)] ;
% 

identifier = ['ValtTauExc',id,num2str(A)] ;
ForIgor.(identifier) = fitExc(3) ;

identifier = ['ValtTauInh',id,num2str(A)] ;
ForIgor.(identifier) = fitInh(3) ;

