function ForIgor = Light2gFilters(Input,Parameters,id,A) ;

% Linear Filter, Transfer function and Nonlinearity from light stimuli to conductance

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stm(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SR(a) = 1/(SI(a) * 1e-6); % Sampling rate in Hz
end

[prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

if Input(A).ITC18flag == 1 ; % check if correct
    SI = SI*1.25 ;
    SR = 1/(SI *1e-6) ;
end

% currents into coductances
if strcmp(id,'Exc')
    data = data/-61.1 ; % change pA into nS
else
    data = data/61.1 ;
end

% rectifier stim
negstmPnts = find(stm<0) ;  % find indicies which would be getting a negative stim voltage
stm(negstmPnts) = 0 ;       % make these points zero

% subtract mean during stimulus from stimulus and data
stm = stm - repmat(mean(stm(:,prePnts+Parameters.STAPnts:end-postPnts),2),1,size(stm,2)) ; 
data = data - repmat(mean(data(:,prePnts+Parameters.STAPnts:end-postPnts),2),1,size(data,2)) ;

%prep stim and response for lin filter
cutfactor = rem(length(stm)-prePnts-postPnts-Parameters.STAPnts,Parameters.ChopPnts) ; % how many pnts remain when we attempt to evenly divide the signal 
reshfactor = floor((length(stm)-prePnts-postPnts-Parameters.STAPnts)/Parameters.ChopPnts)*size(stm,1) ; % number of rows will we have when we divide up the signal further 

signal = stm(:,prePnts+Parameters.STAPnts+1:end-postPnts-cutfactor) ; %take only the time varying stimuli and leave off any possible remainder before reshaping rows 
response = data(:,prePnts+Parameters.STAPnts+1:end-postPnts-cutfactor) ;

signal = reshape(signal',Parameters.ChopPnts,reshfactor)' ; % cut stimuli in rows into more rows
response = reshape(response',Parameters.ChopPnts,reshfactor)' ;

% signal = stm(:,prePnts+Parameters.STAPnts:end-postPnts) ;
% response = data(:,prePnts+Parameters.STAPnts:end-postPnts) ;

[LinearFilter] = LinFilterFinder(signal,response, SR(1), 60) ; % lin filter function (signal, response, samplerate, freqcuttoff) ;
LinearFilterNorm = LinearFilter/max(abs(LinearFilter)) ; % normalize by max amplitude

[xTF, TF] = PowerSpectrumFinder(LinearFilter, SR(1)) ; % get transfer function
TFN = TF/sum(TF) ; % normalize the power spectrum by the total variance

NewPowerSpec = SmoothPowerSpectrum(TFN, xTF, 2,1) ; % smooth the normalized transfer functions
TFNS = NewPowerSpec.PowerSpec ;
xTFNS = NewPowerSpec.Freq ;

% signalPhase = mean(angle(fft(signal,[],2)))*(180/pi) ;
% responsePhase = mean(angle(fft(response,[],2)))*(180/pi) ;
PhaseTransfer = angle(fft(LinearFilter)) ;

% prep structure Igor export
identifier = ['LinF',id,num2str(A)] ;
ForIgor.(identifier) = LinearFilter ;

identifier = ['LinFN',id,num2str(A)] ;
ForIgor.(identifier) = LinearFilterNorm ;

identifier = ['TFun',id,num2str(A)] ;
ForIgor.(identifier) = TF ;

identifier = ['TFunN',id,num2str(A)] ;
ForIgor.(identifier) = TFN ;

identifier = ['xTFun',id,num2str(A)] ;
ForIgor.(identifier) = xTF ;

identifier = ['TFunNS',id,num2str(A)] ;
ForIgor.(identifier) = TFNS ;

identifier = ['xTFunNS',id,num2str(A)] ;
ForIgor.(identifier) = xTFNS ;

% plot linear filter
figure
plot([1:length(LinearFilter)],LinearFilter) 
title(['Linear filter',id,num2str(A)])

% plot transfer function
figure
plot(xTF,TFN)
hold on
plot(xTFNS,TFNS,'r')
title(['transfer function',id,num2str(A)])
a=gca ;
set(a,'xScale','log')
set(a,'yScale','log')

end % end function
