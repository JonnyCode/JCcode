%% Run this script after RodModelProcedures
% Here, for now here I'll visually look at the data

% To first find the Impulse response and it's spectrum, I have to fill out
% the Pulse response with zeros to match the length of the Noise stimuli
clear steps dataMat dataVect Offset handle cond view labelList i j k scaleFactor...
    singleConversion errMatr sqErrMatr modelSqErrMatr modelErrMatr varMatr adjustH bin
if input('Clear all variables?...')
    clear  s h handle sindex hindex cond view tdMeanSqErrVect modelTdMeanSqErrVect truncatedEpochData...
        EpochData truncatedFilterData FilterData H S Y h r s y ans shortH avgOtherImpulseResponse impulseVarianceMatr impulseVarianceVect ...
        avgOtherResponse meanSqErrAcrossTrialsVsModel meanSqErrAcrossTrials meanSqErrAcrossTrialsYCollAreas varVect...
        stimulus poissNoise rndS index IgorStructure collectingArea impulseScalingFactors impulseScalingFactorIndex...
        hMatrix HMatrix YMatrix yMatrix avgOtherPoissNoise extremePoissNoise rBins yBins 
end

clear impulseScalingFactors;
collectingArea = .5;

clear adjustH max min stepsize

view = 1;

labelList = '--';
for cond = 1:length(CellInfo.EpochCondition)
    labelList = strcat(labelList,sprintf('\n%g: %s',cond,CellInfo.EpochCondition(cond).Label{1}));
end
labelList = strcat(labelList,sprintf('\n--'));
disp(labelList);
sindex = input('Enter index no. of stimulus...');
hindex = input('Enter index no. of impulse response...');
disp(sprintf('Convolving stimulus (%s) with impulse (%s)...',...
    CellInfo.EpochCondition(sindex).Label{1},CellInfo.EpochCondition(hindex).Label{1}));

%% Stimulus and scaling
% Stimulus initially given in volts with an offset.
% Rh*/time bin
stimulus = EpochCondition(sindex).StimulusVector...
    (FindSearchPara(EpochCondition(sindex),'PrePoints')+1:length(EpochCondition(sindex).StimulusVector));

% <Rh*/time point> = Stimulus [V] * PhotonFlux [photon/um^2/sec/V] / 
%  sampling interval [points/sec] * Collecting Area [Rh*]/[photon/um^2]
scaleFactor = CellInfo.PhotonFlux / FindSearchPara(EpochCondition(sindex),'SampInterv') * collectingArea;
stimulusMean = EpochCondition(sindex).StimulusVector(1);
s = stimulus - stimulusMean; % subtract mean
stimulusMean = stimulusMean * scaleFactor;
stimulus = stimulus * scaleFactor; %scale the stimulus while maintaining mean
s = s * scaleFactor; %scale the stimulus


%% Impulse response scaling
% Impulse response initially given in pA
% singleResponse (pA / Rh*) = 
%         ImpulseResp [pA] / (integratedFlux [photons/um^2] *
%         CollectingArea [Rh*]/[photon/um^2])
% Collecting Area could be a source of discrepancy (factor of .75-1.5)
h = EpochCondition(hindex).AverageResponse...
    (FindSearchPara(EpochCondition(hindex),'PrePoints')+1:length(EpochCondition(hindex).StimulusVector));

if exist('impulseScalingFactors','var')
    hMatrix = zeros(length(impulseScalingFactors),length(s));
    for i = 1:length(impulseScalingFactors)
        singleConversion = EpochCondition(hindex).UserInfo.StimulusAmp * collectingArea / impulseScalingFactors(i);
        hMatrix(i,1:length(h))= h/singleConversion;
    end
end

singleConversion = EpochCondition(hindex).UserInfo.StimulusAmp * collectingArea; % for pulse on green LED: ~5 Rh*
h = h/singleConversion; % impulse scaled to pA per Rh*
if length(h)<length(s)
    h = [h,zeros(1,length(s)-length(h))];
end


%% Model: Linear filter without noise!

r = EpochCondition(sindex).AverageResponse...
    (FindSearchPara(EpochCondition(sindex),'PrePoints')+1:length(EpochCondition(sindex).AverageResponse));
S = fft(s);
H = fft(h);
Y = H.*S;
y = ifft(Y);


if exist('impulseScalingFactors','var')
    HMatrix = hMatrix;
    YMatrix = hMatrix;
    yMatrix = hMatrix;

    for i = 1:length(impulseScalingFactors)
        HMatrix(i,:) = fft(hMatrix(i,:));
        YMatrix(i,:) = HMatrix(i,:).*S;
        yMatrix(i,:) = ifft(YMatrix(i,:));
    end 
end


%% now comparisons can be made between r (average response) and y (the
% predicted response from the model)

% Error and correlation function

% Calculate mean squared error of each response compared to mean response
% without that response

%% Compare the values of the model to the values of the response.
% In the next cell, analyse the effects of introducing a scaling factor,
% calculated in this cell

maxAmp = max(r);
minAmp = min(r);
steps = 100;
bins = [minAmp:(maxAmp-minAmp)/steps:maxAmp];
yBins = zeros(1,length(bins)-1);
rBins = yBins;
if exist('impulseScalingFactors','var')
    yBinsMatrix = zeros(length(impulseScalingFactors),length(bins)-1);
end
for bin = 1:length(bins)-1
    indices = (r>=bins(bin))==(r<bins(bin+1));
    rBins(bin) = mean(r(indices));
    yBins(bin) = mean(y(indices));
    if exist('impulseScalingFactors','var')
        for impulseScalingFactor = 1:length(impulseScalingFactors)
            yBinsMatrix(impulseScalingFactor,bin) = mean(yMatrix(impulseScalingFactor,indices));
        end
    end
end

clear maxAmp minAmp step bin sat pulse a satindex impulseScalingFactor


%% Calculate the average of impulse response, withholding ith

if ~exist('truncatedFilterData')
    [FilterData,Offset] = GetEpochData(CellInfo,EpochCondition(hindex));
    for i = 1:size(FilterData,1)
        truncatedFilterData(i,:) = FilterData(i,FindSearchPara(EpochCondition(hindex),'PrePoints')+1:end)/singleConversion;
    end
    avgOtherImpulseResponse = truncatedFilterData;
    for i = 1:size(truncatedFilterData,1)
        temp = zeros(1,(size(truncatedFilterData,2)));
        for j = 1:size(truncatedFilterData,1)
            if i~=j
                temp = temp+truncatedFilterData(j,:);
            end
        end
        avgOtherImpulseResponse(i,:) = temp/(j-1);
    end
    clear temp;
end


%% <(impulse response - average of others)^2>, variance in the impulse response

shortH = h(1:size(truncatedFilterData,2));
impulseVarianceMatr = truncatedFilterData - repmat(shortH,size(truncatedFilterData,1),1);
sqErrMatr = (truncatedFilterData - avgOtherImpulseResponse).^2;
impulseVarianceMatr = impulseVarianceMatr.^2;
tdImpulseMeanSqErr  = mean(sqErrMatr, 1); % essentially the variance but comparing to the other averages
tdImpulseVarVect = mean(impulseVarianceMatr,1);


%% Average Filter
% r(t) = F(t)*s(t) = integral(F(tau)s(t-tau)dtau,-inf,inf)
% <r(t)> = <F(t)*s(t)> = <integral(F(tau)s(t-tau)dt,-inf,inf)>
% operate on both sides, inside the exp. value signs by integral(s(t-tau')dt)
% <integral(r(t)s(t-tau')dt)> =
% <integral(F(tau)dtau)*integral(dtau'*s(t-tau)s(t-tau'))
% From Wolfram math: corr(f,g) = conv(conj(f(-t)),g) So 
% <R * S*> = <F> * S(omega) 
% thus:
randomSeed = 0;
randomSeed = input(sprintf('Evaluate Average Filter from random seed noise data? 0 - no, 1 - yes...'));
if (randomSeed)
    R = fft(truncatedEpochData,[],2);
    Smat = R;
    epochNumbers = EpochCondition(sindex).EpochNumbers;
    repeatEpochNumbers = EpochCondition(4).EpochNumbers;
    trimmedEpochNumbers = [];
    for i = 1:length(epochNumbers)
        for j = 1:length(repeatEpochNumbers)
            if repeatEpochNumbers(j) == epochNumbers(i)
                break
            end
        end
        trimmedEpochNumbers(end+1) = epochNumbers(i);
    end
    for epoch = 1:size(R,1)
        [fp, error] = ITCInitializeAnalysis(500000, file);
        [stm, error] = ITCReadEpochStm(EpochNumbers(epoch), 0, fp);
        Smat(i,:) = stm((length(stm)+1-length(r)):end);
    end

    powerSpectrum = dspdata.msspectrum(s);
    F = mean(R.*conj(Smat))./powerSpectrum.data;
    f = ifft(F);
    if(view)
        close(handle)
    end
    clear dataMat dataVect Offset handle cond view labelList i j scaleFactor singleConversion errMatr sqErrMatr modelSqErrMatr modelErrMatr varMatr
    return
end
clear randomSeed 

%% Mean squared error of stimulus response to average, stimulus response to model
% Use this to test the optimum fit of scaling factors

if ~exist('truncatedEpochData','var')
    [EpochData, Offset] = GetEpochData(CellInfo, CellInfo.EpochCondition(sindex));
    for i = 1:size(EpochData,1)
        truncatedEpochData(i,:) = EpochData(i,(size(EpochData,2)-length(s)+1):end);
    end
end

if ~exist('avgOtherResponse','var')
    avgOtherResponse = truncatedEpochData;
    for i = 1:size(truncatedEpochData,1)
        temp = zeros(1,(size(truncatedEpochData,2)));
        for j = 1:size(truncatedEpochData,1)
            if i~=j
                temp = temp+truncatedEpochData(j,:);
            end
        end
        avgOtherResponse(i,:) = temp/(j-1);
    end
    clear temp
end

% Normal calculation of the mean squared error across trials.
if ~exist('meanSqErrAcrossTrials','var')

    meanSqErrAcrossTrials = mean((avgOtherResponse - truncatedEpochData).^2,2);
    meanSqErrAcrossTrialsVsModel = mean((repmat(y,size(truncatedEpochData,1),1) - truncatedEpochData).^2,2);

end

% Now for the scaling data, minimize the error with respect to the average
% if exist('impulseScalingFactors','var')
%     
%     temp = mean((yMatrix-repmat(r,size(yMatrix,1),1)).^2,2); 
%     impulseScalingFactorIndex = find(temp == min(temp));
%     impulseScalingFactor = impulseScalingFactors(impulseScalingFactorIndex);
%     meanSqErrAcrossTrialsVsModelOptimum = meanSqErrAcrossTrialsVsModelSF(:,impulseScalingFactorIndex);
%     yOptimum = yMatrix(impulseScalingFactorIndex,:);
%     
% end

%% <(stim response - average of others)^2>, <(stim response - model)^2>

sqErrMatr = (truncatedEpochData - avgOtherResponse).^2;
modelSqErrMatr = (truncatedEpochData - repmat(y,size(truncatedEpochData,1),1)).^2;
varMatr = (truncatedEpochData - repmat(r,size(truncatedEpochData,1),1)).^2;

modelTdMeanSqErr = mean(modelSqErrMatr,1);
tdMeanSqErr  = mean(sqErrMatr, 1); 
varVect = mean(varMatr,1);

%% Model using poisson statistics
% for this I'm going to need to put the mean back into the stimulus.

% Calculate the fourier transform of a constant stimulus of the mean light
% intensity

stimulusMeanVect = zeros(1,size(truncatedEpochData,2));
stimulusMeanVect(:) = stimulusMean;
temp = fft(stimulusMeanVect).*H;
meanResponse = ifft(temp);
%temp = fft(stimulusMeanVect).*H*impulseScalingFactor;
%meanResponseOptimum = ifft(temp);
clear stimulusMeanVect stimulusMean

% Now run a bunch of random trials of poisson trials and calculate the
% average response

rndS = zeros(size(truncatedEpochData));
poissNoise = rndS;
%poissNoiseOptimum = poissNoise;
for i = 1:size(rndS,1)
    temp = poissrnd(stimulus);
    index = find(temp>0);
    rndS(i,index) = temp(index);
    RNDS = fft(rndS(i,:));
    temp = H.*RNDS;
    poissNoise(i,:) = ifft(temp);
%    temp = H.*RNDS*impulseScalingFactor;
%    poissNoiseOptimum(i,:) = ifft(temp);
end

if ~exist('extremePoissNoise','var')
    disp('Here we go...');
    extremePoissNoise = zeros(2000,size(truncatedEpochData,2));
    for i = 1:size(extremePoissNoise,1)
        temp = poissrnd(stimulus);
        index = find(temp>0);
        extremePoissNoise(i,index) = temp(index);
        RNDS = fft(extremePoissNoise(i,:));
        temp = H.*RNDS;
        extremePoissNoise(i,:) = ifft(temp);
    end
    
    p = mean(extremePoissNoise)-meanResponse;
    poissTdVariance = var(extremePoissNoise);
    disp('ok, finished');
end


poissNoise = poissNoise - repmat(meanResponse,size(poissNoise,1),1);



% poissNoiseOptimum = poissNoiseOptimum - repmat(meanResponseOptimum,size(poissNoiseOptimum,1),1);
% pOptimum = mean(poissNoiseOptimum);

%% Calculate the timedependent mean squared error of the poisson noise, see how it compares

% withhold ith pData vector
if ~exist('avgOtherPoissNoise','var')
    avgOtherPoissNoise = poissNoise;
    for i = 1:size(avgOtherPoissNoise,1)
        temp = zeros(1,size(avgOtherPoissNoise,2));
        for j = 1:size(avgOtherPoissNoise,1)
            if j~=i
                temp = temp+poissNoise(j,:);
            end
        end
        avgOtherPoissNoise(i,:) = temp/(size(avgOtherPoissNoise,1)-1);
    end
end
poissTdMeanSqErr = mean((avgOtherPoissNoise-poissNoise).^2);

clear temp;


%% Compare the poiss noise generated average with the average
for bin = 1:length(bins)-1
    indices = (r>=bins(bin))==(r<bins(bin+1));
    pBins(bin) = mean(p(indices));

    %     pBinsOptimum(bin) = mean(pOptimum(indices));
end

clear poissY RNDS i temp bins indices

%% Create IgorStructure
createIgor = 1;
%createIgor = input(sprintf('Make HDF5 file for data? 0 - no, 1 - yes...'));

% varRCorrected is var(r)-var(r(prestimulus))
varRCorrected =  var(r) - var(EpochCondition(sindex).AverageResponse(400:1999));

IgorStructure.r = r;
IgorStructure.stdR = sqrt(varRCorrected);
IgorStructure.r1 = truncatedEpochData(RandomIndex(size(truncatedEpochData,1)),:);
IgorStructure.r2 = truncatedEpochData(RandomIndex(size(truncatedEpochData,1)),:);
IgorStructure.r3 = truncatedEpochData(RandomIndex(size(truncatedEpochData,1)),:);
IgorStructure.r4 = truncatedEpochData(RandomIndex(size(truncatedEpochData,1)),:);

IgorStructure.y = y;
IgorStructure.stdY = std(y);
IgorStructure.s = s;
IgorStructure.shortH = shortH;
IgorStructure.tdMeanSqErr = tdMeanSqErr;
IgorStructure.modelTdMeanSqErr = modelTdMeanSqErr;
IgorStructure.tdImpulseMeanSqErr = tdImpulseMeanSqErr;
IgorStructure.responseSampleSize = size(truncatedEpochData,1);
IgorStructure.impulseSampleSize = size(truncatedFilterData,1);
IgorStructure.p = p;
IgorStructure.stdP = std(p);
IgorStructure.p1 = poissNoise(RandomIndex(size(poissNoise,1)),:);
IgorStructure.p2 = poissNoise(RandomIndex(size(poissNoise,1)),:);
IgorStructure.p3 = poissNoise(RandomIndex(size(poissNoise,1)),:);
IgorStructure.p4 = poissNoise(RandomIndex(size(poissNoise,1)),:);
IgorStructure.pDist = mean(poissNoise);
IgorStructure.stdPDist = std(IgorStructure.pDist);
IgorStructure.meanResponse = meanResponse;
IgorStructure.stimulus = stimulus;
IgorStructure.rBins = rBins;
IgorStructure.yBins = yBins;
IgorStructure.pBins = pBins;
IgorStructure.poissTdMeanSqErr = poissTdMeanSqErr;
IgorStructure.CellFile = CellInfo.CellFile;
IgorStructure.statNames = {'Unbiased mean squared error of traces scaled by(var(r))', ...
    'Variance in poiss traces scaled by (var(p))'};
IgorStructure.comparisonStatistics = [mean(mean((truncatedEpochData-avgOtherResponse).^2)),...
    mean(poissTdVariance*varRCorrected/var(p))];

%     IgorStructure.tdImpulseVariance = tdImpulseVarVect;
%     IgorStructure.H = H;
%     IgorStructure.S = S;
%     IgorStructure.Y = Y;
%     IgorStructure.h = h;
%     IgorStructure.meanSqErrAcrossTrials = meanSqErrAcrossTrials; % Mean square error across trials.  empirical
%     IgorStructure.meanSqErrAcrossTrialsVsModel =
%     meanSqErrAcrossTrialsVsModel; % Mean square error across trials compared to model.  should be minimized
%    IgorStructure.meanSqErrAcrossTrialsVsModelOptimum = meanSqErrAcrossTrialsVsModelOptimum';
%     IgorStructure.yOptimum = yOptimum;
%        IgorStructure.impulseScalingFactor = impulseScalingFactor;
%         mean(mean((truncatedEpochData-repmat(yOptimum,IgorStructure.responseSampleSize,1)).^2)),...
%         mean(mean((truncatedEpochData-repmat(pOptimum,IgorStructure.responseSampleSize,1)).^2)),...
%         mean(poissTdMeanSqErrOptimum)];
%     IgorStructure.stdP = std(poissNoise,1);
%     IgorStructure.pBinsOptimum = pBinsOptimum;
%     IgorStructure.yBinsOptimum = yBinsMatrix(impulseScalingFactorIndex,:);
%     IgorStructure.pOptimum = mean(poissNoiseOptimum);
%     IgorStructure.meanResponseOptimum = meanResponseOptimum;
%     IgorStructure.stdPOptimum = std(poissNoiseOptimum);
%     IgorStructure.poissonTDMeanSqrErrOptimum =
%     poissonTDMeanSqrErrOptimum;

%% View the results
% errorbar

view = 1;
if (view)
    if ~exist('IgorStructure','var')
        return
    end
    handle = figure;
    subplot(3,1,1);
    plot(TimeAxis(length(h)),h);
    subplot(3,1,2);
    plot(TimeAxis(length(s)),s);
    subplot(3,1,3);
    plot(TimeAxis(length(y)),y);
    disp(sprintf('Displaying lengthened impulse, stimulus, and filtered response...'));
    pause;
    subplot(111);
    plot(TimeAxis(length(IgorStructure.r1)),IgorStructure.r1/IgorStructure.stdR,'y');
    hold on;
    plot(TimeAxis(length(IgorStructure.r2)),IgorStructure.r2/IgorStructure.stdR,'y');
    plot(TimeAxis(length(IgorStructure.r3)),IgorStructure.r3/IgorStructure.stdR,'y');
    plot(TimeAxis(length(IgorStructure.r4)),IgorStructure.r4/IgorStructure.stdR,'y');
    plot(TimeAxis(length(IgorStructure.y)),IgorStructure.y/IgorStructure.stdY,'r');
    plot(TimeAxis(length(IgorStructure.r)),IgorStructure.r/IgorStructure.stdR,'k');
    disp(sprintf('Displaying filtered response and average response...'));
    pause;
    hold off;

    plot(TimeAxis(length(IgorStructure.r)),(IgorStructure.p1)/IgorStructure.stdP,'y');
    hold on;
    plot(TimeAxis(length(IgorStructure.r)),(IgorStructure.p2)/IgorStructure.stdP,'y');
    plot(TimeAxis(length(IgorStructure.r)),(IgorStructure.p3)/IgorStructure.stdP,'y');
    plot(TimeAxis(length(IgorStructure.r)),(IgorStructure.p4)/IgorStructure.stdP,'y');
    plot(TimeAxis(length(IgorStructure.r)),IgorStructure.pDist/IgorStructure.stdPDist,'c');
    plot(TimeAxis(length(IgorStructure.r)),IgorStructure.p/IgorStructure.stdP,'b');    
    plot(TimeAxis(length(IgorStructure.r)),IgorStructure.r/IgorStructure.stdR,'k');
    disp(sprintf('Displaying standard poisson noise, average of noise - scaled mean, and average response...'));
    pause;
    hold off;


    plot(TimeAxis(length(IgorStructure.r)),IgorStructure.tdMeanSqErr/(IgorStructure.stdR.*IgorStructure.stdR),'k');
    disp(sprintf('Displaying tdMeanSqErr...'));
    pause;
    hold on;
    plot(TimeAxis(length(IgorStructure.r)),IgorStructure.poissTdMeanSqErr/(IgorStructure.stdP.*IgorStructure.stdP),'b');

    disp(sprintf('Variance of r (black), variance of p (blue)...'));
    pause;
    hold off;
    plot(TimeAxis(length(r)),IgorStructure.r/IgorStructure.stdR,'k');
    hold on;
    plot(TimeAxis(length(r)),IgorStructure.r/IgorStructure.stdR...
        +sqrt(IgorStructure.tdMeanSqErr/IgorStructure.stdR^2)/sqrt(IgorStructure.responseSampleSize),'k');
    plot(TimeAxis(length(r)),IgorStructure.r/IgorStructure.stdR...
        -sqrt(IgorStructure.tdMeanSqErr/IgorStructure.stdR^2)/sqrt(IgorStructure.responseSampleSize),'k');
    plot(TimeAxis(length(r)),IgorStructure.p/IgorStructure.stdP,'b')
    disp(sprintf('Displaying standard error gaps...'));
    pause;
    hold off;
    bar(IgorStructure.comparisonStatistics);
    disp(sprintf('Number of response samples: %g . Number of Impulse Samples: %g...',...
        IgorStructure.responseSampleSize,IgorStructure.impulseSampleSize));
    pause;
   
    
end 
%% Now create the output structure
if (createIgor)
    condition = EpochCondition(sindex).Label{1};
    condition = condition(~isspace(condition));
    file = strcat(CellFile,'IgorStructure',condition,'.HDF5');
    IgorStructure.name = file;
    IgorStructure.root = condition;
    saveCommand = sprintf('save %s.mat IgorStructure',file);
    exportStructToHDF5(IgorStructure,IgorStructure.name,IgorStructure.root);
    eval(saveCommand);
end
clear condition file saveCommand createIgor

%% Clean Up
if(view)
    close(handle)
end

return
clear  s h handle sindex hindex cond view tdMeanSqErrVect modelTdMeanSqErrVect sindex hindex
