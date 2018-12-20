function[sigPowerPeakFreq, f1Mag, f2Mag, f2f1Ratio] = frequency_analysis(rhods, thetads, angleds, datarun002, ds, maxVal, concat, StimComb, sp, temp, path, binSize)
%preferred direction of all cells needs to come in

%spatial period = 64, temporal = 64, direction provided by theta
close all;

%initialize frequency, power and fft vectors
sigPower = [];
sigFreq = [];
samplingFreq = 1000; %get samples every ms, since spikes are at least few ms apart
stimLength = 8;
sampleTime = 1/samplingFreq;

if(concat == 1)
    sigVec = 1:1:datarun002.stimulus.repetitions*stimLength*samplingFreq;
    sigLength = length(sigVec);
    timeVec = (0:sigLength-1)*sampleTime;
    nextPowTwo = 2^nextpow2(sigLength); %Power of 2 closest to signal length for computational efficiency of fft algorithm
    sigFreq = samplingFreq/2*linspace(0,1,nextPowTwo/2+1); %only looking at half the spectrum since fft duplicates the spectrum

end

if(concat == 0)
    sigVec = 1:1:datarun002.stimulus.repetitions*stimLength*samplingFreq;
    sigLength = length(sigVec)/datarun002.stimulus.repetitions;
    timeVec = (0:sigLength-1)*sampleTime;
    nextPowTwo = 2^nextpow2(sigLength); %Power of 2 closest to signal length for computational efficiency of fft algorithm
    sigFreq = samplingFreq/2*linspace(0,1,nextPowTwo/2+1); %only looking at half the spectrum since fft duplicates the spectrum
end


frameRate = 120; %Our frame rate is 60 %Daniel FR = 120;
f1 = frameRate / temp; %initialize f1 and f2 given temporal period of stimulus
f2 = f1*2;

indexAng = [];
valueAng = [];
angleDir = [];
angleVecAve = [];
angleAll = [];
spikesPerTrial = [];
spikeTimes = [];
spikeVec = [];
sigPowerPeakVal = []; 
sigPowerPeakFreq = []; 
sigFft = [];
f1Mag = [];
f2Mag = [];
f2f1Ratio = [];
spikeVecMean = [];

%find row (temporal period) and column (spatial period)
spatPer = find(unique(StimComb(:,1))==sp);
tempPer = find(unique(StimComb(:,2))==temp);

%use angle that has maximum total average firing rate
if(maxVal == 1)
    [valueAng indexAng] = max(rhods{tempPer,spatPer}'); 
    angleDir = thetads{tempPer,spatPer}(1,indexAng);
    angleDir = angleDir*180/pi;
end

%use angle closest to that shown by vector sum
if(maxVal == 0)
    angleVecAve = angleds{tempPer,spatPer}*180/pi;%angle in which vector sum is pointing
    angleVecAve(angleVecAve<0) = angleVecAve(angleVecAve<0)+360;%remove negative angles
    angleAll = sort(thetads{tempPer,spatPer}(1,:))*180/pi;
    angleAll(1,length(angleAll)+1) = 360; %all angles of DG stimulus
    angleDir = interp1(angleAll,angleAll,angleVecAve, 'nearest'); %find angle closest to angleVec from 9 angles in angleAll
    angleDir(angleDir==360) = 0;
end


%Calculate frequency spectrums for all the cells from spike raster in
%direction found above
for j = 1:length(angleDir)
    stimLength = 8;
    for i = 1:length(datarun002.stimulus.combinations) %find trial number to get spike times from for each trial
        if((datarun002.stimulus.combinations(1,i).SPATIAL_PERIOD == sp) && (datarun002.stimulus.combinations(1,i).TEMPORAL_PERIOD == temp)&& (datarun002.stimulus.combinations(1,i).DIRECTION == angleDir(1,j)))
            trialnum = i;
        end
    end
    dsCellIndex = get_cell_indices(datarun002, ds(1,j)); %cell index from 1 to 400+
    [spikesPerTrial, psth, bins] = get_psth_sr(datarun002.spikes{dsCellIndex,1}, datarun002.stimulus.triggers(datarun002.stimulus.trial_list == trialnum), 'stop', stimLength, 'bin_size', binSize);  
    
    for h = 2:length(spikesPerTrial)
        spikesPerTrial{h,1} = spikesPerTrial{h,1}+stimLength;
        stimLength = stimLength + 8; %shift spike times by 8 for every trial so they can be concatenated into one vector of spikes
    end
    
    spikeTimes = cell2mat(spikesPerTrial); %concatenate all spike times to one vector
    spikeTimes = round(spikeTimes*(1000)); %since spike times are precise up to 1/20000s and sampling frequency is only every ms, round the spike times to nearest ms, but keep them as whole numbers because that's the only way ismember function finds them
    spikeVec = double(ismember(sigVec, spikeTimes)); %vector of zeroes and ones of spike instances
    stimLength = 8;
    
    if(concat == 0)
        lengthVec =  stimLength *samplingFreq;
        startVec = 1;
        for g = 1:(length(spikeVec))/lengthVec
            spikeVecMean(g,:) = spikeVec(1, startVec:lengthVec);
            startVec = startVec +  stimLength *samplingFreq;
            lengthVec = lengthVec +  stimLength *samplingFreq;
        end
        spikeVecMean = sum(spikeVecMean);
        spikeVec = [];
        spikeVec = spikeVecMean;
    end
        
spikeVec = spikeVec - mean(spikeVec); %subtract mean to remove DC component in Fourier transform
sigFft = fft(spikeVec,nextPowTwo)/sigLength; %calculate fft, scale by signal length
sigPower(j,:) = 2*abs(sigFft(1:nextPowTwo/2+1)); %only looking at 1/2 the spectrum - single sided spectrum so multiply sigFft power values by 2, also only look at absolute value (real part /magnitude of fft)
%plot(sigFreq,sigPower(j,:));
%plot(sigFreq(1:400),sigPower(j,1:400)); %plots a lower range of frequencies so you can see f1, f2, f3
% 
% 
% %path = '/Users/sravi/matlab/DS cell analysis/2012-10-15-0/Frequency Analysis/Concat PSTH/64ave/'
% %save_figure_pdf(path, num2str(ds(1,j)), gcf);
% 
%pause
end

%Where to take f1 and f2 from?? closest value to that in frequency vector
%or closest to peaks found by max function

[sigPowerPeakVal sigPowerPeakFreq] = max(sigPower');
sigPowerPeakFreq(2,:) = sigFreq(1, sigPowerPeakFreq(1,:));

[f1Diff f1Ind] = min(abs(f1 - sigPowerPeakFreq(2,:)));
[f2Diff f2Ind] = min(abs(f2 - sigPowerPeakFreq(2,:)));
f1Mag = sigPower(:,sigPowerPeakFreq(1,f1Ind));
f2Mag = sigPower(:,sigPowerPeakFreq(1,f2Ind));
f2f1Ratio = f2Mag./f1Mag;



%%%%%%%%%%%%%%%previous code using PSTH and not spike times%%%%%%%%%%%%%%%%

% % B = reshape(binned.', 1,[]);
% % psth = B./binSize; %?? necessary to divide by bin size or not?
% psthoff = psth - mean(psth); %need to subtract mean or not?
% n = nn;
% SigFFT = fft(psthoff, n); %Concatenated psth
% sigPower(j,:) = fftshift(abs(SigFFT)); % Rearranges output by moving zero frequency component to middle of spectrum
% F(j,:) = [-n/2:n/2-1];
% plot(F(j,:),sigPower(j,:));
% xlabel('frequency / fs');
% ylabel('Power');
% title(ds(1,j));
% 
% totalP(1,j) = trapz(F(j,:),sigPower(j,:));
% 

% 
% % bins = 0:0.1:64.7;
% % plot(bins, psth)
% % title(ds(1,j))
% 
% n = nn;
% %SigFFT = fft(p, n); %Mean PSTH
% SigFFT = fft(psthoff, n); %Concatenated psth
% %plot(abs(SigFFT));
% sigPower(j,:) = fftshift(abs(SigFFT)); % Rearranges output by moving zero frequency component to middle of spectrum
% F(j,:) = [-n/2:n/2-1];
% plot(F(j,:),sigPower(j,:));
% xlabel('frequency / fs');
% ylabel('Power');
% title(ds(1,j));
% % hold off;
% 
% totalP(1,j) = trapz(F(j,:),sigPower(j,:));

% for i = 1:size(sigPower,1)
%     k = 1;
%     for j = 1:size(sigPower,2)/10
%         auc(i,j) = sum(sigPower(i, k:k+9));
%         k = k+10;
%     end
% end
% 
% [sigPowerPeakVal sigPowerPeakFreq] = max(sigPower');
% sigPowerPeakFreq(2,:) = F(1, sigPowerPeakFreq(1,:));
% 
% f1 = ceil(mode(sigPowerPeakFreq(1,:))/10);
% f2 = ceil(find(F(1,:) == F(1,mode(sigPowerPeakFreq(1,:)))*2)/10);
% ff1 = mode(sigPowerPeakFreq(2,:));
% 
% if(isempty(f2))
%     F2F1 = 0;
% else
%     F2F1 = auc(:,f2)./auc(:,f1); 
% end
% 
% end



% e = sigPower(j,:);
% Y = e(e ~= a(1,j));
% [c(1,j) d(1,j)] = max(Y);
% 
% A(:,2) = auc(:,95)./auc(:,39);















