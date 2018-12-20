function[StimComb, f1Mag, f2Mag] = crg_analysis(datarun006, cellids)
close all;
StimComb = zeros(length(datarun006.stimulus.combinations), 4);
for i = 1:length(datarun006.stimulus.combinations)
    StimComb(i,1) = datarun006.stimulus.combinations(1,i).TEMPORAL_PERIOD;
    StimComb(i,2) = datarun006.stimulus.combinations(1,i).SPATIAL_PERIOD;
    StimComb(i,3) = datarun006.stimulus.combinations(1,i).SPATIAL_PHASE;
    StimComb(i,4) = datarun006.stimulus.combinations(1,i).SPATIAL_PHASE./datarun006.stimulus.combinations(1,i).SPATIAL_PERIOD;
end

% f1Mag = zeros(45,length(datarun006.stimulus.combinations));
% f2Mag = zeros(45,length(datarun006.stimulus.combinations));

f1Mag = zeros(length(cellids),length(datarun006.stimulus.combinations));
f2Mag = zeros(length(cellids),length(datarun006.stimulus.combinations));
sigFreq = [];
sigVec = [];
timeVec = [];
samplingFreq = 1000; %get samples every ms, since spikes are at least few ms apart
sampleTime = 1/samplingFreq;
stimLength = 8;
sigVec = 1:1:datarun006.stimulus.repetitions*stimLength*samplingFreq;
sigLength = length(sigVec);
timeVec = (0:sigLength-1)*sampleTime;
nextPowTwo = 2^nextpow2(sigLength); %Power of 2 closest to signal length for computational efficiency of fft algorithm
sigFreq = samplingFreq/2*linspace(0,1,nextPowTwo/2+1); %only looking at half the spectrum since fft duplicates the spectrum

frameRate = 60; %Our frame rate is 60 %Daniel FR = 120;
f1 = frameRate / datarun006.stimulus.params.TEMPORAL_PERIOD; %initialize f1 and f2 given temporal period of stimulus
f2 = f1*2;

[f1Val f1Freq] = min(abs(f1 -sigFreq));
[f2Val f2Freq] = min(abs(f2 -sigFreq));

binSize = 0.05;

for j = 1:length(cellids)
    for i = 1:length(datarun006.stimulus.combinations) %find trial number to get spike times from for each trial
        stimLength = 8;
        sigPower = [];
        spikeTimes = [];
        spikeVec = [];
        sigFft = [];
        sigPowerPeakFreq = [];
        spikesPerTrial = [];
        trialnum = i;
        dsCellIndex = get_cell_indices(datarun006, cellids(1,j)); %cell index from 1 to 400+

        [spikesPerTrial, psth, bins] = get_psth_sr(datarun006.spikes{dsCellIndex,1}, datarun006.stimulus.triggers(datarun006.stimulus.trial_list == trialnum), 'stop', stimLength, 'bin_size', binSize);  
        for h = 2:length(spikesPerTrial)
            spikesPerTrial{h,1} = spikesPerTrial{h,1}+stimLength;
            stimLength = stimLength + 8; %shift spike times by 8 for every trial so they can be concatenated into one vector of spikes
        end
    
    spikeTimes = cell2mat(spikesPerTrial); %concatenate all spike times to one vector
    spikeTimes = round(spikeTimes*(1000)); %since spike times are precise up to 1/20000s and sampling frequency is only every ms, round the spike times to nearest ms, but keep them as whole numbers because that's the only way ismember function finds them
    spikeVec = double(ismember(sigVec, spikeTimes)); %vector of zeroes and ones of spike instances
    
    spikeVec = spikeVec - mean(spikeVec); %subtract mean to remove DC component in Fourier transform
    sigFft = fft(spikeVec,nextPowTwo)/sigLength; %calculate fft, scale by signal length
    sigPower = 2*abs(sigFft(1:nextPowTwo/2+1)); %only look at first half of fft 
    
     %plot(sigFreq(1,1:400),sigPower(1,1:400));
     %pause;
    
    f1Mag(j,i) = sigPower(1,f1Freq);
    f2Mag(j,i) = sigPower(1,f2Freq);
    
    end
end

end

%% crg plots

all = 1:length(datarun006.stimulus.combinations);
[an indic] = sortrows(StimComb2, 4);%sort according to spatial phase
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);

F1F2s = cell(length(cellids),2);
for a = 1:length(cellids)
     c = 1;
     valsf1 = [];
     valsf2 = [];
     for b = 1:length(unique(StimComb2(:,2)));
        valsf1(b,:) = f1Mag(a, all(1,c:c+7));
         valsf2(b,:) = f2Mag(a, all(1,c:c+7));
         c = c+8;
     end
     F1F2s{a,1} = valsf1;
     F1F2s{a,2} = valsf2;
end

 
 x = unique(StimComb2(:,4));
 xx = unique(StimComb2(:,2));
 for a = 1:length(cellids)
     figure;
     for b = 1:11
     subplot(3,4,b)
     plot(x, F1F2s{a,1}(b,:), 'r')
     hold on;
     plot(x, F1F2s{a,2}(b,:), 'g')
     hold off;
     end
     f1SF = max(F1F2s{a,1}');
     f2SF = mean(F1F2s{a,2}');
     subplot(3,4,12)
     plot(xx, f1SF, 'r')
     hold on;
     plot(xx, f2SF, 'g')
     hold off;
     legend('F1', 'F2');
     save_figure_pdf('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/CRGPlots/', num2str(cellids(1,a)), gcf);
     close;
 end
 

 
 
 
 c = 1;
   for b = 1:length(unique(StimComb(:,2)));
        valsf1(b,:) = f1Mag(1, c:c+7);
         valsf2(b,:) = f2Mag(1,c:c+7);
         c = c+8;
     end
  F1F2s{1,1} = valsf1;
     F1F2s{1,2} = valsf2;
 
 
    x = unique(StimComb(:,4));
 xx = unique(StimComb(:,2));
 for a = 1
     figure;
     for b = 1:11
     subplot(3,4,b)
     plot(x, F1F2s{a,1}(b,:), 'r')
     hold on;
     plot(x, F1F2s{a,2}(b,:), 'g')
     hold off;
     end
     f1SF = max(F1F2s{a,1}');
     f2SF = mean(F1F2s{a,2}');
     subplot(3,4,12)
     plot(xx, f1SF, 'r')
     hold on;
     plot(xx, f2SF, 'g')
     hold off;
     legend('F1', 'F2');
     save_figure_pdf('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/CRGPlots/Rasters/4666/', num2str('4398-F1F2'), gcf);
     close;
 end 
     
     
     
     
     
     
     
     
     
     
     
     
 
 
 
 all = 1:length(datarun006.stimulus.combinations);
[an indic] = sortrows(StimComb2, 4);%sort according to spatial phase
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);


%subpin = [12 5 4 3 10 17 18 19]; %Places on subplot for each angle from 0 to 315 deg
get_rf_fit_radius(datarun000, 7352)
for a = 1:length(ss)
in1 = ismember(an(:,2), ss(1,a))'; % spatial period
SC2 = an(ismember(an(:,2), ss(1,a)), :);% spatial period
A1 = all(in1);%

h3 = figure;
for j = 1:length(A1)
    trigpre =  ismember(datarun006.stimulus.trial_list,A1(j));
    destaxes = subplot(2,4,j,'Parent',h3);
    spikesbytrials = get_raster(datarun006.spikes{get_cell_indices(datarun006, 7352),1}, datarun006.stimulus.triggers(trigpre), 0, 8, 0, 4, 'stop', 10, 'foa', destaxes, 'tic_color', [0 0 0]);
end

    save_figure_pdf('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/CRGPlots/Rasters/7352/', num2str(ss(1,a)), gcf);
    close;
end



for a = 1:length(cellids)
    plot(max(max(F1F2s{a,1}')), max(mean(F1F2s{a,2}')), 'o');
    hold on;
    %f1f2Ratio(1,a) = max(max(F1F2s{a,1}'))/max(mean(F1F2s{a,2}'));
end

for a = 1:length(cellids)
 f1f2Ratio(1,a) = max(mean(F1F2s{a,2}')/max(F1F2s{a,1}'));
end

 