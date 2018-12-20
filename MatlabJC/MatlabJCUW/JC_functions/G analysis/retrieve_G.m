FileName = 'Data/rotation/013107Hc3'

% open file
cd ~  %makes directory my own
[fp, error] = ITCInitializeAnalysis(95000, FileName);

% choose epoch
epoch = 5;

% read and plot measured response
figure(1);
[dat, error] = ITCReadEpoch(epoch, 0, fp);
subplot(1, 2, 1);
plot(dat);

% read and plot conductance stimuli used (raw conductances scaled by amplitudes)
[excg, inhg, error] = ITCReadEpochStmGClamp(epoch, 0, fp);
subplot(1, 2, 2);
plot(excg, 'r');
hold on;
plot(inhg);

