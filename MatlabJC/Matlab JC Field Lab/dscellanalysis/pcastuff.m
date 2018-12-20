cellids = intersect(intersect(datarun000.cell_ids, datarun001.cell_ids),datarun002.cell_ids);
addpath('/Users/sravi/matlab/Classification/');
global tcnormvr tcnormauc tcnormnorm tcnormmx tcnormminn tcnormsd minn mx mn vr normval auc sd mxt minnt isinormsd isinormauc isinormnorm isinormvr isinormmx isinormmn isinormminn acnormsd acnormauc acnormnorm acnormmx acnormminn magin anglein mxovmn mxtmintdiff aveSpTrig times amp pulsenormsd pulsenormnorm pulsenormvr pulsenormmx pulsenormmn pulsenormauc pulses pulsenormsdPSTH pulsenormaucPSTH pulsenormnormPSTH pulsenormvrPSTH pulsenormmxPSTH pulsenormmnPSTH pulsesPSTH A Aauc Anorm Amn Amx Asd Avr f1Mag f2Mag;



%Autocorrelation Functions, Receptive Field Diameters, Time Courses
%%%%%%%%%%%%%%%%%%%%%%%
wh = datarun001.triggers(1:4:length(datarun001.triggers), 1); %For 2012-10-10-1 dataset, that is how the triggers are arranges
gr = datarun001.triggers(2:4:length(datarun001.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001, get_cell_indices(datarun001,cellids), 0, '/0', wh, gr, 10, false);

binSize = 0.1:0.1:10; %change depending on length of trial
for a = 1:length(cellids)
 pulsePSTH(:,a) = sum(nhist{a,1})./50; %change depending on num of trials
end
for i = 1:length(cellids)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
 aucpulsePSTH(1,i) = trapz(abs(binSize), abs(pulsePSTH(:,i)));
end

sdpulsePSTH = std(pulsePSTH);
vrpulsePSTH = var(pulsePSTH);
meanpulsePSTH = mean(pulsePSTH);
[minpulsePSTH minpulsePSTHindex] = min(pulsePSTH);
[maxpulsePSTH maxpulsePSTHindex] = max(pulsePSTH);
mnovmxpulsePSTH = minpulsePSTH./maxpulsePSTH;

normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
aucpulsePSTH = repmat(aucpulsePSTH, size(pulsePSTH, 1), 1);
sdpulsePSTH = repmat(sdpulsePSTH, size(pulsePSTH, 1), 1);
vrpulsePSTH = repmat(vrpulsePSTH, size(pulsePSTH, 1), 1);
maxpulsePSTH = repmat(maxpulsePSTH, size(pulsePSTH, 1), 1);
meanpulsePSTH = repmat(meanpulsePSTH, size(pulsePSTH, 1), 1);
minpulsePSTH = repmat(minpulsePSTH, size(pulsePSTH, 1), 1);

pulsenormsdPSTH = pulsePSTH./sdpulsePSTH;
pulsenormaucPSTH = pulsePSTH./aucpulsePSTH;
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;
pulsenormvrPSTH = pulsePSTH./vrpulsePSTH;
pulsenormmxPSTH = pulsePSTH./maxpulsePSTH;
pulsenormmnPSTH = pulsePSTH./meanpulsePSTH;

pulsesPSTH = [minpulsePSTH(1,:); maxpulsePSTH(1,:); minpulsePSTHindex(1,:); maxpulsePSTHindex(1,:); mnovmxpulsePSTH(1,:);];

%4D PCA

%total number of spikes in each pulse in each trial
%calculate average number of spikes per pulse for each cell

aveSpTrig = zeros(4, length(sumSpTrTrig)); % 4 triggers - change if there is diferent trigger no
for n = 1:length(sumSpTrTrig)
 aveSpTrig(:, n) = sum(sumSpTrTrig{n,1})./length(sumSpTrTrig{n,1});
end

%Calculate all values we will be using to normalize
pulsenormnorm = zeros(size(sumSpTrTrig{1,1},2), length(sumSpTrTrig)); % 4 triggers - change if there is diferent trigger no
pulsenormsd = zeros(size(sumSpTrTrig{1,1},2), length(sumSpTrTrig)); % 4 triggers - change if there is diferent trigger no
pulsenormvr = zeros(size(sumSpTrTrig{1,1},2), length(sumSpTrTrig)); % 4 triggers - change if there is diferent trigger no
pulsenormmn = zeros(size(sumSpTrTrig{1,1},2), length(sumSpTrTrig));
pulsenormmx = zeros(size(sumSpTrTrig{1,1},2), length(sumSpTrTrig));
pulsenormauc = zeros(size(sumSpTrTrig{1,1},2), length(sumSpTrTrig));

for q = 1:length(aveSpTrig)
 pulsenormnorm(:,q) = aveSpTrig(:,q)./norm(aveSpTrig(:,q));
 pulsenormsd(:,q) = aveSpTrig(:,q)./std(aveSpTrig(:,q));
 pulsenormvr(:,q) = aveSpTrig(:,q)./var(aveSpTrig(:,q));
 pulsenormmn(:,q) = aveSpTrig(:,q)./mean(aveSpTrig(:,q));
 pulsenormmx(:,q) = aveSpTrig(:,q)./max(aveSpTrig(:,q));
 pulsenormauc(:,q) = aveSpTrig(:,q)./trapz(aveSpTrig(:,1));
end

[minpulse minpulseindex] = min(aveSpTrig);
[maxpulse maxpulseindex] = max(aveSpTrig);
mnovmxpulse = minpulse./maxpulse;

pulses = [minpulse(1,:); maxpulse(1,:); minpulseindex(1,:); maxpulseindex(1,:); mnovmxpulse(1,:);];


%%%%%%%%%%%%%%%%%%%%%%%%
[tc nontc] = get_time_courses_matrix(datarun000, cellids); %or cellids
x = 1:1:30;

for i = 1:length(cellids) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

sd = std(tc);
vr = var(tc);
[mx mxt] = max(tc);
mn = mean(tc);
[minn minnt] = min(tc);

mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
sd = repmat(sd, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
mn = repmat(mn, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);


%normalize time courses by each value
tcnormsd = tc./sd;
tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormmn = tc./mn;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

% %calculate time of zero crossing parameter
% timeZeroCross = cell(length(tc),1);
% 
% for j = 1:length(tc)
%  alltimes = [];
%  m = 1;
%  for k = 1:(size(tc,1)-1)
%  if((tc(k,j) >= 0 && tc(k+1,j) < 0) || (tc(k,j) < 0 && tc(k+1,j) >= 0))
%  alltimes(1,m) = k;
%  m = m+1;
%  end
%  end
%  timeZeroCross{j,1} = alltimes;
% end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datarun000 = get_autocorrelations(datarun000, cellids);
x1 = 0:1/2000:0.1; 
%nonds - cells not ds, tc - their time courses
%x1 = datarun000.autocorrelation{1, 1}.bins;

for i = 1:length(cellids) %or nonds
 ac(:,i) = datarun000.autocorrelation{get_cell_indices(datarun000, cellids(1,i)), 1}.probabilities;
 normvalacf(1, i) = norm( ac(:,i));
 aucacf(1,i) = trapz(abs(x1), abs(ac(:,i)));
end 

sdacf = std(ac);
vracf = var(ac);
mxacf= max(ac);
mnacf = mean(ac);
minnacf = min(ac);

aucacf = repmat(aucacf, size(ac, 1), 1);
normvalacf = repmat(normvalacf, size(ac, 1), 1);
sdacf = repmat(sdacf, size(ac, 1), 1);
vracf = repmat(vracf, size(ac, 1), 1);
mxacf = repmat(mxacf, size(ac, 1), 1);
mnacf = repmat(mnacf, size(ac, 1), 1);
minnacf = repmat(minnacf, size(ac, 1), 1);

acnormsd = ac./sdacf;
acnormauc = ac./aucacf;
acnormnorm = ac./normvalacf;
acnormvr = ac./vracf;
acnormmx = ac./mxacf;
acnormmn = ac./mnacf;
acnormminn = ac./minnacf;

%%%%%%%%%%%%%%%%%%%%%
%Interspikeintervals

datarun000 = get_interspikeinterval(datarun000, cellids);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses

for i = 1:length(cellids) %or nonds
 isi(:,i) = datarun000.interspikeinterval{get_cell_indices(datarun000, cellids(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
 aucisi(1,i) = trapz(abs(x2), abs(isi(:,i)));
end 

sdisi = std(isi);
vrisi = var(isi);
mxisi= max(isi);
mnisi = mean(isi);
minnisi = min(isi);

aucisi = repmat(aucisi, size(isi, 1), 1);
normvalisi = repmat(normvalisi, size(isi, 1), 1);
sdisi = repmat(sdisi, size(isi, 1), 1);
vrisi = repmat(vrisi, size(isi, 1), 1);
mxisi = repmat(mxisi, size(isi, 1), 1);
mnisi = repmat(mnisi, size(isi, 1), 1);
minnisi = repmat(minnisi, size(isi, 1), 1);

isinormsd = isi./sdisi;
isinormauc = isi./aucisi;
isinormnorm = isi./normvalisi;
isinormvr = isi./vrisi;
isinormmx = isi./mxisi;
isinormmn = isi./mnisi;
isinormminn = isi./minnisi;


%%%%%%%%%%%%%%%%%%%%%%%

[sigPowerPeakFreq, f1Mag, f2Mag, f2f1Ratio] = frequency_analysis(rhoin, thetain, anglein, datarun002, cellids, 1, 1, StimComb, 64, 256, 0, 1);
plot(f1Mag, f2Mag, 'o')

%%%%%%%%%%%%%%%
global tcnormvr tcnormauc tcnormnorm tcnormmx tcnormminn tcnormsd minn mx mn vr normval auc sd mxt minnt isinormsd isinormauc isinormnorm isinormvr isinormmx isinormmn isinormminn acnormsd acnormauc acnormnorm acnormmx acnormminn magin anglein mxovmn mxtmintdiff aveSpTrig times amp pulsenormsd pulsenormnorm pulsenormvr pulsenormmx pulsenormmn pulsenormauc pulses pulsenormsdPSTH pulsenormaucPSTH pulsenormnormPSTH pulsenormvrPSTH pulsenormmxPSTH pulsenormmnPSTH pulsesPSTH A Aauc Anorm Amn Amx Asd Avr f1Mag f2Mag;

[Classes_All] = classification(datarun000, cellids, get_cell_indices(datarun000, cellids));


[Classes_All] = classification(datarun000, nonds, get_cell_indices(datarun000, nonds));
[Classes_All] = classification(datarun000, normal, get_cell_indices(datarun000, normal));

[Classes_All] = classification(datarun000, on, get_cell_indices(datarun000, on));


{tcnormauc,acnormauc,tcnormnorm,acnormnorm,tcnormmx,acnormmx, amp, times, times, acnormsd, tcnormnorm, acnormsd}
% 
% {A A Aauc Aauc Anorm Anorm Amn Amn Amx Amx Asd Asd}
% {'Pulse addition'; 'Pulse addition-AUC-PC1,2'; 'Pulse addition-NORM-PC1,2'; 'Pulse addition-MN-PC1,2';'Pulse addition-MAX-PC1,2'; 'Pulse addition-SD-PC1,2'}
% [1 2 3 4 5 6 7 8 9 10 11 12]
% [1 2 1 2 1 2 1 2 1 2 1 2]

% {Aauc tcnormauc Anorm tcnormnorm Amx tcnormmx Asd tcnormsd Aauc acnormauc tcnormauc acnormauc}
% {'PulseAUC-TCAUC'; 'PulseNorm-TCNORM'; 'Pulsemx - TCmx'; 'PulseSD-TCsd'; 'PulseAUC-ACAUC'; 'TCAUC-ACAUC'}
% [1 2 3 4 5 6 7 8 9 10 11 12]
% [1 1 1 1 1 1 1 1 1 1 1 1]


{tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,tcnormmx,tcnormmx,amp,amp, times,times, magin{1,1}',magin{2,1}'}
{tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,tcnormmx,tcnormmx,amp,amp, times,times, magin{2,1}',magin{3,1}'}


{tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,tcnormmx,tcnormmx,tcnormauc,tcnormauc,tcnormnorm,tcnormnorm, magin{1,1}',magin{2,1}'}
{'TC-AUC-PC1,2'; 'TC-NORM-PC1,2'; 'TC-MX-PC1,2'; 'AMP-PC1,2';'TIMES-PC1,2'; 'MAG-1-2'}
[1 2 3 4 5 6 7 8 9 10]
[1 2 1 2 1 2 1 2 1 2]

{tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,tcnormmx,tcnormmx,amp,amp, times,times,tcnormvr, tcnormvr}
{'TC-AUC-PC1,2'; 'TC-NORM-PC1,2'; 'TC-MX-PC1,2'; 'AMP-PC1,2';'TIMES-PC1,2'; 'TC-NORM-VR-1,2'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 2 1 2 1 2 1 2 1 2]

{amp,amp,amp,amp,amp,amp, times,times,times,times,times,times}
{'AMP-PC1,2';'AMP-PC1,3';'AMP-PC2,3';'TIMES-PC1,2';'TIMES-PC1,3';'TIMES-PC2,3';}

[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 3 2 3 1 2 1 3 2 3]

{tcnormmx, tcnormmx, tcnormvr, tcnormmx, tcnormvr, tcnormauc, tcnormauc, tcnormauc, tcnormauc, tcnormnorm,tcnormvr, tcnormmx}
{'TC-MX-PC1,3'; 'TC-VR-MX-PC1,1';'TC-VR-AUC-PC2,2';'TC-AUC-PC1,3'; 'TC-AUC-NORM-PC2,2'; 'TC-VR-MX-PC2,2'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 3 1 1 2 2 1 3 2 2 2 2];

{tcnormauc,tcnormauc,tcnormauc,tcnormauc,tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,tcnormnorm,tcnormnorm,tcnormnorm,tcnormnorm}
{'TC-AUC-PC1,2'; 'TC-AUC-PC1,3'; 'TC-AUC-PC2,3'; 'TC-NORM-PC1,2'; 'TC-NORM-PC1,3'; 'TC-NORM-PC2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 3 2 3 1 2 1 3 2 3]

{tcnormvr,tcnormvr,tcnormvr,tcnormvr,tcnormvr,tcnormvr,tcnormmx, tcnormmx, tcnormmx,tcnormmx,tcnormmx,tcnormmx}
{'TC-VR-PC1,2'; 'TC-VR-PC1,3'; 'TC-VR-PC2,3'; 'TC-MX-PC1,2'; 'TC-MX-PC1,3'; 'TC-MX-PC2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 3 2 3 1 2 1 3 2 3]

{tcnormvr, tcnormmx, tcnormvr, tcnormauc, tcnormvr, tcnormnorm, tcnormmx, tcnormauc, tcnormmx, tcnormnorm, tcnormauc, tcnormnorm}
{'TC-VR-MX-PC1,1'; 'TC-VR-AUC-PC1,1'; 'TC-VR-NORM-PC1,1'; 'TC-MX-AUC-PC1,1'; 'TC-MX-NORM-PC1,1'; 'TC-AUC-NORM-PC1,1'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 1 1 1 1]


{tcnormvr, tcnormmx, tcnormvr, tcnormauc, tcnormvr, tcnormnorm, tcnormmx, tcnormauc, tcnormmx, tcnormnorm, tcnormauc, tcnormnorm}
{'TC-VR-MX-PC2,2'; 'TC-VR-AUC-PC2,2'; 'TC-VR-NORM-PC2,2'; 'TC-MX-AUC-PC2,2'; 'TC-MX-NORM-PC2,2'; 'TC-AUC-NORM-PC2,2'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[2 2 2 2 2 2 2 2 2 2 2 2]

{tcnormvr, tcnormmx, tcnormvr, tcnormauc, tcnormvr, tcnormnorm, tcnormmx, tcnormauc, tcnormmx, tcnormnorm, tcnormauc, tcnormnorm}
{'TC-VR-PC1,2'; 'TC-VR-PC1,3'; 'TC-VR-PC2,3'; 'TC-MX-PC1,2'; 'TC-MX-PC1,3'; 'TC-MX-PC2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[3 3 3 3 3 3 3 3 3 3 3 3]


{acnormsd,acnormsd, acnormsd,acnormsd, acnormsd,acnormsd, acnormauc,acnormauc, acnormauc,acnormauc,acnormauc,acnormauc}
{'AC-SD-PC1,2'; 'AC-SD-PC1,3'; 'AC-SD-PC2,3'; 'AC-AUC-PC1,2'; 'ACAUC-PC1,3'; 'AC-AUC-PC2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 3 2 3 1 2 1 3 2 3]

{acnormnorm,acnormnorm, acnormnorm,acnormnorm, acnormnorm,acnormnorm, acnormmx,acnormmx, acnormmx,acnormmx,acnormmx,acnormmx}
{'AC-NORM-PC1,2'; 'AC-NORM-PC1,3'; 'AC-NORM-PC2,3'; 'AC-MX-PC1,2'; 'ACMX-PC1,3'; 'AC-MX-PC2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 3 2 3 1 2 1 3 2 3]

{acnormsd,acnormauc, acnormauc,acnormnorm, acnormnorm,acnormmx, acnormmx,acnormauc,acnormnorm, acnormsd, acnormmx, acnormsd}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 1 1 1 1]

{acnormsd,acnormauc, acnormauc,acnormnorm, acnormnorm,acnormmx, acnormmx,acnormauc,acnormnorm, acnormsd, acnormmx, acnormsd}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[2 2 2 2 2 2 2 2 2 2 2 2]

{acnormsd,acnormauc, acnormauc,acnormnorm, acnormnorm,acnormmx, acnormmx,acnormauc,acnormnorm, acnormsd, acnormmx, acnormsd}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[3 3 3 3 3 3 3 3 3 3 3 3]

{tcnormvr, acnormsd, tcnormvr, acnormauc, tcnormvr, acnormnorm, tcnormvr, acnormmx, tcnormmx, acnormsd, tcnormmx, acnormauc}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 1 1 1 1]

{tcnormvr, acnormsd, tcnormvr, acnormauc, tcnormvr, acnormnorm, tcnormvr, acnormmx, tcnormmx, acnormsd, tcnormmx, acnormauc}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[2 2 2 2 2 2 2 2 2 2 2 2]

{tcnormvr, acnormsd, tcnormvr, acnormauc, tcnormvr, acnormnorm, tcnormvr, acnormmx, tcnormmx, acnormsd, tcnormmx, acnormauc}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[3 3 3 3 3 3 3 3 3 3 3 3]


{tcnormmx, acnormnorm, tcnormmx, acnormmx, tcnormauc, acnormsd, tcnormauc, acnormauc, tcnormauc, acnormnorm, tcnormauc, acnormmx}
{'TC-MX-AC-NORM-PC1,1'; 'TC-MX-AC-MX-PC1,1'; 'TC-AUC-AC-SD-PC1,1'; 'TC-AUC-AC-AUC-PC1,1';'TC-AUC-AC-NORM-PC1,1'; 'TC-AUC-AC-MX-PC1,1'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 1 1 1 1]

{tcnormmx, acnormnorm, tcnormmx, acnormmx, tcnormauc, acnormsd, tcnormauc, acnormauc, tcnormauc, acnormnorm, tcnormauc, acnormmx}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[2 2 2 2 2 2 2 2 2 2 2 2]

{tcnormmx, acnormnorm, tcnormmx, acnormmx, tcnormauc, acnormsd, tcnormauc, acnormauc, tcnormauc, acnormnorm, tcnormauc, acnormmx}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[3 3 3 3 3 3 3 3 3 3 3 3]

{tcnormnorm, acnormsd, tcnormnorm, acnormauc, tcnormnorm, acnormnorm, tcnormnorm, acnormmx, tcnormnorm, acnormsd, tcnormnorm, acnormauc}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 2 2 2 2]

{tcnormnorm, acnormnorm, tcnormnorm, acnormmx,tcnormnorm, acnormsd, tcnormnorm, acnormauc, tcnormnorm, acnormnorm, tcnormnorm, acnormmx}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'AC-AUC-PC2,3'; 'AC-MAX2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12] 
[2 2 2 2 3 3 3 3 3 3 3 3]
 








{isinormauc,isinormauc,isinormauc,isinormauc,isinormauc,isinormauc,isinormnorm,isinormnorm,isinormnorm,isinormnorm,isinormnorm,isinormnorm}
{'ISI-AUC-PC1,2'; 'ISI-AUC-PC1,3'; 'ISI-AUC-PC2,3'; 'ISI-NORM-PC1,2'; 'ISI-NORM-PC1,3'; 'ISI-NORM-PC2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 3 2 3 1 2 1 3 2 3]

{isinormvr,isinormvr,isinormvr,isinormvr,isinormvr,isinormvr,isinormmx, isinormmx, isinormmx,isinormmx,isinormmx,isinormmx}
{'ISI-VR-PC1,2'; 'ISI-VR-PC1,3'; 'ISI-VR-PC2,3'; 'ISI-MX-PC1,2'; 'ISI-MX-PC1,3'; 'ISI-MX-PC2,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 3 2 3 1 2 1 3 2 3]

{isinormvr, isinormmx, isinormvr, isinormauc, isinormvr, isinormnorm, isinormmx, isinormauc, isinormmx, isinormnorm, isinormauc, isinormnorm}
{'ISI-VR-MX-PC1,1'; 'ISI-VR-AUC-PC1,1'; 'ISI-VR-NORM-PC1,1'; 'ISI-MX-AUC-PC1,1'; 'ISI-MX-NORM-PC1,1'; 'ISI-AUC-NORM-PC1,1'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 1 1 1 1]

{isinormvr, isinormmx, isinormvr, isinormauc, isinormvr, isinormnorm, isinormmx, isinormauc, isinormmx, isinormnorm, isinormauc, isinormnorm}
{'ISI-VR-MX-PC2,2'; 'ISI-VR-AUC-PC2,2'; 'ISI-VR-NORM-PC2,2'; 'ISI-MX-AUC-PC2,2'; 'ISI-MX-NORM-PC2,2'; 'ISI-AUC-NORM-PC2,2'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[2 2 2 2 2 2 2 2 2 2 2 2]

{isinormvr, isinormmx, isinormvr, isinormauc, isinormvr, isinormnorm, isinormmx, isinormauc, isinormmx, isinormnorm, isinormauc, isinormnorm}
{'ISI-VR-MX-PC3,3'; 'ISI-VR-AUC-PC3,3'; 'ISI-VR-NORM-PC3,3'; 'ISI-MX-AUC-PC3,3'; 'ISI-MX-NORM-PC3,3'; 'ISI-AUC-NORM-PC3,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[3 3 3 3 3 3 3 3 3 3 3 3]




{tcnormmx, isinormnorm, tcnormmx, isinormmx, tcnormauc, isinormsd, tcnormauc, isinormauc, tcnormauc, isinormnorm, tcnormauc, isinormmx}
{'TC-MX-ISI-NORM-PC1,1'; 'TC-MX-ISI-MX-PC1,1'; 'TC-AUC-ISI-SD-PC1,1'; 'TC-AUC-ISI-AUC-PC1,1';'TC-AUC-ISI-NORM-PC1,1'; 'TC-AUC-ISI-MX-PC1,1'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 1 1 1 1]

{tcnormmx, isinormnorm, tcnormmx, isinormmx, tcnormauc, isinormsd, tcnormauc, isinormauc, tcnormauc, isinormnorm, tcnormauc, isinormmx}
{'TC-MX-ISI-NORM-PC2,2'; 'TC-MX-ISI-MX-PC2,2'; 'TC-AUC-ISI-SD-PC2,2'; 'TC-AUC-ISI-AUC-PC2,2';'TC-AUC-ISI-NORM-PC2,2'; 'TC-AUC-ISI-MX-PC2,2'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[2 2 2 2 2 2 2 2 2 2 2 2]

{tcnormmx, isinormnorm, tcnormmx, isinormmx, tcnormauc, isinormsd, tcnormauc, isinormauc, tcnormauc, isinormnorm, tcnormauc, isinormmx}
{'TC-MX-ISI-NORM-PC3,3'; 'TC-MX-ISI-MX-PC3,3'; 'TC-AUC-ISI-SD-PC3,3'; 'TC-AUC-ISI-AUC-PC3,3';'TC-AUC-ISI-NORM-PC3,3'; 'TC-AUC-ISI-MX-PC3,3'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[3 3 3 3 3 3 3 3 3 3 3 3]

{tcnormnorm, isinormsd, tcnormnorm, isinormauc, tcnormnorm, isinormnorm, tcnormnorm, isinormmx, tcnormnorm, isinormsd, tcnormnorm, isinormauc}
{'TC-NORM-ISI-SD-PC1,1'; 'TC-NORM-ISI-AUC-PC1,1'; 'TC-NORM-ISI-NORM-PC1,1'; 'TC-NORM-ISI-MX-PC1,1';'TC-NORM-ISI-SD-PC2,2'; 'TC-NORM-ISI-AUC-PC2,2'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 2 2 2 2]

{tcnormnorm, isinormnorm, tcnormnorm, isinormmx,tcnormnorm, isinormsd, tcnormnorm, isinormauc, tcnormnorm, isinormnorm, tcnormnorm, isinormmx}
{'TC-NORM-ISI-NORM-PC2,2'; 'TC-NORM-ISI-MX-PC2,2'; 'TC-NORM-ISI-SD-PC3,3'; 'TC-NORM-ISI-AUC-PC3,3';'TC-NORM-ISI-NORM-PC3,3'; 'TC-NORM-ISI-MX-PC3,3'}
[1 2 3 4 5 6 7 8 9 10 11 12] 
[2 2 2 2 3 3 3 3 3 3 3 3]













{tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,tcnormmx,tcnormmx,amp,amp, times,times, f1Mag, f2Mag}
{'TC-AUC-PC1,2'; 'TC-NORM-PC1,2'; 'TC-MX-PC1,2'; 'AMP-PC1,2';'TIMES-PC1,2'; 'MAG-2-3'}
[1 2 3 4 5 6 7 8 9 10]
[1 2 1 2 1 2 1 2 1 2]



% {acnormsd,acnormsd, acnormauc,acnormauc, acnormnorm,acnormnorm, acnormmx,acnormmx,acnormauc, tcnormauc, acnormmx, tcnormauc}
% {'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'PULSE-PC1,2'; 'TC-AUC1-2'}
% [1 2 3 4 5 6 7 8 9 10 11 12]
% [1 2 1 2 1 2 1 2 1 1 1 1]


% {acnormnorm,acnormnorm, acnormmx,acnormmx,acnormnorm, tcnormnorm, acnormmx, tcnormmx, acnormauc,tcnormauc, acnormnorm, tcnormmx}
% {'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'PULSE-PC1,2'; 'TC-AUC1-2'}
% [1 2 3 4 5 6 7 8 9 10 11 12]
% [1 2 1 2 1 1 1 1 1 1 1 1]


{tcnormvr tcnormvr,tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,tcnormmx,tcnormmx,tcnormminn, tcnormminn, mxt(1,:)',minnt(1,:)'};
[1 2 3 4 5 6 7 8 9 10]
[1 2 1 2 1 2 1 2 1 2]





{acnormnorm,tcnormnorm, acnormmx,tcnormmx,acnormnorm, pulsenormsdPSTH, acnormmx, pulsesPSTH, acnormauc,pulsenormsdPSTH, acnormnorm, pulsesPSTH}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'PULSE-PC1,2'; 'TC-AUC1-2'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 1 1 1 1]

{pulsenormsdPSTH, tcnormnorm, pulsesPSTH, tcnormnorm, pulsenormsdPSTH , tcnormmx, pulsesPSTH, tcnormmx, pulsenormsdPSTH, pulsenormsdPSTH, pulsesPSTH, pulsesPSTH}







{acnormsd,acnormsd, acnormauc,acnormauc, acnormnorm,acnormnorm, acnormmx,acnormmx,aveSpTrig, aveSpTrig, tcnormauc, tcnormauc}
{'AC-SD-PC1,2'; 'AC-AUC-PC1,2'; 'AC-NORM-PC1,2'; 'AC-MAX1,2';'PULSE-PC1,2'; 'TC-AUC1-2'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 2 1 2 1 2 1 2 1 2]

% pulsenormauc pulsenormnormPSTH pulsenormmnPSTH 




{pulsenormsd, pulsenormsd, pulsenormnorm, pulsenormnorm, pulsenormvr, pulsenormvr, pulsenormmx, pulsenormmx, pulsenormmn, pulsenormmn, pulses, pulses}
{'PULSE-SD-PC1,2'; 'PULSE-NORM-PC1,2'; 'PULSE-VR-PC1,2'; 'PULSE-MAX1,2';'PULSE-NORM1,2'; 'PULSES-MAX-MINS'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 2 1 2 1 2 1 2 1 2]

{pulsenormsdPSTH, pulsenormsdPSTH, pulsenormnormPSTH, pulsenormnormPSTH, pulsenormvrPSTH, pulsenormvrPSTH, pulsenormmxPSTH, pulsenormmxPSTH, pulsenormaucPSTH, pulsenormaucPSTH, pulsesPSTH, pulsesPSTH}
{'PULSE-SD-PC1,2'; 'PULSE-NORM-PC1,2'; 'PULSE-VR-PC1,2'; 'PULSE-MAX1,2';'PULSE-NORM1,2'; 'PULSES-MAX-MINS'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 2 1 2 1 2 1 2 1 2 1 2]

{pulsenormsdPSTH, tcnormnorm, pulsesPSTH, tcnormnorm, pulsenormsdPSTH , tcnormmx, pulsesPSTH, tcnormmx, pulsenormsdPSTH, pulsenormsdPSTH, pulsesPSTH, pulsesPSTH}
{'PULSE-SD-PC1,2'; 'PULSE-NORM-PC1,2'; 'PULSE-VR-PC1,2'; 'PULSE-MAX1,2';'PULSE-NORM1,2'; 'PULSES-MAX-MINS'}
[1 2 3 4 5 6 7 8 9 10 11 12]
[1 1 1 1 1 1 1 1 1 2 1 2]

% {tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,tcnormmx,tcnormmx,amp,amp, times,times,tcnormvr, tcnormvr}

%NEXT THING TO TRY AFTER ARVO: COMBINATIONS OF PARAMS
%MORE PARAMS
%STA FIT
%NON LINEARITY
%STV
%RF FIT

{tcnormvr tcnormvr,tcnormauc,tcnormauc,tcnormnorm,tcnormnorm,anglein{2,1}', anglein{3,1}',magin{2,:}',magin{3,:}',minn(1,:)', mx(1,:)'}
{'TC-VR-PC1,2'; 'TC-AUC-PC1,2'; 'TC-NORM-PC1,2'; 'ANGLE2-3';'MAG2-3'; 'MIN-VS-MAX'};
[1 2 3 4 5 6];
[1 2 1 2 1 2];
 
{mxovmn(1,:)',mxtmintdiff(1,:)',tcnormauc,tcnormauc,aveSpTrig,aveSpTrig,tcnormmx,tcnormmx,magin{2,:}',magin{3,:}',acnormnorm, acnormnorm}
{'MX/MN - VS MXT-MINT'; 'TC-AUC-PC1,2'; 'PULSE1-2'; 'TCNORMMX1-2';'MAG-2-3'; 'ACF-NORM1-2'}
[3 4 5 6 7 8 11 12]
[1 2 1 2 1 2 1 2]


{acnormsd,acnormsd, acnormauc,acnormauc, acnormnorm,acnormnorm, acnormmx,acnormmx,aveSpTrig, aveSpTrig, tcnormauc, tcnormauc}

[minn(1,:)', mn(1,:)', minn(1,:)', vr(1,:)', minn(1,:)', normval(1,:)', minn(1,:)', auc(1,:)', sd(1,:)', mn(1,:)', sd(1,:)', mx(1,:)']


% [a b d] = intersect(c, cellids);
% {acnormsd(:,d),acnormsd(:,d), acnormauc(:,d),acnormauc(:,d), acnormnorm(:,d),acnormnorm(:,d), acnormmx(:,d),acnormmx(:,d),aveSpTrig(:,d), aveSpTrig(:,d), tcnormauc(:,d), tcnormauc(:,d)}




%define PCA variables as global variables outside here and in function!!!!!
 
 
 %other variables to plot:
 
 minn(1,:), mn(1,:)
 minn(1,:), vr(1,:)
 minn(1,:), normval(1,:)
 minn(1,:), auc(1,:)
 sd(1,:), mn(1,:)
 sd(1,:), mx(1,:)
 sd(1,:), auc(1,:)
 
 mn(1,:), mx(1,:)
mn(1,:), vr(1,:)
mn(1,:), auc(1,:)

mx(1,:), vr(1,:)
mx(1,:), normval(1,:)
mx(1,:), auc(1,:)
vr(1,:),auc(1,:)
auc(1,:),normval(1,:)
mxt(1,:), minnt(1,:)

 c = 1;
b = [];

for a = 1:length(Classes_All)
 if(strfind(Classes_All{a,1}, 'on'))
 b(c,1) = a;
 c = c+1;
 end
end

A = cellids(1,b)
 A = nonds(1,b)


ont3 = cellids(1,b)
A = cellids(1,b)

plot_rf_summaries(datarun000, A,'label', true)
save_figure_pdf('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/AllPlots/ON/t3/', 'RF', gcf)
plot_rf_portraits(datarun000,A)
save_figure_pdf('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/AllPlots/ON/t3/', 'RF Portraits', gcf)
plot_time_courses(datarun000, A, 1, 1);
save_figure_pdf('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/AllPlots/ON/t3/', 'Time Courses', gcf)
plot_autocorrelograms(datarun000, A, 'foa', 0);
save_figure_pdf('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/AllPlots/ON/t3/', 'ACFs', gcf)


A = [93 485 661 2297 2418 3151 4295 4879 5147 5297 6077 6079 6333 7339];

%2013-04-03-0 Rat WT dataset
ds = [50 541 980 1006 1592 1773 1924 3935 4231 4791 5150 5312 5358 5509 6034 6453 6618];
offt1 = [32 107 211 304 408 438 631 736 887 1129 1143 1186 1292 1443 1458 1654 1668 1741 1848 2090 2553 2732 2928 3049 3094 3573 3843 4051 4055 4202 4263 4307 4442 4457 4696 4742 4831 5102 5133 5239 5298 5446 5672 5791 5822 5883 5987 6169 6425 6571 6632 6721 6767 6827 6856 6916 7006 7010 7098 7174 7489 7518 7651];
offt2 = [272 482 1426 1817 1893 2582 3185 3422 3766 3858 4172 4366 5116 5283 5389 5600 5824 5836 5957 6211 6241 6391 6483 6665 7606];
offt3 = [691 1070 1264 1861 2029 3874 5296 5716 5911 5956 6274 6814];
offt4 = [407 571 1970 4277 4531 4966 5840 5928 7159 7442];
offt5 = [453 556 634 1101 1414 1700 3183 4039 4112 4128 4381 4456 4578 4609 4756 4846 4938 5090 5149 5611 6032 6858];
offt6 = [1786 2088 2822 2851 3752 7069 7201 7577];

ont1 = [3 197 226 348 378 436 601 646 1160 1411 1520 1831 2027 2104 2491 2866 3693 4036 4187 4261 4308 4351 4352 4504 4847 4880 4952 5028 5118 5596 6544 7070 7083 7336 7338 7503 7595];
ont2 = [76 93 196 486 904 934 1021 1847 1863 2416 2749 3526 4416 4488 5208 5583 5747 6257 6410 6560 6620 7607];
ont3 = [ 122 258 273 542 587 616 678 721 1053 1068 1100 1295 1459 1699 1726 1833 1877 2058 2134 2194 2644 2793 2957 3226 3558 3871 3887 3889 3962 4099 4113 4490 4577 4622 4623 4713 4924 4937 4968 5042 5073 5371 5462 5476 5659 5839 5913 6002 6106 6181 6228 6467 6481 6542 6587 6800 6828 6889 6979 7413 7428 7441 7474 7563 7638]; %big group
ont4 = [64 287 1475 2059 2476 2613 3468 6980]; %questionable
ont5 = [1325 1698 1922 2763 2882 3244 3246 3828 4561]; %questionable

%2012-10-15-0 dataset - all intersecting cells only
ds = [257,301,333,438,467,708,1097,1595,1683,1685,1895,2042,2449,2898,3512,3636,3736,3815,3842,3934,4069,4157,4173,4353,4731,4789,4846,4985,5150,5632,5702,5719,6229,6321,6332,6751,6797,7308];
nonds = [4,31,46,62,76,94,154,182,226,272,347,407,424,454,496,514,528,586,649,692,751,753,766,768,782,783,843,857,860,872,888,889,901,979,991,1051,1067,1081,1084,1098,1111,1127,1130,1156,1191,1246,1310,1339,1381,1384,1385,1398,1411,1427,1486,1532,1549,1564,1578,1581,1591,1637,1652,1684,1786,1817,1876,1877,1892,1893,1908,1921,1966,1969,1999,2011,2026,2136,2161,2177,2192,2206,2208,2236,2253,2311,2326,2328,2343,2357,2373,2401,2417,2461,2462,2478,2521,2522,2555,2597,2686,2716,2732,2746,2794,2809,2851,2868,2881,2896,2929,3002,3061,3091,3139,3152,3198,3200,3226,3241,3244,3258,3286,3287,3303,3319,3422,3452,3530,3559,3586,3589,3601,3634,3635,3648,3679,3691,3692,3695,3721,3767,3812,3813,3857,3859,3889,3917,3933,3935,3946,3991,3994,4006,4021,4022,4096,4097,4098,4130,4145,4188,4231,4234,4235,4246,4278,4279,4294,4324,4427,4442,4459,4486,4487,4501,4503,4562,4591,4668,4697,4732,4771,4774,4788,4864,4892,4941,4997,4998,4999,5071,5073,5088,5104,5116,5148,5179,5223,5359,5405,5433,5446,5464,5567,5569,5641,5645,5657,5658,5672,5703,5705,5733,5791,5836,5851,5853,5896,5898,5927,5941,6034,6064,6093,6106,6125,6139,6140,6152,6155,6170,6196,6257,6260,6286,6304,6361,6363,6376,6380,6391,6392,6422,6439,6451,6455,6512,6542,6589,6721,6722,6737,6752,6811,6826,6828,6886,6903,6931,6976,6980,6992,7021,7040,7067,7069,7096,7157,7186,7203,7234,7261,7278,7306,7354,7442,7471,7475,7487,7503,7517,7520,7532,7562,7667];
all = [4,31,46,62,76,94,154,182,226,257,272,301,333,347,407,424,438,454,467,496,514,528,586,649,692,708,751,753,766,768,782,783,843,857,860,872,888,889,901,979,991,1051,1067,1081,1084,1097,1098,1111,1127,1130,1156,1191,1246,1310,1339,1381,1384,1385,1398,1411,1427,1486,1532,1549,1564,1578,1581,1591,1595,1637,1652,1683,1684,1685,1786,1817,1876,1877,1892,1893,1895,1908,1921,1966,1969,1999,2011,2026,2042,2136,2161,2177,2192,2206,2208,2236,2253,2311,2326,2328,2343,2357,2373,2401,2417,2449,2461,2462,2478,2521,2522,2555,2597,2686,2716,2732,2746,2794,2809,2851,2868,2881,2896,2898,2929,3002,3061,3091,3139,3152,3198,3200,3226,3241,3244,3258,3286,3287,3303,3319,3422,3452,3512,3530,3559,3586,3589,3601,3634,3635,3636,3648,3679,3691,3692,3695,3721,3736,3767,3812,3813,3815,3842,3857,3859,3889,3917,3933,3934,3935,3946,3991,3994,4006,4021,4022,4069,4096,4097,4098,4130,4145,4157,4173,4188,4231,4234,4235,4246,4278,4279,4294,4324,4353,4427,4442,4459,4486,4487,4501,4503,4562,4591,4668,4697,4731,4732,4771,4774,4788,4789,4846,4864,4892,4941,4985,4997,4998,4999,5071,5073,5088,5104,5116,5148,5150,5179,5223,5359,5405,5433,5446,5464,5567,5569,5632,5641,5645,5657,5658,5672,5702,5703,5705,5719,5733,5791,5836,5851,5853,5896,5898,5927,5941,6034,6064,6093,6106,6125,6139,6140,6152,6155,6170,6196,6229,6257,6260,6286,6304,6321,6332,6361,6363,6376,6380,6391,6392,6422,6439,6451,6455,6512,6542,6589,6721,6722,6737,6751,6752,6797,6811,6826,6828,6886,6903,6931,6976,6980,6992,7021,7040,7067,7069,7096,7157,7186,7203,7234,7261,7278,7306,7308,7354,7442,7471,7475,7487,7503,7517,7520,7532,7562,7667];

on = [4,31,154,182,272,424,454,496,528,692,766,783,857,860,888,979,1067,1111,1339,1385,1398,1486,1532,1549,1578,1591,1786,1876,1877,1892,1969,2011,2026,2136,2161,2208,2311,2328,2417,2461,2597,2716,2851,2929,3002,3152,3198,3241,3244,3303,3319,3422,3452,3530,3559,3586,3634,3691,3721,3812,3859,3917,3946,3991,3994,4022,4096,4130,4145,4231,4246,4279,4427,4501,4562,4668,4774,4892,4941,4997,4998,5073,5088,5104,5179,5359,5405,5567,5657,5791,5853,5898,5927,5941,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6542,6722,6737,6752,6826,6903,6980,7067,7096,7157,7186,7203,7261,7354,7442,7487,7520,7532,7562];
off = [46,62,76,94,226,347,407,514,586,649,751,753,768,782,843,872,889,901,991,1051,1081,1084,1098,1127,1130,1156,1191,1246,1310,1381,1384,1411,1427,1564,1581,1637,1652,1684,1817,1893,1908,1921,1966,1999,2177,2192,2206,2236,2253,2326,2343,2357,2373,2401,2462,2478,2521,2522,2555,2686,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3287,3589,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3933,3935,4006,4021,4097,4098,4188,4234,4235,4278,4294,4324,4442,4459,4486,4487,4503,4591,4697,4732,4771,4788,4864,4999,5071,5116,5148,5223,5433,5446,5464,5569,5641,5645,5658,5672,5703,5705,5733,5836,5851,5896,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6391,6422,6451,6589,6721,6811,6828,6886,6931,6976,6992,7021,7040,7069,7234,7278,7306,7471,7475,7503,7517,7667];

onwithsnrcutoff = [4,31,154,182,272,454,496,528,692,766,783,857,860,888,979,1067,1111,1339,1385,1398,1486,1532,1549,1578,1591,1786,1876,1877,1892,1969,2011,2026,2136,2161,2208,2311,2328,2417,2461,2597,2716,2851,2929,3002,3152,3198,3241,3244,3303,3319,3422,3452,3530,3559,3586,3634,3691,3721,3812,3859,3917,3946,3991,3994,4022,4096,4130,4145,4231,4246,4427,4501,4562,4668,4774,4892,4941,4997,4998,5073,5088,5104,5179,5359,5405,5567,5657,5791,5853,5898,5927,5941,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6542,6722,6737,6752,6826,6903,6980,7067,7096,7157,7186,7203,7261,7354,7442,7487,7520,7532,7562];
offwithsnrcutoff = [46,62,76,94,226,347,407,514,586,649,751,753,768,782,843,872,889,901,991,1051,1081,1084,1098,1127,1130,1156,1191,1246,1310,1381,1384,1411,1427,1564,1581,1637,1652,1684,1817,1893,1908,1921,1966,1999,2177,2192,2206,2236,2253,2326,2343,2357,2373,2401,2462,2478,2521,2522,2555,2686,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3287,3589,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3933,3935,4006,4021,4097,4098,4188,4234,4235,4278,4294,4324,4442,4459,4486,4487,4503,4591,4697,4732,4771,4788,4864,4999,5071,5116,5148,5223,5433,5446,5464,5569,5641,5645,5658,5672,5703,5705,5733,5836,5851,5896,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6391,6422,6451,6589,6721,6811,6828,6886,6931,6976,6992,7021,7040,7069,7234,7278,7306,7471,7503,7517,7667];

%upated on and off only on june 18th

ont21ststage = [272,424,783,888,979,1067,1084,1398,1486,1876,2026,2136,2417,2597,2929,3198,3244,3303,3559,3634,3917,3946,3991,4246,4562,4668,4997,4998,5657,5705,5898,5927,6106,6170,6455,6512,6722,6903,7067,7069,7203,7354,7562,454,496]; %45 cells
ont11ststage = [4,31,154,454,496,528,692,860,1111,1339,1578,1786,1877,1892,1969,2011,2161,2208,2461,2716,2851,3002,3152,3319,3422,3530,3691,3721,3994,4096,4145,4427,4501,4774,5073,5088,5104,5359,5567,5791,5853,5941,6155,6542,6826,7096,7261,7442,7487,7532];
ont1getrid1ststage = [7096, 5791, 31, 454 496 528 1578 2716 3152 3530 4145 6155];
%91 cells = [182,272,424,649,766,783,857,888,979,1067,1084,1385,1398,1486,1532,1549,1591,1876,2026,2136,2236,2311,2328,2417,2597,2929,3198,3241,3244,3303,3452,3559,3586,3634,3812,3859,3917,3946,3991,4022,4130,4231,4246,4562,4668,4892,4941,4997,4998,5179,5405,5657,5705,5898,5927,6093,6106,6125,6152,6170,6304,6363,6392,6439,6455,6512,6722,6737,6752,6903,6980,7067,7069,7157,7186,7203,7354,7520,7562,7096,5791,31,454,496,528,1578,2716,3152,3530,4145,6155];
%53cells: [182,424,649,766,857,1084,1385,1532,1549,1591,2136,2236,2311,2328,3241,3452,3586,3812,3859,4022,4130,4231,4892,4941,4997,5179,5405,5705,6093,6125,6152,6304,6363,6392,6439,6737,6752,6980,7069,7157,7186,7520,7096,5791,31,496,528,1578,2716,3152,3530,4145,6155]


%123 cells = [62,76,94,226,347,407,514,586,753,843,872,889,901,991,1098,1127,1130,1156,1191,1246,1310,1384,1411,1427,1564,1581,1684,1817,1893,1908,1921,1999,2192,2206,2253,2326,2343,2357,2373,2462,2478,2521,2522,2555,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3935,4097,4098,4188,4234,4235,4278,4279,4294,4442,4459,4486,4487,4697,4732,4771,4788,4864,4999,5071,5116,5148,5433,5464,5569,5645,5672,5703,5733,5851,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6422,6589,6721,6828,6886,6931,6992,7040,7234,7278,7475,7503,7517,7667];
offt1 = [46,751,768,782,1051,1081,1381,1637,1652,1966,2177,2401,2686,3287,3589,3933,4006,4021,4324,4503,4591,5223,5446,5641,5658,5836,5896,6391,6451,6811,6976,7021,7306,7471];
ont2 = [1398,2929,3559,3991,4998,6106,6722,7354,979,2026,3198,3634,4246,5657,6170,6903,7562,2417,3244,3917,454,5898,6455,7067,783,1067,272,3303,3946,4562,5927,6512,7203,888,1486,1876,2597,4668];

offt2hand = [4442 4459 4486 1098 4864 1310 4999 1384 5071 1581 514 1893 5433 1908 5464 1999 5569 2192 5851 2357 586 2462 6034 2522 6140 2746 6196 2809 6422 2868 6589 2881 6721 3061 6828 3139 7234 3226 7503 3258 7667 347 872 3601 94 3692 3813 4098];
offt1hand = [46 5223 5446 1051 5641 1081 5658 1381 5836 1637 5896 1652 6391 1966 6451 2177 6811 2401 6976 2686 7021 3287 7306 3589 7471 3933 751 4006 768 4324 782 4503 4591];
offt3 = [62 991 1156 4234 4278 4487 5733 6286 6931];
offt4 = [889 1684 3200 3679 3767 3935 4188 4294 5645 6064];
offt5 = [1246 2253 2373 3695 5116];
ont1hand = [154 2161 3319 4774 5853 7442 1786 2208 3422 5073 6542 7487 1877 2461 3691 5088 6826 7532 1111 1892 2851 3994 5359 692 860 1339 2011 3002 4 5567 7261];
ont2hand = [1398 2929 3559 3991 4998 6106 6722 7354 979 2026 3198 3634 4246 5657 6170 6903 7562 2417 3244 3917 454 5898 6455 7067 783 1067 272 3303 3946 4562 5927 6512 7203 888];
ont3hand = [1549 3452 4941 5405 6980 7520 1385 2136 4892 5179 6093 7157];





%2012-05-31-1 dataset
ds = [92 153 258 335 602 635 828 980 1247 1384 1412 1581 1639 1668 1848 1969 2252 2298 2704 2717 2854 2975 3019 3047 3181 3440 3693 3871 3933 3992 4324 4354 4564 4606 4714 5044 5057 5089 5221 5299 5584 5795 5988 6468 6530 6665 6889 6965];
nonds = [4 32 61 107 109 123 197 198 272 287 347 406 437 467 468 481 482 558 573 621 647 662 679 709 736 784 871 902 946 977 979 992 1006 1022 1041 1066 1174 1231 1246 1249 1261 1262 1264 1265 1295 1296 1307 1396 1397 1400 1426 1429 1459 1460 1501 1506 1563 1591 1606 1607 1608 1624 1625 1637 1741 1756 1774 1787 1789 1802 1817 1831 1834 1846 1996 2029 2041 2042 2071 2086 2210 2266 2267 2269 2326 2328 2356 2402 2403 2416 2446 2476 2671 2702 2746 2747 2750 2761 2763 2764 2808 2810 2838 2957 2960 2987 2989 3002 3005 3018 3031 3061 3091 3092 3122 3151 3184 3214 3229 3242 3256 3271 3273 3303 3363 3406 3453 3468 3470 3511 3516 3556 3601 3618 3619 3676 3722 3753 3841 3886 3961 3991 4038 4051 4066 4083 4096 4201 4231 4232 4246 4247 4276 4277 4338 4383 4429 4445 4486 4532 4546 4561 4562 4668 4697 4700 4713 4730 4833 4846 4864 4891 4906 4940 4951 4983 5000 5028 5086 5251 5252 5282 5326 5356 5390 5433 5448 5478 5524 5537 5585 5612 5627 5686 5732 5791 5807 5808 5836 5852 5882 5896 5971 5986 6032 6033 6121 6152 6226 6227 6271 6272 6301 6316 6318 6376 6409 6421 6436 6453 6481 6482 6557 6586 6631 6677 6680 6707 6752 6890 6901 6918 6946 6947 6949 6963 7006 7082 7097 7111 7156 7186 7231 7381 7441 7456 7458 7486 7533 7576 7608 7621 7636 7640 7652 7667];
on = [4 32 109 123 287 437 481 573 621 647 662 679 709 784 902 977 979 992 1006 1041 1174 1246 1262 1264 1265 1295 1396 1397 1400 1429 1459 1460 1591 1607 1624 1756 1774 1787 1789 1817 1834 1846 2029 2041 2086 2269 2326 2328 2403 2671 2747 2750 2763 2764 2808 2810 2838 2957 2960 3002 3005 3018 3092 3122 3151 3184 3214 3229 3273 3363 3406 3453 3511 3516 3618 3676 3841 3886 3961 3991 4201 4231 4232 4247 4338 4383 4429 4445 4486 4546 4562 4668 4697 4846 4940 4951 4983 5028 5086 5252 5282 5448 5478 5612 5732 5791 5807 5852 5896 6033 6226 6227 6272 6316 6318 6409 6421 6453 6482 6680 6752 6890 6901 6918 6947 6963 7097 7156 7441 7456 7458 7486 7533 7640 7652];
off = [61 107 197 198 272 347 406 467 468 482 558 736 871 946 1022 1066 1231 1249 1261 1296 1307 1426 1501 1506 1563 1606 1608 1625 1637 1741 1802 1831 1996 2042 2071 2210 2266 2267 2356 2402 2416 2446 2476 2702 2746 2761 2987 2989 3031 3061 3091 3242 3256 3271 3303 3468 3470 3556 3601 3619 3722 3753 4038 4051 4066 4083 4096 4246 4276 4277 4532 4561 4700 4713 4730 4833 4864 4891 4906 5000 5251 5326 5356 5390 5433 5524 5537 5585 5627 5686 5808 5836 5882 5971 5986 6032 6121 6152 6271 6301 6376 6436 6481 6557 6586 6631 6677 6707 6946 6949 7006 7082 7111 7186 7231 7381 7576 7608 7621 7636 7667];
ont1 = [4 902 1006 1246 1262 1396 1624 1787 2086 2269 2957 3184 3214 3406 4232 4940 5448 5612 6421 6918 6963 7441 7652];
ont2 = [123 287 437 481 573 679 992 1264 1459 1607 1756 2029 2328 2403 2808 3092 3122 3511 3841 3991 4247 4445 4846 5252 5732 5896 6272 6453 6752 6947 7097 7156 7456 7486 7533 7640];
%3 questionable cells in on t2:
%7640 6453 4445  
ont3 = [109 784 1174 1429 1817 1834 2810 3018 3229 3363];
offt1 = [197 406 468 736 871 946 1231 1426 1741 2042 2267 2446 2761 3061 3753 4246 4276 4906 5326 5356 5585 5686 5971 6032 6557 6631 6946 7111 7381 7636];
offt2 = [61 107 272 467 482 558 1022 1066 1261 1501 1608 1831 2071 2266 2416 2702 2746 3031 3256 3303 4051 4083 4561 4891 5251 5627 5808 5986 6121 6301 6586 6677 7186 7231]; %removed 1625
offt3 = [1506 4730 5537 6481 7082];
offt4 = [198 1249 3242 4038 4700 4833 5836 6152 6376];
offt5 = [4864 5524 7006]; 



%2012-10-31 dataset
ds = [92 619 1487 1996 2074 2613 3077 6242 7503 316 1037 1683 3125 5148 6602 484 559 1442 2118 2433 2942 3019 3260 4159 4548 4625 5000 5210 6290 6708 6781 7083 7159 7518];
cellids = [4,49,78,92,108,259,272,286,316,362,378,391,393,421,469,484,543,556,559,572,591,618,619,647,681,724,781,783,785,903,917,962,976,991,1006,1037,1067,1081,1083,1126,1172,1203,1248,1277,1306,1368,1381,1430,1442,1471,1487,1501,1576,1577,1595,1606,1670,1683,1712,1726,1731,1741,1756,1772,1773,1816,1862,1878,1922,1952,1954,1996,2028,2074,2086,2087,2101,2118,2146,2176,2177,2178,2356,2371,2419,2433,2494,2506,2536,2539,2581,2582,2597,2613,2656,2659,2687,2719,2747,2825,2851,2856,2867,2868,2884,2896,2897,2899,2942,2971,2973,3019,3046,3049,3076,3077,3121,3125,3215,3244,3260,3274,3289,3317,3482,3571,3616,3661,3843,3871,3887,3905,3931,3946,4036,4038,4066,4112,4113,4126,4142,4156,4159,4171,4172,4186,4204,4248,4277,4351,4366,4383,4384,4413,4518,4548,4578,4625,4681,4685,4697,4712,4713,4714,4727,4730,4771,4786,4846,4876,4892,4983,4999,5000,5026,5028,5042,5058,5072,5089,5117,5132,5148,5210,5240,5267,5281,5297,5328,5386,5401,5431,5462,5495,5506,5552,5596,5629,5642,5746,5748,5765,5791,5792,5794,5841,5853,5866,5897,5926,5943,5956,5961,6002,6003,6031,6033,6062,6091,6106,6122,6211,6212,6213,6242,6271,6289,6290,6376,6422,6483,6497,6499,6541,6542,6587,6590,6602,6617,6631,6646,6662,6692,6708,6752,6753,6781,6784,6796,6812,6827,6888,6901,6933,6962,6964,6993,7051,7067,7083,7098,7114,7156,7157,7159,7203,7291,7324,7336,7442,7444,7472,7502,7503,7518,7562,7593,7640,7668,7669];
nonds = [4,49,78,108,259,272,286,362,378,391,393,421,469,543,556,572,591,618,647,681,724,781,783,785,903,917,962,976,991,1006,1067,1081,1083,1126,1172,1203,1248,1277,1306,1368,1381,1430,1471,1501,1576,1577,1595,1606,1670,1712,1726,1731,1741,1756,1772,1773,1816,1862,1878,1922,1952,1954,2028,2086,2087,2101,2146,2176,2177,2178,2356,2371,2419,2494,2506,2536,2539,2581,2582,2597,2656,2659,2687,2719,2747,2825,2851,2856,2867,2868,2884,2896,2897,2899,2971,2973,3046,3049,3076,3121,3215,3244,3274,3289,3317,3482,3571,3616,3661,3843,3871,3887,3905,3931,3946,4036,4038,4066,4112,4113,4126,4142,4156,4171,4172,4186,4204,4248,4277,4351,4366,4383,4384,4413,4518,4578,4681,4685,4697,4712,4713,4714,4727,4730,4771,4786,4846,4876,4892,4983,4999,5026,5028,5042,5058,5072,5089,5117,5132,5240,5267,5281,5297,5328,5386,5401,5431,5462,5495,5506,5552,5596,5629,5642,5746,5748,5765,5791,5792,5794,5841,5853,5866,5897,5926,5943,5956,5961,6002,6003,6031,6033,6062,6091,6106,6122,6211,6212,6213,6271,6289,6376,6422,6483,6497,6499,6541,6542,6587,6590,6617,6631,6646,6662,6692,6752,6753,6784,6796,6812,6827,6888,6901,6933,6962,6964,6993,7051,7067,7098,7114,7156,7157,7203,7291,7324,7336,7442,7444,7472,7502,7562,7593,7640,7668,7669];
on = [4,78,108,286,362,543,618,647,781,785,903,917,962,991,1067,1081,1203,1248,1306,1381,1471,1595,1606,1712,1756,1878,1952,2086,2087,2146,2371,2419,2494,2536,2582,2687,2719,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4038,4113,4156,4172,4186,4248,4277,4384,4685,4697,4713,4714,4876,4892,4983,5028,5042,5132,5240,5267,5297,5328,5401,5552,5642,5746,5765,5792,5794,5866,5897,5956,5961,6002,6031,6211,6289,6497,6590,6646,6662,6692,6752,6753,6784,6812,6962,7067,7098,7114,7336,7442,7472,7640,7669];
off = [49,259,272,378,391,393,421,469,556,572,591,681,724,783,976,1006,1083,1126,1172,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2101,2176,2177,2178,2356,2506,2539,2581,2597,2656,2659,2825,2867,2868,2884,2899,2971,2973,3046,3049,3121,3274,3317,3482,3571,3616,3661,3871,3905,3946,4036,4066,4112,4126,4142,4171,4204,4351,4366,4383,4413,4518,4578,4681,4712,4727,4730,4771,4786,4846,4999,5026,5058,5072,5089,5117,5281,5386,5431,5462,5495,5506,5596,5629,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6122,6212,6213,6271,6376,6422,6483,6499,6541,6542,6587,6617,6631,6796,6827,6888,6901,6933,6964,6993,7051,7156,7157,7203,7291,7324,7444,7502,7562,7593,7668];

onwithsnrcutoff = [4,78,108,286,362,543,618,647,781,785,903,917,962,991,1067,1203,1248,1306,1381,1471,1595,1606,1712,1756,1878,1952,2086,2087,2146,2371,2419,2494,2536,2582,2687,2719,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4113,4156,4172,4186,4248,4277,4384,4685,4697,4713,4714,4876,4892,4983,5028,5042,5132,5240,5267,5297,5328,5401,5552,5642,5746,5765,5792,5794,5866,5897,5956,5961,6002,6031,6211,6289,6497,6590,6646,6662,6692,6752,6753,6784,6812,6962,7067,7098,7114,7336,7442,7472,7640,7669];
offwithsnrcutoff = [49,259,272,378,391,393,421,469,556,572,591,681,724,783,976,1006,1083,1126,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2101,2176,2177,2178,2356,2506,2539,2581,2597,2656,2659,2825,2867,2868,2884,2899,2971,2973,3046,3049,3274,3317,3482,3571,3616,3661,3871,3905,3946,4036,4066,4112,4142,4171,4204,4351,4366,4383,4413,4518,4578,4681,4712,4727,4730,4771,4786,4846,4999,5026,5058,5072,5089,5117,5281,5386,5431,5462,5495,5506,5596,5629,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6122,6212,6213,6271,6376,6422,6483,6499,6541,6542,6587,6617,6631,6796,6827,6888,6901,6933,6964,6993,7051,7156,7157,7203,7291,7324,7444,7502,7562,7593,7668];

ont11ststage = [4,78,286,618,917,1248,1606,1712,2087,2719,4113,4156,4172,4248,4277,4697,4713,4786,4983,5401,5552,5897,6289,6590,6888,7067,7472];
ont1getridfrom1ststage = [1606 4113 4156 4786 6888 ];
%84 cells=[108,362,543,647,781,785,903,962,991,1067,1081,1172,1203,1306,1381,1471,1595,1756,1878,1952,2086,2146,2371,2419,2494,2536,2582,2687,2747,2851,2856,2897,3076,3215,3244,3289,3843,3887,3931,4038,4186,4384,4685,4714,4876,4892,5028,5042,5132,5240,5267,5297,5328,5642,5746,5792,5794,5866,5956,5961,6002,6031,6211,6497,6646,6662,6692,6752,6753,6784,6812,6962,7098,7114,7336,7442,7444,7640,7669,1606,4113,4156,4786,6888];
%50 cells =[543,781,785,1067,1081,1172,1306,1381,1471,1595,1878,1952,2146,2494,2536,2687,2747,2851,3215,3244,3931,4038,4186,4714,4892,5042,5132,5240,5267,5328,5792,5794,5866,5961,6211,6646,6662,6692,6753,6784,6962,7114,7336,7444,7640,1606,4113,4156,4786,6888];
ont21ststge= [108,362,647,903,962,991,1081,1203,1756,2086,2371,2419,2582,2856,2897,3076,3289,3843,3887,3931,4384,4685,4876,5028,5267,5297,5642,5746,5956,6002,6031,6497,6662,6752,6812,6962,7098,7442,7444,7669,4113]; %41 cells
ont2 = [108 362 647 903 962 991 1203 1756 2086 2371 2419 2582 2856 2897 3076 3289 3843 3887 4384 4685 4876 5028 5297 5642 5746 5956 6002 6031 6497 6752 6812 7098 7442 7669]; 
nonds = [4,49,78,108,259,272,286,362,378,391,393,421,469,543,556,572,591,618,647,681,724,781,783,785,903,917,962,976,991,1006,1067,1081,1083,1126,1172,1203,1248,1277,1306,1368,1381,1430,1471,1501,1576,1577,1595,1606,1670,1712,1726,1731,1741,1756,1772,1773,1816,1862,1878,1922,1952,1954,2028,2086,2087,2101,2146,2176,2177,2178,2356,2371,2419,2494,2506,2536,2539,2581,2582,2597,2656,2659,2687,2719,2747,2825,2851,2856,2867,2868,2884,2896,2897,2899,2971,2973,3046,3049,3076,3121,3215,3244,3274,3289,3317,3482,3571,3616,3661,3843,3871,3887,3905,3931,3946,4036,4038,4066,4112,4113,4126,4142,4156,4171,4172,4186,4204,4248,4277,4351,4366,4383,4384,4413,4518,4578,4681,4685,4697,4712,4713,4714,4727,4730,4771,4786,4846,4876,4892,4983,4999,5026,5028,5042,5058,5072,5089,5117,5132,5240,5267,5281,5297,5328,5386,5401,5431,5462,5495,5506,5552,5596,5629,5642,5746,5748,5765,5791,5792,5794,5841,5853,5866,5897,5926,5943,5956,5961,6002,6003,6031,6033,6062,6091,6106,6122,6211,6212,6213,6271,6289,6376,6422,6483,6497,6499,6541,6542,6587,6590,6617,6631,6646,6662,6692,6752,6753,6784,6796,6812,6827,6888,6901,6933,6962,6964,6993,7051,7067,7098,7114,7156,7157,7203,7291,7324,7336,7442,7444,7472,7502,7562,7593,7640,7668,7669];

ont3hand = [1067 1306 1595 1878 1952 2494 2747 2851 5240 5328 6692 6753];
ont2hand = [108 362 647 903 962 991 1203 1756 2419 2582 2856 2897 3076 3289 3843 3887 3931 4384 4685 4876 5028 5297 5642 5746 5956 6002 6031 6211 6497 6752 6812 7098 7442]; 

offt1 = [offt1hand 2506 4066];
offt11ststage = [49,378,572,1083,1126,1501,1741,2101,2176,2506,2825,3049,3317,3616,3661,4066,4112,4351,4413,4518,4578,4846,5058,5462,5596,5629,6062,6106,6122,6213,6422,6483,6541,6617,6933,6993,7051,7156,7668];
offt1hand = [378 572 1083 1126 2101 2176 2825 3049 3317 3616 4112 4351 4413 4518 4846 5058 5462 5629 6122 6422 6617 6933 7051 7156];
%115 cells = [259,272,391,393,421,469,556,591,681,724,783,976,1006,1277,1368,1430,1576,1577,1670,1726,1731,1772,1773,1816,1862,1922,1954,2028,2177,2178,2356,2539,2581,2597,2656,2659,2867,2868,2884,2896,2899,2971,2973,3046,3121,3274,3482,3571,3871,3905,3946,4036,4126,4142,4171,4204,4366,4383,4681,4712,4727,4730,4771,4999,5026,5072,5089,5117,5281,5386,5431,5495,5506,5748,5765,5791,5841,5853,5926,5943,6003,6033,6091,6212,6271,6376,6499,6542,6587,6631,6796,6827,6901,6964,7157,7203,7291,7324,7502,7562,7593,49,1501,1741,3661,4578,5596,6062,6106,6213,6483,6541,6993,7668,4518];
offt2hand = [421 469 556 976 1006 1277 1576 1670 1726 1741 1922 2028 2356 2581 2656 2868 3046 3274 3482 3571 3871 4142 4171 4366 4578 4712 4999 5026 5281 5386 5748 5791 5943 6003 6106 6271 6587 6631 6827 6964 7291 7502 5853];


%2012-10-10-1
cellids = [17,31,77,91,211,242,244,272,317,362,395,437,454,481,496,512,558,571,631,646,647,677,695,706,783,796,800,813,917,932,995,1023,1036,1081,1100,1126,1142,1144,1189,1231,1246,1263,1265,1279,1281,1368,1369,1396,1397,1416,1442,1459,1486,1517,1518,1531,1594,1595,1624,1666,1681,1696,1697,1701,1712,1726,1759,1801,1803,1863,1876,1879,1894,1895,1938,1939,1983,2029,2057,2073,2117,2118,2134,2146,2177,2193,2401,2402,2435,2461,2462,2464,2491,2507,2536,2555,2597,2631,2701,2703,2716,2748,2750,2794,2836,2975,3002,3016,3017,3031,3076,3093,3121,3182,3213,3287,3306,3316,3331,3363,3421,3422,3423,3452,3470,3496,3512,3515,3574,3587,3617,3618,3736,3766,3769,3811,3828,3856,3874,3931,3947,4024,4038,4055,4068,4126,4127,4128,4174,4188,4261,4262,4307,4321,4323,4336,4366,4367,4384,4443,4471,4487,4488,4519,4534,4548,4549,4577,4580,4606,4621,4711,4712,4728,4771,4833,4877,4953,4981,4983,4996,5011,5087,5131,5146,5176,5328,5358,5375,5403,5419,5449,5465,5496,5506,5508,5536,5614,5627,5642,5671,5702,5732,5746,5761,5777,5821,5836,5868,5957,5959,6031,6047,6110,6136,6138,6170,6186,6196,6197,6226,6241,6272,6275,6319,6331,6378,6391,6438,6483,6511,6544,6556,6572,6588,6617,6633,6676,6722,6751,6752,6754,6767,6796,6813,6842,6856,6932,6949,6977,6992,7007,7052,7096,7114,7189,7232,7276,7336,7352,7353,7354,7381,7382,7426,7472,7561,7564,7608,7610,7621,7668];
ds = [395,695,706,995,1023,1265,1369,1416,1595,1666,1681,1712,1879,2555,2631,2975,3306,3470,3574,3769,3828,3931,4323,4366,4549,4577,4833,5375,5419,5465,5627,6110,6138,6186,6544,6633,6754,7114,7381];
nonds = [17,31,77,91,211,242,244,272,317,362,437,454,481,496,512,558,571,631,646,647,677,783,796,800,813,917,932,1036,1081,1100,1126,1142,1144,1189,1231,1246,1263,1279,1281,1368,1396,1397,1442,1459,1486,1517,1518,1531,1594,1624,1696,1697,1701,1726,1759,1801,1803,1863,1876,1894,1895,1938,1939,1983,2029,2057,2073,2117,2118,2134,2146,2177,2193,2401,2402,2435,2461,2462,2464,2491,2507,2536,2597,2701,2703,2716,2748,2750,2794,2836,3002,3016,3017,3031,3076,3093,3121,3182,3213,3287,3316,3331,3363,3421,3422,3423,3452,3496,3512,3515,3587,3617,3618,3736,3766,3811,3856,3874,3947,4024,4038,4055,4068,4126,4127,4128,4174,4188,4261,4262,4307,4321,4336,4367,4384,4443,4471,4487,4488,4519,4534,4548,4580,4606,4621,4711,4712,4728,4771,4877,4953,4981,4983,4996,5011,5087,5131,5146,5176,5328,5358,5403,5449,5496,5506,5508,5536,5614,5642,5671,5702,5732,5746,5761,5777,5821,5836,5868,5957,5959,6031,6047,6136,6170,6196,6197,6226,6241,6272,6275,6319,6331,6378,6391,6438,6483,6511,6556,6572,6588,6617,6676,6722,6751,6752,6767,6796,6813,6842,6856,6932,6949,6977,6992,7007,7052,7096,7189,7232,7276,7336,7352,7353,7354,7382,7426,7472,7561,7564,7608,7610,7621,7668];
on = [31,77,211,242,244,272,558,631,646,677,783,813,1036,1081,1100,1144,1189,1281,1397,1442,1517,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2402,2435,2464,2491,2703,2716,2748,3002,3031,3213,3287,3316,3331,3422,3512,3515,3587,3618,3811,4024,4055,4126,4127,4174,4262,4443,4488,4534,4548,4711,4712,4877,5087,5146,5176,5496,5506,5642,5671,5702,5777,5821,5868,5957,6136,6197,6272,6275,6319,6572,6722,6751,6752,6767,6796,6977,7096,7189,7232,7354,7426,7472,7564,7610,7621,7668];
off = [17,91,317,362,437,454,481,496,512,571,647,796,800,917,932,1126,1142,1231,1246,1263,1279,1368,1396,1459,1486,1518,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2146,2193,2461,2462,2507,2536,2597,2701,2750,2794,2836,3016,3017,3076,3093,3121,3182,3363,3421,3423,3452,3496,3617,3736,3766,3856,3874,3947,4038,4068,4128,4188,4261,4307,4321,4336,4367,4384,4471,4487,4519,4580,4606,4621,4728,4771,4953,4981,4983,4996,5011,5131,5328,5358,5403,5449,5508,5536,5614,5732,5746,5761,5836,5959,6031,6047,6170,6196,6226,6241,6331,6378,6391,6438,6483,6511,6556,6588,6617,6676,6813,6842,6856,6932,6949,6992,7007,7052,7276,7336,7352,7353,7382,7561,7608];
onwithsnrcutoff = [31,77,211,242,244,272,558,631,646,677,813,1036,1081,1100,1144,1189,1281,1397,1442,1517,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2402,2435,2464,2491,2703,2716,2748,3002,3031,3213,3287,3316,3331,3422,3512,3515,3587,3618,3811,4024,4126,4127,4174,4262,4443,4488,4534,4548,4711,4712,4877,5087,5146,5176,5496,5506,5642,5671,5702,5777,5821,5868,5957,6136,6197,6272,6275,6572,6722,6751,6752,6767,6796,6977,7096,7189,7232,7426,7472,7564,7610,7621,7668];
offwithsnrcutoff = [17,91,317,362,437,454,481,496,512,571,647,796,800,917,932,1126,1142,1231,1246,1263,1279,1368,1396,1459,1486,1518,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2146,2193,2461,2462,2507,2536,2597,2701,2750,2794,2836,3016,3017,3076,3093,3121,3182,3363,3421,3423,3452,3496,3617,3736,3766,3856,3874,3947,4038,4068,4128,4188,4261,4307,4321,4336,4367,4384,4471,4487,4519,4580,4606,4621,4728,4771,4953,4981,4983,4996,5011,5131,5328,5358,5403,5449,5508,5536,5614,5732,5746,5761,5836,5959,6031,6047,6170,6196,6226,6241,6378,6391,6438,6483,6511,6556,6588,6617,6676,6813,6842,6856,6932,6949,6992,7007,7052,7276,7336,7352,7353,7382,7561,7608];
ont2 = [31 272  646  1036 1081 1397 1442 1696 2134  2177 2401 2435 2491 2703 3002 3031 3316 3512 3587 3811 4174 4262 4711 5146 5702 5957 6136 6272 6572 6722 6751 6796 7426 7472];



















cellids_10 = [17,31,77,91,211,242,244,272,317,362,395,437,454,481,496,512,558,571,631,646,647,677,695,706,783,796,800,813,917,932,995,1023,1036,1081,1100,1126,1142,1144,1189,1231,1246,1263,1265,1279,1281,1368,1369,1396,1397,1416,1442,1459,1486,1517,1518,1531,1594,1595,1624,1666,1681,1696,1697,1701,1712,1726,1759,1801,1803,1863,1876,1879,1894,1895,1938,1939,1983,2029,2057,2073,2117,2118,2134,2146,2177,2193,2401,2402,2435,2461,2462,2464,2491,2507,2536,2555,2597,2631,2701,2703,2716,2748,2750,2794,2836,2975,3002,3016,3017,3031,3076,3093,3121,3182,3213,3287,3306,3316,3331,3363,3421,3422,3423,3452,3470,3496,3512,3515,3574,3587,3617,3618,3736,3766,3769,3811,3828,3856,3874,3931,3947,4024,4038,4055,4068,4126,4127,4128,4174,4188,4261,4262,4307,4321,4323,4336,4366,4367,4384,4443,4471,4487,4488,4519,4534,4548,4549,4577,4580,4606,4621,4711,4712,4728,4771,4833,4877,4953,4981,4983,4996,5011,5087,5131,5146,5176,5328,5358,5375,5403,5419,5449,5465,5496,5506,5508,5536,5614,5627,5642,5671,5702,5732,5746,5761,5777,5821,5836,5868,5957,5959,6031,6047,6110,6136,6138,6170,6186,6196,6197,6226,6241,6272,6275,6319,6331,6378,6391,6438,6483,6511,6544,6556,6572,6588,6617,6633,6676,6722,6751,6752,6754,6767,6796,6813,6842,6856,6932,6949,6977,6992,7007,7052,7096,7114,7189,7232,7276,7336,7352,7353,7354,7381,7382,7426,7472,7561,7564,7608,7610,7621,7668]
cellids_15 = [4,31,46,62,76,94,154,182,226,257,272,301,333,347,407,424,438,454,467,496,514,528,586,649,692,708,751,753,766,768,782,783,843,857,860,872,888,889,901,979,991,1051,1067,1081,1084,1097,1098,1111,1127,1130,1156,1191,1246,1310,1339,1381,1384,1385,1398,1411,1427,1486,1532,1549,1564,1578,1581,1591,1595,1637,1652,1683,1684,1685,1786,1817,1876,1877,1892,1893,1895,1908,1921,1966,1969,1999,2011,2026,2042,2136,2161,2177,2192,2206,2208,2236,2253,2311,2326,2328,2343,2357,2373,2401,2417,2449,2461,2462,2478,2521,2522,2555,2597,2686,2716,2732,2746,2794,2809,2851,2868,2881,2896,2898,2929,3002,3061,3091,3139,3152,3198,3200,3226,3241,3244,3258,3286,3287,3303,3319,3422,3452,3512,3530,3559,3586,3589,3601,3634,3635,3636,3648,3679,3691,3692,3695,3721,3736,3767,3812,3813,3815,3842,3857,3859,3889,3917,3933,3934,3935,3946,3991,3994,4006,4021,4022,4069,4096,4097,4098,4130,4145,4157,4173,4188,4231,4234,4235,4246,4278,4279,4294,4324,4353,4427,4442,4459,4486,4487,4501,4503,4562,4591,4668,4697,4731,4732,4771,4774,4788,4789,4846,4864,4892,4941,4985,4997,4998,4999,5071,5073,5088,5104,5116,5148,5150,5179,5223,5359,5405,5433,5446,5464,5567,5569,5632,5641,5645,5657,5658,5672,5702,5703,5705,5719,5733,5791,5836,5851,5853,5896,5898,5927,5941,6034,6064,6093,6106,6125,6139,6140,6152,6155,6170,6196,6229,6257,6260,6286,6304,6321,6332,6361,6363,6376,6380,6391,6392,6422,6439,6451,6455,6512,6542,6589,6721,6722,6737,6751,6752,6797,6811,6826,6828,6886,6903,6931,6976,6980,6992,7021,7040,7067,7069,7096,7157,7186,7203,7234,7261,7278,7306,7308,7354,7442,7471,7475,7487,7503,7517,7520,7532,7562,7667];
cellids_31 = [4,49,78,92,108,259,272,286,316,362,378,391,393,421,469,484,543,556,559,572,591,618,619,647,681,724,781,783,785,903,917,962,976,991,1006,1037,1067,1081,1083,1126,1172,1203,1248,1277,1306,1368,1381,1430,1442,1471,1487,1501,1576,1577,1595,1606,1670,1683,1712,1726,1731,1741,1756,1772,1773,1816,1862,1878,1922,1952,1954,1996,2028,2074,2086,2087,2101,2118,2146,2176,2177,2178,2356,2371,2419,2433,2494,2506,2536,2539,2581,2582,2597,2613,2656,2659,2687,2719,2747,2825,2851,2856,2867,2868,2884,2896,2897,2899,2942,2971,2973,3019,3046,3049,3076,3077,3121,3125,3215,3244,3260,3274,3289,3317,3482,3571,3616,3661,3843,3871,3887,3905,3931,3946,4036,4038,4066,4112,4113,4126,4142,4156,4159,4171,4172,4186,4204,4248,4277,4351,4366,4383,4384,4413,4518,4548,4578,4625,4681,4685,4697,4712,4713,4714,4727,4730,4771,4786,4846,4876,4892,4983,4999,5000,5026,5028,5042,5058,5072,5089,5117,5132,5148,5210,5240,5267,5281,5297,5328,5386,5401,5431,5462,5495,5506,5552,5596,5629,5642,5746,5748,5765,5791,5792,5794,5841,5853,5866,5897,5926,5943,5956,5961,6002,6003,6031,6033,6062,6091,6106,6122,6211,6212,6213,6242,6271,6289,6290,6376,6422,6483,6497,6499,6541,6542,6587,6590,6602,6617,6631,6646,6662,6692,6708,6752,6753,6781,6784,6796,6812,6827,6888,6901,6933,6962,6964,6993,7051,7067,7083,7098,7114,7156,7157,7159,7203,7291,7324,7336,7442,7444,7472,7502,7503,7518,7562,7593,7640,7668,7669];
ds_10 = [395,695,706,995,1023,1265,1369,1416,1595,1666,1681,1712,1879,2555,2631,2975,3306,3470,3574,3769,3828,3931,4323,4366,4549,4577,4833,5375,5419,5465,5627,6110,6138,6186,6544,6633,6754,7114,7381];
ds_31 = [92,316,484,559,619,1037,1442,1487,1683,1996,2074,2118,2433,2613,2942,3019,3077,3125,3260,4159,4548,4625,5000,5148,5210,6242,6290,6602,6708,6781,7083,7159,7503,7518];
nonds_10 = [17,31,77,91,211,242,244,272,317,362,437,454,481,496,512,558,571,631,646,647,677,783,796,800,813,917,932,1036,1081,1100,1126,1142,1144,1189,1231,1246,1263,1279,1281,1368,1396,1397,1442,1459,1486,1517,1518,1531,1594,1624,1696,1697,1701,1726,1759,1801,1803,1863,1876,1894,1895,1938,1939,1983,2029,2057,2073,2117,2118,2134,2146,2177,2193,2401,2402,2435,2461,2462,2464,2491,2507,2536,2597,2701,2703,2716,2748,2750,2794,2836,3002,3016,3017,3031,3076,3093,3121,3182,3213,3287,3316,3331,3363,3421,3422,3423,3452,3496,3512,3515,3587,3617,3618,3736,3766,3811,3856,3874,3947,4024,4038,4055,4068,4126,4127,4128,4174,4188,4261,4262,4307,4321,4336,4367,4384,4443,4471,4487,4488,4519,4534,4548,4580,4606,4621,4711,4712,4728,4771,4877,4953,4981,4983,4996,5011,5087,5131,5146,5176,5328,5358,5403,5449,5496,5506,5508,5536,5614,5642,5671,5702,5732,5746,5761,5777,5821,5836,5868,5957,5959,6031,6047,6136,6170,6196,6197,6226,6241,6272,6275,6319,6331,6378,6391,6438,6483,6511,6556,6572,6588,6617,6676,6722,6751,6752,6767,6796,6813,6842,6856,6932,6949,6977,6992,7007,7052,7096,7189,7232,7276,7336,7352,7353,7354,7382,7426,7472,7561,7564,7608,7610,7621,7668];
nonds_15 = [4,31,46,62,76,94,154,182,226,272,347,407,424,454,496,514,528,586,649,692,751,753,766,768,782,783,843,857,860,872,888,889,901,979,991,1051,1067,1081,1084,1098,1111,1127,1130,1156,1191,1246,1310,1339,1381,1384,1385,1398,1411,1427,1486,1532,1549,1564,1578,1581,1591,1637,1652,1684,1786,1817,1876,1877,1892,1893,1908,1921,1966,1969,1999,2011,2026,2136,2161,2177,2192,2206,2208,2236,2253,2311,2326,2328,2343,2357,2373,2401,2417,2461,2462,2478,2521,2522,2555,2597,2686,2716,2732,2746,2794,2809,2851,2868,2881,2896,2929,3002,3061,3091,3139,3152,3198,3200,3226,3241,3244,3258,3286,3287,3303,3319,3422,3452,3530,3559,3586,3589,3601,3634,3635,3648,3679,3691,3692,3695,3721,3767,3812,3813,3857,3859,3889,3917,3933,3935,3946,3991,3994,4006,4021,4022,4096,4097,4098,4130,4145,4188,4231,4234,4235,4246,4278,4279,4294,4324,4427,4442,4459,4486,4487,4501,4503,4562,4591,4668,4697,4732,4771,4774,4788,4864,4892,4941,4997,4998,4999,5071,5073,5088,5104,5116,5148,5179,5223,5359,5405,5433,5446,5464,5567,5569,5641,5645,5657,5658,5672,5703,5705,5733,5791,5836,5851,5853,5896,5898,5927,5941,6034,6064,6093,6106,6125,6139,6140,6152,6155,6170,6196,6257,6260,6286,6304,6361,6363,6376,6380,6391,6392,6422,6439,6451,6455,6512,6542,6589,6721,6722,6737,6752,6811,6826,6828,6886,6903,6931,6976,6980,6992,7021,7040,7067,7069,7096,7157,7186,7203,7234,7261,7278,7306,7354,7442,7471,7475,7487,7503,7517,7520,7532,7562,7667];
nonds_31 = [4,49,78,108,259,272,286,362,378,391,393,421,469,543,556,572,591,618,647,681,724,781,783,785,903,917,962,976,991,1006,1067,1081,1083,1126,1172,1203,1248,1277,1306,1368,1381,1430,1471,1501,1576,1577,1595,1606,1670,1712,1726,1731,1741,1756,1772,1773,1816,1862,1878,1922,1952,1954,2028,2086,2087,2101,2146,2176,2177,2178,2356,2371,2419,2494,2506,2536,2539,2581,2582,2597,2656,2659,2687,2719,2747,2825,2851,2856,2867,2868,2884,2896,2897,2899,2971,2973,3046,3049,3076,3121,3215,3244,3274,3289,3317,3482,3571,3616,3661,3843,3871,3887,3905,3931,3946,4036,4038,4066,4112,4113,4126,4142,4156,4171,4172,4186,4204,4248,4277,4351,4366,4383,4384,4413,4518,4578,4681,4685,4697,4712,4713,4714,4727,4730,4771,4786,4846,4876,4892,4983,4999,5026,5028,5042,5058,5072,5089,5117,5132,5240,5267,5281,5297,5328,5386,5401,5431,5462,5495,5506,5552,5596,5629,5642,5746,5748,5765,5791,5792,5794,5841,5853,5866,5897,5926,5943,5956,5961,6002,6003,6031,6033,6062,6091,6106,6122,6211,6212,6213,6271,6289,6376,6422,6483,6497,6499,6541,6542,6587,6590,6617,6631,6646,6662,6692,6752,6753,6784,6796,6812,6827,6888,6901,6933,6962,6964,6993,7051,7067,7098,7114,7156,7157,7203,7291,7324,7336,7442,7444,7472,7502,7562,7593,7640,7668,7669];
off_10 = [17,91,317,362,437,454,481,496,512,571,647,796,800,917,932,1126,1142,1231,1246,1263,1279,1368,1396,1459,1486,1518,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2146,2193,2461,2462,2507,2536,2597,2701,2750,2794,2836,3016,3017,3076,3093,3121,3182,3363,3421,3423,3452,3496,3617,3736,3766,3856,3874,3947,4038,4068,4128,4188,4261,4307,4321,4336,4367,4384,4471,4487,4519,4580,4606,4621,4728,4771,4953,4981,4983,4996,5011,5131,5328,5358,5403,5449,5508,5536,5614,5732,5746,5761,5836,5959,6031,6047,6170,6196,6226,6241,6378,6391,6438,6483,6511,6556,6588,6617,6676,6813,6842,6856,6932,6949,6992,7007,7052,7276,7336,7352,7353,7382,7561,7608];
off_15 = [46,62,76,94,226,347,407,514,586,649,751,753,768,782,843,872,889,901,991,1051,1081,1084,1098,1127,1130,1156,1191,1246,1310,1381,1384,1411,1427,1564,1581,1637,1652,1684,1817,1893,1908,1921,1966,1999,2177,2192,2206,2236,2253,2326,2343,2357,2373,2401,2462,2478,2521,2522,2555,2686,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3287,3589,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3933,3935,4006,4021,4097,4098,4188,4234,4235,4278,4294,4324,4442,4459,4486,4487,4503,4591,4697,4732,4771,4788,4864,4999,5071,5116,5148,5223,5433,5446,5464,5569,5641,5645,5658,5672,5703,5705,5733,5836,5851,5896,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6391,6422,6451,6589,6721,6811,6828,6886,6931,6976,6992,7021,7040,7069,7234,7278,7306,7471,7503,7517,7667];
off_31 = [49,259,272,378,391,393,421,469,556,572,591,681,724,783,976,1006,1083,1126,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2101,2176,2177,2178,2356,2506,2539,2581,2597,2656,2659,2825,2867,2868,2884,2899,2971,2973,3046,3049,3274,3317,3482,3571,3616,3661,3871,3905,3946,4036,4066,4112,4142,4171,4204,4351,4366,4383,4413,4518,4578,4681,4712,4727,4730,4771,4786,4846,4999,5026,5058,5072,5089,5117,5281,5386,5431,5462,5495,5506,5596,5629,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6122,6212,6213,6271,6376,6422,6483,6499,6541,6542,6587,6617,6631,6796,6827,6888,6901,6933,6964,6993,7051,7156,7157,7203,7291,7324,7444,7502,7562,7593,7668];
off_other10 = [17,91,362,437,454,496,512,647,796,800,932,1231,1246,1263,1279,1368,1396,1459,1486,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2193,2462,2507,2536,2597,2701,2750,2794,3016,3076,3093,3121,3182,3363,3423,3452,3496,3617,3766,3874,3947,4038,4068,4188,4307,4321,4367,4384,4471,4487,4519,4580,4606,4621,4953,4981,4983,4996,5011,5328,5358,5449,5536,5614,5732,5746,5836,5959,6031,6170,6226,6241,6378,6391,6483,6511,6556,6588,6617,6842,6856,6932,6949,7007,7276,7336,7352,7353,7561,7608];
off_other15 = [62,76,94,226,347,407,514,586,649,753,843,872,889,901,991,1084,1098,1127,1130,1156,1191,1246,1310,1384,1411,1427,1564,1581,1684,1817,1893,1908,1921,1999,2192,2206,2236,2253,2326,2343,2357,2373,2462,2478,2521,2522,2555,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3935,4097,4098,4188,4234,4235,4278,4294,4442,4459,4486,4487,4697,4732,4771,4788,4864,4999,5071,5116,5148,5433,5464,5569,5645,5672,5703,5705,5733,5851,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6422,6589,6721,6828,6886,6931,6992,7040,7069,7234,7278,7503,7517,7667];
off_other31 = [49,259,272,391,393,421,469,556,591,681,724,783,976,1006,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2177,2178,2356,2539,2581,2597,2656,2659,2867,2868,2884,2899,2971,2973,3046,3274,3482,3571,3661,3871,3905,3946,4036,4066,4142,4171,4204,4366,4383,4518,4578,4681,4712,4727,4730,4771,4786,4999,5026,5072,5089,5117,5281,5386,5431,5495,5506,5596,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6212,6213,6271,6376,6483,6499,6541,6542,6587,6631,6796,6827,6888,6901,6964,6993,7157,7203,7291,7324,7444,7502,7562,7593];
off_otherother10 = [17,454,496,512,800,1231,1279,1368,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2117,2462,2507,2597,2750,2794,3016,3076,3093,3363,3423,3617,3874,4068,4188,4321,4367,4384,4487,4519,4580,4606,4983,4996,5011,5449,5536,5614,5732,5746,5836,5959,6170,6241,6378,6483,6511,6588,6617,6856,6949,7336,7353,7561,7608];
off_otherother15 = [62,226,407,649,843,889,901,991,1084,1130,1156,1191,1246,1411,1427,1564,1684,1817,2206,2236,2253,2326,2343,2373,2478,2521,2555,2732,2794,2896,3091,3200,3286,3635,3648,3679,3695,3767,3857,3889,3935,4097,4188,4234,4235,4278,4294,4487,4697,4732,4771,4788,5116,5148,5645,5672,5705,5733,6064,6139,6257,6260,6286,6361,6376,6380,6886,6931,6992,7040,7069,7278,7517];
off_otherother31 = [49,259,272,391,393,591,681,724,783,1368,1430,1501,1577,1731,1741,1772,1773,1816,1862,1954,2177,2178,2539,2597,2659,2867,2884,2899,2971,2973,3661,3905,3946,4036,4066,4204,4383,4518,4578,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5853,5926,5943,6033,6062,6091,6106,6212,6213,6376,6483,6499,6541,6542,6796,6888,6901,6993,7157,7203,7324,7444,7562,7593];
off_otherotherother10 = [17,454,496,800,1231,1279,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2117,2462,2507,2597,2794,3016,3076,3093,3363,3423,3617,3874,4068,4321,4367,4384,4519,4580,4606,4983,4996,5011,5449,5536,5614,5746,5836,5959,6241,6378,6483,6617,6856,6949,7336,7353,7561,7608];
off_otherotherother15 = [226,407,649,843,889,901,1084,1130,1191,1246,1411,1427,1564,1684,1817,2206,2236,2253,2326,2343,2373,2478,2521,2555,2732,2794,2896,3091,3200,3286,3635,3648,3679,3695,3767,3857,3889,3935,4097,4188,4235,4294,4697,4732,4771,4788,5116,5148,5645,5672,5705,6064,6139,6257,6260,6361,6376,6380,6886,6992,7040,7069,7278,7517];
off_otherotherother31 = [49,259,272,391,393,591,681,724,783,1368,1430,1501,1731,1741,1772,1773,1816,1862,2177,2178,2539,2597,2867,2899,2971,2973,3661,3905,3946,4036,4066,4204,4518,4578,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5853,5926,5943,6062,6091,6106,6212,6213,6376,6483,6499,6541,6542,6796,6888,6901,6993,7157,7203,7444,7562,7593];

off_otherotherotherother10 = [17,454,496,800,1231,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2462,2597,2794,3093,3363,3423,3617,3874,4321,4367,4384,4519,4580,4606,4996,5449,5536,5614,5746,5836,5959,6378,6483,6856,6949,7336,7353,7561];
off_otherotherotherother15 = [226,407,649,843,901,1130,1246,1411,1427,1564,1817,2206,2236,2253,2326,2343,2373,2521,2555,2732,2794,3091,3286,3695,3857,3889,4097,4235,4697,4732,4771,5116,5148,5672,6139,6257,6260,6376,6380,6886,6992,7517];
off_otherotherotherother31 = [49,259,272,391,393,724,783,1368,1430,1501,1731,1741,1773,1816,1862,2177,2539,2597,2867,2971,3661,3905,4066,4204,4518,4578,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5853,5926,5943,6062,6091,6106,6212,6213,6376,6483,6541,6542,6796,6888,6901,6993,7157,7203,7562,7593];


offt1_10 = [317,481,571,917,1126,1142,1518,2146,2461,2836,3017,3421,3736,3856,4128,4261,4336,4728,4771,5131,5403,5508,5761,6047,6196,6438,6676,6813,6992,7052,7382];
offt1_15 = [46,751,768,782,1051,1081,1381,1637,1652,1966,2177,2401,2686,3287,3589,3933,4006,4021,4324,4503,4591,5223,5446,5641,5658,5836,5896,6391,6451,6811,6976,7021,7306,7471];
offt1_31 = [378,572,1083,1126,2101,2176,2506,2825,3049,3317,3616,4112,4351,4413,4846,5058,5462,5629,6122,6422,6617,6933,7051,7156,7668];
offt2_10 = [91,362,437,647,796,932,1246,1263,1396,1726,1876,1894,1938,1983,2193,2536,2701,3121,3182,3452,3496,3766,3947,4038,4307,4471,4621,4953,4981,5328,5358,6031,6226,6391,6556,6842,6932,7007,7276,7352];
offt2_15 =[76,753,1127,1921,5703,4442,4459,4486,1098,4864,1310,4999,1384,5071,1581,514,1893,5433,1908,5464,1999,5569,2192,5851,2357,586,2462,6034,2522,6140,2746,6196,2809,6422,2868,6589,2881,6721,3061,6828,3139,7234,3226,7503,3258,7667,347,872,3601,94,3692,3813,4098];
offt2_31 =[5791,4142,1670,421,469,556,976,1006,1277,1576,1726,1922,2028,2356,2581,2656,2868,3046,3274,3482,3571,3871,4171,4366,4712,4999,5026,5281,5386,5748,6003,6271,6587,6631,6827,6964,7291,7502];
offt3_10 = [512,1368,2750,4188,4487,5732,6170,6511,6588];
offt3_15 =[62,991,1156,4234,4278,4487,5733,6286,6931];
offt3_31 =[1577,1954,2659,7324,6033,4383,2884];

offt4_10 = [2117,3076,4068,4983,5011,6241,7608,6617,2507,1279,3016];
offt4_15 = [889,1191, 1684,2478,2896,3200,3648,3679,3767,3935,4188,4788, 4294,5645,6064,6361,7278, 1084, 3635,7040,7069,5705];
offt4_31 =[681,1772,2178,2899,3946,6499,2973,7444,591,4036];

offt5_31 = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];
offt5_10 = [800,2597,3617,4321,4519,5536,5614,7353,3423,2029,6949];
offt5_15 = [[1246,2253,3695,5116,6260,4771,6380,226]];

off_otherotherotherotherother31 = [49,259,272,391,393,783,1430,1501,1741,1862,2177,2539,2597,2971,3661,3905,4066,4518,4578,4681,4727,4730,4771,4786,5072,5089,5495,5506,5596,5841,5853,5926,5943,6062,6091,6106,6212,6213,6541,6888,6901,7562,7593]
off_otherotherotherotherother10 = [17,454,496,1231,1459,1486,1624,1701,1801,1803,1895,2057,2073,2462,2794,3093,3363,3874,4367,4384,4580,4606,4996,5449,5746,5836,5959,6378,6483,6856,7336,7561];
off_otherotherotherotherother15 = [407,649,843,901,1130,1411,1427,1564,1817,2206,2236,2326,2343,2373,2521,2555,2732,2794,3091,3286,3857,3889,4097,4235,4697,4732,5148,5672,6139,6257,6376,6886,6992,7517];



on_10 = [31,77,211,242,244,272,558,631,646,677,813,1036,1081,1100,1144,1189,1281,1397,1442,1517,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2402,2435,2464,2491,2703,2716,2748,3002,3031,3213,3287,3316,3331,3422,3512,3515,3587,3618,3811,4024,4126,4127,4174,4262,4443,4488,4534,4548,4711,4712,4877,5087,5146,5176,5496,5506,5642,5671,5702,5777,5821,5868,5957,6136,6197,6272,6275,6572,6722,6751,6752,6767,6796,6977,7096,7189,7232,7426,7472,7564,7610,7621,7668];
on_15 = [4,31,154,182,272,454,496,528,692,766,783,857,860,888,979,1067,1111,1339,1385,1398,1486,1532,1549,1578,1591,1786,1876,1877,1892,1969,2011,2026,2136,2161,2208,2311,2328,2417,2461,2597,2716,2851,2929,3002,3152,3198,3241,3244,3303,3319,3422,3452,3530,3559,3586,3634,3691,3721,3812,3859,3917,3946,3991,3994,4022,4096,4130,4145,4231,4246,4427,4501,4562,4668,4774,4892,4941,4997,4998,5073,5088,5104,5179,5359,5405,5567,5657,5791,5853,5898,5927,5941,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6542,6722,6737,6752,6826,6903,6980,7067,7096,7157,7186,7203,7261,7354,7442,7487,7520,7532,7562];
on_31 = [4,78,108,286,362,543,618,647,781,785,903,917,962,991,1067,1203,1248,1306,1381,1471,1595,1606,1712,1756,1878,1952,2086,2087,2146,2371,2419,2494,2536,2582,2687,2719,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4113,4156,4172,4186,4248,4277,4384,4685,4697,4713,4714,4876,4892,4983,5028,5042,5132,5240,5267,5297,5328,5401,5552,5642,5746,5765,5792,5794,5866,5897,5956,5961,6002,6031,6211,6289,6497,6590,6646,6662,6692,6752,6753,6784,6812,6962,7067,7098,7114,7336,7442,7472,7640,7669];
on_other10 = [31,242,272,558,631,646,677,813,1036,1081,1100,1189,1281,1397,1442,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2435,2464,2491,2703,3002,3031,3213,3316,3422,3512,3515,3587,3618,3811,4024,4126,4174,4262,4443,4488,4534,4548,4711,4877,5146,5496,5506,5702,5777,5957,6136,6272,6275,6572,6722,6751,6752,6796,7096,7189,7232,7426,7472,7564,7610,7621];
on_other15 = [31,182,272,454,496,528,766,783,857,888,979,1067,1385,1398,1486,1532,1549,1578,1591,1876,2026,2136,2311,2328,2417,2597,2716,2929,3152,3198,3241,3244,3303,3452,3530,3559,3586,3634,3812,3859,3917,3946,3991,4022,4130,4145,4231,4246,4562,4668,4892,4941,4997,4998,5179,5405,5657,5791,5898,5927,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6722,6737,6752,6903,6980,7067,7096,7157,7186,7203,7354,7520,7562];
on_other31 = [108,362,543,647,781,785,903,962,991,1067,1203,1306,1381,1471,1595,1606,1756,1878,1952,2086,2146,2371,2419,2494,2536,2582,2687,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4113,4156,4186,4384,4685,4714,4876,4892,5028,5042,5132,5240,5267,5297,5328,5642,5746,5765,5792,5794,5866,5956,5961,6002,6031,6211,6497,6646,6662,6692,6752,6753,6784,6812,6962,7098,7114,7336,7442,7640,7669];
on_otherother10 = [242,558,631,677,813,1100,1189,1281,1594,1697,1759,1863,1939,2118,2464,3213,3422,3515,3618,4024,4126,4443,4488,4534,4548,4877,5496,5506,5777,6275,6752,7096,7189,7232,7564,7610,7621];
on_otherother15 = [31,182,496,528,766,857,1385,1532,1549,1578,1591,2136,2311,2328,2716,3152,3241,3452,3530,3586,3812,3859,4022,4130,4145,4231,4892,4941,4997,5179,5405,5791,6093,6125,6152,6155,6304,6363,6392,6439,6737,6752,6980,7096,7157,7186,7520];
on_otherother31 = [543,781,785,1067,1306,1381,1471,1595,1606,1878,1952,2146,2494,2536,2687,2747,2851,2896,3215,3244,3931,4113,4156,4186,4714,4892,5042,5132,5240,5267,5328,5765,5792,5794,5866,5961,6211,6646,6662,6692,6753,6784,6962,7114,7336,7640];
ont1_10 = [77,211,244,1144,1517,1531,2402,2716,2748,3287,3331,4127,4712,5087,5176,5642,5671,5821,5868,6197,6767,6977,7668];
ont1_15 = [4,154,692,860,1111,1339,1786,1877,1892,1969,2011,2161,2208,2461,2851,3002,3319,3422,3691,3721,3994,4096,4427,4501,4774,5073,5088,5104,5359,5567,5853,5941,6542,6826,7261,7442,7487,7532];
ont1_31 = [4,78,286,618,917,1248,1712,2087,2719,4172,4248,4277,4697,4713,4983,5401,5552,5897,6289,6590,7067,7472];

ont2_10 = [31,272,646,1036,1081,1397,1442,1696,2134,2177,2401,2435,2491,2703,3002,3031,3316,3512,3587,3811,4174,4262,4711,5146,5702,5957,6136,6272,6572,6722,6751,6796,7426,7472];
ont2_15 = [1398,2929,3559,3991,4998,6106,6722,7354,979,2026,3198,3634,4246,5657,6170,6903,7562,2417,3244,3917,454,5898,6455,7067,783,1067,272,3303,3946,4562,5927,6512,7203,888,1486,1876,2597,4668];
ont2_31 = [108,362,647,903,962,991,1203,1756,2086,2371,2419,2582,2856,2897,3076,3289,3843,3887,4384,4685,4876,5028,5297,5642,5746,5956,6002,6031,6497,6752,6812,7098,7442,7669];
% Initial conditions
ont3_31 = [1067 1306 1595 1878 2851 5328 6692 6753 ];
ont3_15 = [1385 1549 2136 3452 4892 4941 5179 5405 6093 6980 7157 7520];
ont3_10 = [1759 3213 4443 4534 6275 6752 7564];














