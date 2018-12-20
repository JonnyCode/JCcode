%% load dataruns - drifting grating, white noise, pulses

addpath('/Users/sneharavi/Documents/MATLAB/Classification/');
addpath('/Users/sneharavi/Documents/MATLAB/DS cell analysis/');

[datarun002_10] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-10-1/data001-3600-7200s/', 'data002-map/data002-map', 1,'stimuli/s02',0);
[datarun000_10] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-10-1/data001-3600-7200s/', 'data001-map/data001-map', 0,0,1);
[datarun001_10] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-10-1/data001-3600-7200s/', 'data003-map/data003-map', 0,0,0);

%% DS - 64 / 32 256

cellids_10 = intersect((intersect(datarun000_10.cell_ids, datarun001_10.cell_ids)), datarun002_10.cell_ids);
%cellids_10 = [17,31,77,91,211,242,244,272,317,362,395,437,454,481,496,512,558,571,631,646,647,677,695,706,783,796,800,813,917,932,995,1023,1036,1081,1100,1126,1142,1144,1189,1231,1246,1263,1265,1279,1281,1368,1369,1396,1397,1416,1442,1459,1486,1517,1518,1531,1594,1595,1624,1666,1681,1696,1697,1701,1712,1726,1759,1801,1803,1863,1876,1879,1894,1895,1938,1939,1983,2029,2057,2073,2117,2118,2134,2146,2177,2193,2401,2402,2435,2461,2462,2464,2491,2507,2536,2555,2597,2631,2701,2703,2716,2748,2750,2794,2836,2975,3002,3016,3017,3031,3076,3093,3121,3182,3213,3287,3306,3316,3331,3363,3421,3422,3423,3452,3470,3496,3512,3515,3574,3587,3617,3618,3736,3766,3769,3811,3828,3856,3874,3931,3947,4024,4038,4055,4068,4126,4127,4128,4174,4188,4261,4262,4307,4321,4323,4336,4366,4367,4384,4443,4471,4487,4488,4519,4534,4548,4549,4577,4580,4606,4621,4711,4712,4728,4771,4833,4877,4953,4981,4983,4996,5011,5087,5131,5146,5176,5328,5358,5375,5403,5419,5449,5465,5496,5506,5508,5536,5614,5627,5642,5671,5702,5732,5746,5761,5777,5821,5836,5868,5957,5959,6031,6047,6110,6136,6138,6170,6186,6196,6197,6226,6241,6272,6275,6319,6331,6378,6391,6438,6483,6511,6544,6556,6572,6588,6617,6633,6676,6722,6751,6752,6754,6767,6796,6813,6842,6856,6932,6949,6977,6992,7007,7052,7096,7114,7189,7232,7276,7336,7352,7353,7354,7381,7382,7426,7472,7561,7564,7608,7610,7621,7668];
[tc nontc] = get_time_courses_matrix(datarun000_10, cellids_10); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(cellids_10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

%DS cells
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_10, cellids_10, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
ds_init10 = [395,695,706,995,1023,1265,1369,1416,1595,1666,1681,1712,1879,2555,2631,2975,3306,3470,3574,3769,3828,3931,4323,4366,4549,4577,4833,5375,5419,5465,5627,6110,6138,6186,6544,6633,6754,7114,7381];

[C ia ib] = intersect(ds_init10, cellids_10);
vc = ones(length(cellids_10),1);
vc(ib) = 2;

close all;
X = [];
N = [];
p = [];
X(:,1) = log(mag{1,1})';
X(:,2) = log(mag{2,1})';
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_10, cellids_10, tcnormnorm,0, vc);

ds_10 = [];
ds_10 = cellids_10(idx==2);
nonds_10 = cellids_10(idx==1);
%lr = sum(p(idx==2,2))/length(ds_10)
%lr = sum(p(idx==1,1))/length(nonds_10)

%ds_10 = [395,695,706,995,1023,1265,1369,1416,1595,1666,1681,1712,1879,2555,2631,2975,3306,3470,3574,3769,3828,3931,4323,4366,4549,4577,4833,5375,5419,5465,5627,6110,6138,6186,6544,6633,6754,7114,7381]
%nonds_10 = [17,31,77,91,211,242,244,272,317,362,437,454,481,496,512,558,571,631,646,647,677,783,796,800,813,917,932,1036,1081,1100,1126,1142,1144,1189,1231,1246,1263,1279,1281,1368,1396,1397,1442,1459,1486,1517,1518,1531,1594,1624,1696,1697,1701,1726,1759,1801,1803,1863,1876,1894,1895,1938,1939,1983,2029,2057,2073,2117,2118,2134,2146,2177,2193,2401,2402,2435,2461,2462,2464,2491,2507,2536,2597,2701,2703,2716,2748,2750,2794,2836,3002,3016,3017,3031,3076,3093,3121,3182,3213,3287,3316,3331,3363,3421,3422,3423,3452,3496,3512,3515,3587,3617,3618,3736,3766,3811,3856,3874,3947,4024,4038,4055,4068,4126,4127,4128,4174,4188,4261,4262,4307,4321,4336,4367,4384,4443,4471,4487,4488,4519,4534,4548,4580,4606,4621,4711,4712,4728,4771,4877,4953,4981,4983,4996,5011,5087,5131,5146,5176,5328,5358,5403,5449,5496,5506,5508,5536,5614,5642,5671,5702,5732,5746,5761,5777,5821,5836,5868,5957,5959,6031,6047,6136,6170,6196,6197,6226,6241,6272,6275,6319,6331,6378,6391,6438,6483,6511,6556,6572,6588,6617,6676,6722,6751,6752,6767,6796,6813,6842,6856,6932,6949,6977,6992,7007,7052,7096,7189,7232,7276,7336,7352,7353,7354,7382,7426,7472,7561,7564,7608,7610,7621,7668]

%% ON - OFF Cells
temp_tcs = get_time_courses_matrix(datarun000_10, nonds_10);
tc_fit = [];
final_params  =[];
for i = 1:length(nonds_10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(nonds_10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(nonds_10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[minn minnt] = min(tcfittednormnorm);
[eval extt] = max(abs(tcfittednormnorm));
extrval = [];
for i = 1:length(extt)
extrval(i) = tcfittednormnorm(extt(i),i);
end


[tc nontc] = get_time_courses_matrix(datarun000_10, nonds_10); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(nonds_10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

on_init10 = [31,77,211,242,244,272,558,631,646,677,783,813,1036,1081,1100,1144,1189,1281,1397,1442,1517,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2402,2435,2464,2491,2703,2716,2748,3002,3031,3213,3287,3316,3331,3422,3512,3515,3587,3618,3811,4024,4055,4126,4127,4174,4262,4443,4488,4534,4548,4711,4712,4877,5087,5146,5176,5496,5506,5642,5671,5702,5777,5821,5868,5957,6136,6197,6272,6275,6319,6572,6722,6751,6752,6767,6796,6977,7096,7189,7232,7354,7426,7472,7564,7610,7621,7668];
[C ia ib] = intersect(on_init10, nonds_10);
vc = ones(length(nonds_10),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_10, nonds_10, tcnormnorm,0, vc);
on_allsnr10 = nonds_10(idx ==2)
off_allsnr10 = nonds_10(idx ==1)
%on_allsnr10 = [31,77,211,242,244,272,558,631,646,677,783,813,1036,1081,1100,1144,1189,1281,1397,1442,1517,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2402,2435,2464,2491,2703,2716,2748,3002,3031,3213,3287,3316,3331,3422,3512,3515,3587,3618,3811,4024,4055,4126,4127,4174,4262,4443,4488,4534,4548,4711,4712,4877,5087,5146,5176,5496,5506,5642,5671,5702,5777,5821,5868,5957,6136,6197,6272,6275,6319,6572,6722,6751,6752,6767,6796,6977,7096,7189,7232,7354,7426,7472,7564,7610,7621,7668];
%off_allsnr10 = [17,91,317,362,437,454,481,496,512,571,647,796,800,917,932,1126,1142,1231,1246,1263,1279,1368,1396,1459,1486,1518,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2146,2193,2461,2462,2507,2536,2597,2701,2750,2794,2836,3016,3017,3076,3093,3121,3182,3363,3421,3423,3452,3496,3617,3736,3766,3856,3874,3947,4038,4068,4128,4188,4261,4307,4321,4336,4367,4384,4471,4487,4519,4580,4606,4621,4728,4771,4953,4981,4983,4996,5011,5131,5328,5358,5403,5449,5508,5536,5614,5732,5746,5761,5836,5959,6031,6047,6170,6196,6226,6241,6331,6378,6391,6438,6483,6511,6556,6588,6617,6676,6813,6842,6856,6932,6949,6992,7007,7052,7276,7336,7352,7353,7382,7561,7608];
%% SNR CUTOFF
c = get_cell_indices(datarun000_10, on_allsnr10);
snronall = [];
for i = 1:length(c)
    r1 = sort(datarun000_10.stas.rfs{c(1,i),1}(:)', 'descend');
    snronall(1,i) = mean(r1(1:4))./std(r1);
end
on_10 = on_allsnr10(snronall > (mean(snronall) - 2.5*std(snronall)));

% hax=axes; 
% hold on;
% hist(snronall)
% SP= mean(snronall) - 2.5*std(snronall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('On cutoff')


c = get_cell_indices(datarun000_10, off_allsnr10);
snroffall = [];
for i = 1:length(c)
    r1 = sort(datarun000_10.stas.rfs{c(1,i),1}(:)', 'descend');
    snroffall(1,i) = mean(r1(1:4))./std(r1);
end
off_10 = off_allsnr10(snroffall > (mean(snroffall) - 2.5*std(snroffall)));

% hax=axes; 
% hold on;
% hist(snroffall)
% SP= mean(snroffall) - 2.5*std(snroffall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('Off cutoff')

snroff_10 = off_allsnr10(snroffall < (mean(snroffall) - 2.5*std(snroffall)));
snron_10 = on_allsnr10(snronall < (mean(snronall) - 2.5*std(snronall)));

% snron_10 = [783 4055 6319 7354];
% snroff_10 = [6331]
% on_10 = [31,77,211,242,244,272,558,631,646,677,813,1036,1081,1100,1144,1189,1281,1397,1442,1517,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2402,2435,2464,2491,2703,2716,2748,3002,3031,3213,3287,3316,3331,3422,3512,3515,3587,3618,3811,4024,4126,4127,4174,4262,4443,4488,4534,4548,4711,4712,4877,5087,5146,5176,5496,5506,5642,5671,5702,5777,5821,5868,5957,6136,6197,6272,6275,6572,6722,6751,6752,6767,6796,6977,7096,7189,7232,7426,7472,7564,7610,7621,7668]
% off_10 = [17,91,317,362,437,454,481,496,512,571,647,796,800,917,932,1126,1142,1231,1246,1263,1279,1368,1396,1459,1486,1518,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2146,2193,2461,2462,2507,2536,2597,2701,2750,2794,2836,3016,3017,3076,3093,3121,3182,3363,3421,3423,3452,3496,3617,3736,3766,3856,3874,3947,4038,4068,4128,4188,4261,4307,4321,4336,4367,4384,4471,4487,4519,4580,4606,4621,4728,4771,4953,4981,4983,4996,5011,5131,5328,5358,5403,5449,5508,5536,5614,5732,5746,5761,5836,5959,6031,6047,6170,6196,6226,6241,6378,6391,6438,6483,6511,6556,6588,6617,6676,6813,6842,6856,6932,6949,6992,7007,7052,7276,7336,7352,7353,7382,7561,7608];
%% ON - OFF after SNR check

nondssnr_10 = [on_10 off_10];

temp_tcs = get_time_courses_matrix(datarun000_10, nondssnr_10);
tc_fit = [];
final_params  =[];
for i = 1:length(nondssnr_10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(nondssnr_10) %fit time course
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(nondssnr_10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all fitted time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[minn minnt] = min(tcfittednormnorm);
[eval extt] = max(abs(tcfittednormnorm));
extrval = [];
for i = 1:length(extt)
extrval(i) = tcfittednormnorm(extt(i),i);
end


[tc nontc] = get_time_courses_matrix(datarun000_10, nondssnr_10); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(nondssnr_10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

on_initsnr10 = on_10;
[C ia ib] = intersect(on_initsnr10, nondssnr_10);
vc = ones(length(nondssnr_10),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_10, nondssnr_10, tcnormnorm,0, vc);
onon_10 = nondssnr_10(idx ==2); %idx might change so be careful - on might be idx 2 and off idx 1
offoff_snr10 = nondssnr_10(idx ==1);


    plot(X(idx==2,1),X(idx==2,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
    hold on;
        plot(X(idx==1,1),X(idx==1,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]); 
            legend('Cluster 1 - On Cells','Cluster 2 - Off Cells', 'Location','NW')

set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON-OFF/', 'classfaftersnr', gcf)



%% ON T1 ------

%2012-10-10-0: 19 t1, 4 not there
on_10 = [31,77,211,242,244,272,558,631,646,677,813,1036,1081,1100,1144,1189,1281,1397,1442,1517,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2402,2435,2464,2491,2703,2716,2748,3002,3031,3213,3287,3316,3331,3422,3512,3515,3587,3618,3811,4024,4126,4127,4174,4262,4443,4488,4534,4548,4711,4712,4877,5087,5146,5176,5496,5506,5642,5671,5702,5777,5821,5868,5957,6136,6197,6272,6275,6572,6722,6751,6752,6767,6796,6977,7096,7189,7232,7426,7472,7564,7610,7621,7668];
ont1_10_init = [77,244,1144,1517,2402,2716,3287,3331,4712,5087,5176,5642,5671,5821,6197,6977,7668,211, 1531,2748,5868, 6767,4127]; 
%  unsure t1s within hand IC: 1531 2748 5868; 211 tc is off   6767 4127 are ok but dg might be lower - doesn't change cluster with or without these cells

% pulses?
% Pm vs 1 - ratio 

temp_tcs = get_time_courses_matrix(datarun000_10, on_10);
tc_fit = [];
final_params  =[];
for i = 1:length(on_10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)


[tc nontc] = get_time_courses_matrix(datarun000_10, on_10); %or cellids
x = 1:1:30;
auc = [];
mx = [];
normval = [];
tcnormnorm = [];
tcnormauc = [];
tcnormmx = [];

for i = 1:length(on_10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
  auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 
[mx mxt] = max(tc);
auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
tcnormnorm = tc./normval;
tcnormauc = tc./auc;
tcnormmx = tc./mx;

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_10, on_10, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 32),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
close all;

vc = [];
[C ia ib] = intersect(ont1_10_init, on_10);
vc = ones(length(on_10),1);
vc(ib) = 2;

[COEFF,SCORE] = princomp(tcnormmx');


X = [];
X(:,1) = A32';
X(:,2) = A256';
%X(:,3) = SCORE(:,3);
%X(:,3) = (TCParams.maxval./TCParams.minval)';
%X(:,1) = TCParams.maxtim';
X(:,3) = TCParams.dot';
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, on_10, tcnormnorm,0, vc);
on_10(idx==2)
ont1_10 = on_10(idx==2);
on_other10 = on_10(idx==1);
ont1_10 =[77,211,244,1144,1517,2402,2716,3287,3331,4712,5087,5176,5642,5671,5821,6197,6767,6977,7668];
%on_other10 =[31,242,272,558,631,646,677,813,1036,1081,1100,1189,1281,1397,1442,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2435,2464,2491,2703,2748,3002,3031,3213,3316,3422,3512,3515,3587,3618,3811,4024,4126,4127,4174,4262,4443,4488,4534,4548,4711,4877,5146,5496,5506,5702,5777,5868,5957,6136,6272,6275,6572,6722,6751,6752,6796,7096,7189,7232,7426,7472,7564,7610,7621];
%%
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),49, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
hold on;
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),49, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('Average Spikes Fired - Fast Speed Drifting Grating')
ylabel('Average Spikes Fired - Slow Speed Drifting Grating')
zlabel('Degree of Transience of Time Course')

legend('Cluster 1 - ON T1 Cells','Cluster 2 - Other ON Cells', 'Location','NW')
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);

%%
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1classf2', gcf)

%%
plot_rf_summaries(datarun000_10, ont1_10, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1rf', gcf)

%%
plot_time_courses(datarun000_10,ont1_10, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1tc', gcf)

%%
  plot(x, tcnormnorm(:,idx==1), 'b')
  hold on;
 plot(x, tcnormnorm(:,idx==2), 'r')
 h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('TC of ON T1', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1tccomp', gcf)

%%
datarun000_10 = get_interspikeinterval(datarun000_10, ont1_10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont1_10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, ont1_10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
 maxisi(1,i) = max(isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
maxisi = repmat(maxisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;
isimax = isi./maxisi;

figure();
shadedErrorBar(x2(1:50),mean(isinormnorm(1:50, :)'),std(isinormnorm(1:50, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1isifull', gcf)

%%
wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); 
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,ont1_10), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:24; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont1_10)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end

b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

spikesbytrials{1,1} = get_raster(datarun001_10.spikes{get_cell_indices(datarun001_10, ont1_10(1,3)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 24 0 4]);
hold on;
stairs([0 3 12 15 24],[w g b g g], 'Color', 'r', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1psthcell244', gcf);
close;


figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
plot(binSize,psthnorm(3,:), 'b');
stairs([0 3 12 15 24],[0.4 0.35 0.3 0.35 0.35], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1psth2', gcf)

%% ON T2
%dg, rf, time course params dot area zc time value, tc pc
%later: isi, pulses

on_other10 =[31,242,272,558,631,646,677,813,1036,1081,1100,1189,1281,1397,1442,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2435,2464,2491,2703,2748,3002,3031,3213,3316,3422,3512,3515,3587,3618,3811,4024,4126,4127,4174,4262,4443,4488,4534,4548,4711,4877,5146,5496,5506,5702,5777,5868,5957,6136,6272,6275,6572,6722,6751,6752,6796,7096,7189,7232,7426,7472,7564,7610,7621];
ont2_10_init = [31,272,646,1036,1081,1397,1442,1696,2134,2177,2401,2435,2491,2703,3002,3031,3316,3512,3587,3811,4174,4262,4711,5146,5702,5957,6136,6272,6572,6722,6751,6796,7426,7472];

 temp_tcs = get_time_courses_matrix(datarun000_10, on_other10);
tc_fit = [];
final_params  =[];
for i = 1:length(on_other10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_other10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_other10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)
 

 [tc nontc] = get_time_courses_matrix(datarun000_10,  on_other10); %or cellids
x = 1:1:30;
auc = [];
mx = [];
normval = [];
tcnormnorm = [];
tcnormauc = [];
tcnormmx = [];
for i = 1:length(on_other10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 
[mx mxt] = max(tc);
auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormmx = tc./mx;
% 
% [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_10, on_other10, 0);
% [mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
% NS2 = [];
% A32 = [];
% A256 = [];
% NS2 = NumSpikesCell';
% A32 = sum(NS2(find(StimComb(:,2) == 32),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
% A256 = sum(NS2(find(StimComb(:,2) == 256),:));
% close all;


datarun000_10 = get_interspikeinterval(datarun000_10, on_other10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
for i = 1:length(on_other10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, on_other10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

% close all;
% rstd = [];
% meanpix = [];
% cellind = get_cell_indices(datarun000_10, on_other10);
% stamat = cell(length(cellind),1);
%  
% for i = 1:length(cellind)
%     B = [];
%     B = datarun000_10.stas.rfs{cellind(i), 1}(:)';
%     meanpix(i) = mean(B);
%     rstd(i) = robust_std(B, [1]);
%     stamat{i,1} = zeros(size(datarun000_10.stas.rfs{cellind(i),1},1),size(datarun000_10.stas.rfs{cellind(i),1},2));
%     for j = 1:size(datarun000_10.stas.rfs{cellind(i),1},1)
%             for k = 1:size(datarun000_10.stas.rfs{cellind(i),1},2)
%                 if (datarun000_10.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
%                     stamat{i,1}(j,k) = 1;
%                 end
%             end
%     end
%  
%     
% end
%     
%  
%  
%  
% close all;
% v = [];
% pm = [];
% for i = 1:length(cellind)
%     rm = [];
%     rn = [];
%     rmm = [];
%     rnn = [];
%     DT = [];
%     kr = [];
%     points = [];
%     [rm,rn] = find(stamat{i,1}); %finds all the significant pixels
%             if ~isempty(rm)
%                 [rmm,rnn] = pix_border(rm,rn); %adding everything by .5 and subtracting by 0.5 - i think it is calculating the 4 coordinate values for each pixel - edge correction
%                 DT = delaunayTriangulation(rmm,rnn); %prediction is still correct with inputs rm and rn, but edges need correcting
%                 if ~isempty(DT.ConnectivityList)
%                     [kr v(1,i)] = convexHull(DT);
%                 end
%                 %spy(stamat{i,1}, 'LineSpec', 'r');
%                 [x,y] = find(stamat{i,1});
%                 clr = stamat{i,1}(stamat{i,1}~=0);
%                 %scatter(y,x,20,'MarkerFaceColor',[rand(1) rand(1) rand(1)],'LineWidth',0.05)
%                 %set(gca,'Xdir','reverse');%'Ydir','reverse')
%                 %plot(DT.Points(:,2),DT.Points(:,1), '.','markersize',3);
%                 hold on;
%                 %plot(DT.Points(kr,2), DT.Points(kr,1), 'Color',[rand(1) rand(1) rand(1)]);
%                 %hold off;
%                 %pause;
%                 points(:,1) = DT.Points(kr,2);
%                 points(:,2) = DT.Points(kr,1);
%                 perimeter = 0;
%                 for j = 1:size(points, 1)-1
%                     perimeter = perimeter + norm(points(j, :) - points(j+1, :));
%                 end
%                 perimeter = perimeter + norm(points(end, :) - points(1, :)); % Last point to first
%                 pm(i) = perimeter;
%             else
%             end
%             v(2,i) = length(rm);
% end
% 
% radius = [];
% radius = get_rf_fit_radius(datarun000_10, on_other10);

vc = [];
[C ia ib] = intersect(ont2_10_init, on_other10);
vc = ones(length(on_other10),1);
vc(ib) = 2;

 
[COEFF,SCORE] = princomp(tcnormnorm');
[COEFF1,SCORE1] = princomp(isinormnorm');
%  minval maxtim maxval zc tc auc 1-2-3  pm 
 X = [];
X(:,1) = TCParams.minval;
X(:,2) = TCParams.zerocrossing;
X(:,3) =  SCORE1(:,1);
% % X(:,1) = SCORE(:,3);
% % X(:,2) = SCORE(:,3);
% X(:,3) = TCParams.zerocrossing;
%X(:,3) = rf;
%X(:,4) = pm;
%X(:,5) = 1-(v(2,:)./v(1,:));
 [idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, on_other10, tcnormnorm,0, vc);
 on_other10(idx==2)
ont2_10 = on_other10(idx==2);
on_otherother10 = on_other10(idx==1);
%zc score 3; zc pm ;
%minval pm; minval score 3
%dot SRF has good separation - look into 3d - with score 1 or 3
%score 1 - pm; 1 - ratio
ont2_10 =[31,272,646,1036,1081,1397,1442,1696,2177,2401,2435,2491,2703,3002,3031,3316,3512,3587,3811,4174,4262,4711,5146,5702,5957,6136,6272,6572,6722,6751,6796,7426,7472];
%on_otherother10 =[242,558,631,677,813,1100,1189,1281,1531,1594,1697,1759,1863,1939,2118,2134,2464,2748,3213,3422,3515,3618,4024,4126,4127,4443,4488,4534,4548,4877,5496,5506,5777,5868,6275,6752,7096,7189,7232,7564,7610,7621];

 %% on t2 extra cells
  plot(x, tcnormnorm(:,idx==1), 'b')
  hold on;
 plot(x, tcnormnorm(:,idx==2), 'r')
 h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('ISI of ON T2', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2tccomp', gcf)

%%
on_extra10 = [ 1291 241 1802 4681 5222 5237 5867];
datarun000_10 = get_interspikeinterval(datarun000_10, on_extra10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
for i = 1:length(on_extra10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, on_extra10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

 plot(x2, isinormnorm(:,:), 'k');
 h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('ISI of ON T2 - extra', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2isiextra', gcf)

%%
 plot_rf_summaries(datarun000_10, [ 1291 241 1802 4681 5222 5237 5867], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2extrarf', gcf)
 
 
 plot_time_courses(datarun000_10, [ 1291 241 1802 4681 5222 5237 5867],'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2extratc', gcf)

plot_time_courses(datarun000_10,  [ont1_10 242 2118 3752 4216 4921 7381 7501 7576],'all', true, 'bw', true);
ismember([4216 7381 7501 7576], cellids_10)
ismember([5191 5868 6571], cellids_10)

%%
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),49, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
hold on;
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),49, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('Minimum Value - Time Course')
ylabel('Zero Crossing - Time Course ')
zlabel('PC-1 - Interspike interval')

legend('Cluster 1 - ON T2 Cells','Cluster 2 - Other ON Cells', 'Location','NW')
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);

%%
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2classf2', gcf)

%%
plot_rf_summaries(datarun000_10, ont2_10, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2rf', gcf)

%%
plot_time_courses(datarun000_10,ont2_10, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2tc', gcf)

%%
datarun000_10 = get_interspikeinterval(datarun000_10, ont2_10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont2_10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, ont2_10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
 maxisi(1,i) = max(isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
maxisi = repmat(maxisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;
isimax = isi./maxisi;

figure();
shadedErrorBar(x2(1:50),mean(isinormnorm(1:50, :)'),std(isinormnorm(1:50, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2isifull', gcf)


%%
wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); 
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,ont2_10), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:24; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont2_10)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end

b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
plot(binSize,psthnorm(5,:), 'b');
stairs([0 3 12 15 24],[0.55 0.5 0.45 0.5 0.5], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2psth2', gcf)
close;


spikesbytrials{1,1} = get_raster(datarun001_10.spikes{get_cell_indices(datarun001_10, ont2_10(1,5)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 24 0 4]);
hold on;
stairs([0 3 12 15 24],[w g b g g], 'Color', 'r', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T2/', 'ont2psthcell1081', gcf);
close;




%%
on_otherother10 =[242,558,631,677,813,1100,1189,1281,1531,1594,1697,1759,1863,1939,2118,2134,2464,2748,3213,3422,3515,3618,4024,4126,4127,4443,4488,4534,4548,4877,5496,5506,5777,5868,6275,6752,7096,7189,7232,7564,7610,7621];
ont3_10_init = [1759 3213 4443 4534 6275 6752 7564];


datarun000_10 = get_interspikeinterval(datarun000_10, on_otherother10);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(on_otherother10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, on_otherother10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


  temp_tcs = get_time_courses_matrix(datarun000_10, on_otherother10);
tc_fit = [];
final_params  =[];
for i = 1:length(on_otherother10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_otherother10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_otherother10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)


[tc nontc] = get_time_courses_matrix(datarun000_10, on_otherother10); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(on_otherother10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

% [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_10, on_otherother10, 0);
% [mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
% NS2 = [];
% A32 = [];
% A256 = [];
% NS2 = NumSpikesCell';
% A32 = sum(NS2(find(StimComb(:,2) == 32),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
% A256 = sum(NS2(find(StimComb(:,2) == 256),:));
% close all;

% close all;
% rstd = [];
% meanpix = [];
% cellind = get_cell_indices(datarun000_10, on_otherother10);
% stamat = cell(length(cellind),1);
%  
% for i = 1:length(cellind)
%     B = [];
%     B = datarun000_10.stas.rfs{cellind(i), 1}(:)';
%     meanpix(i) = mean(B);
%     rstd(i) = robust_std(B, [1]);
%     stamat{i,1} = zeros(size(datarun000_10.stas.rfs{cellind(i),1},1),size(datarun000_10.stas.rfs{cellind(i),1},2));
%     for j = 1:size(datarun000_10.stas.rfs{cellind(i),1},1)
%             for k = 1:size(datarun000_10.stas.rfs{cellind(i),1},2)
%                 if (datarun000_10.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
%                     stamat{i,1}(j,k) = 1;
%                 end
%             end
%     end
%  
%     
% end
%     
%  
%  
%  
% close all;
% v = [];
% pm = [];
% for i = 1:length(cellind)
%     rm = [];
%     rn = [];
%     rmm = [];
%     rnn = [];
%     DT = [];
%     kr = [];
%     points = [];
%     [rm,rn] = find(stamat{i,1}); %finds all the significant pixels
%             if ~isempty(rm)
%                 [rmm,rnn] = pix_border(rm,rn); %adding everything by .5 and subtracting by 0.5 - i think it is calculating the 4 coordinate values for each pixel - edge correction
%                 DT = delaunayTriangulation(rmm,rnn); %prediction is still correct with inputs rm and rn, but edges need correcting
%                 if ~isempty(DT.ConnectivityList)
%                     [kr v(1,i)] = convexHull(DT);
%                 end
%                 %spy(stamat{i,1}, 'LineSpec', 'r');
%                 [x,y] = find(stamat{i,1});
%                 clr = stamat{i,1}(stamat{i,1}~=0);
%                 %scatter(y,x,20,'MarkerFaceColor',[rand(1) rand(1) rand(1)],'LineWidth',0.05)
%                 %set(gca,'Xdir','reverse');%'Ydir','reverse')
%                 %plot(DT.Points(:,2),DT.Points(:,1), '.','markersize',3);
%                 hold on;
%                 %plot(DT.Points(kr,2), DT.Points(kr,1), 'Color',[rand(1) rand(1) rand(1)]);
%                 %hold off;
%                 %pause;
%                 points(:,1) = DT.Points(kr,2);
%                 points(:,2) = DT.Points(kr,1);
%                 perimeter = 0;
%                 for j = 1:size(points, 1)-1
%                     perimeter = perimeter + norm(points(j, :) - points(j+1, :));
%                 end
%                 perimeter = perimeter + norm(points(end, :) - points(1, :)); % Last point to first
%                 pm(i) = perimeter;
%             else
%             end
%             v(2,i) = length(rm);
% end
% 
% radius = [];
% radius = get_rf_fit_radius(datarun000_10, on_otherother10);

[C ia ib] = intersect(ont3_10_init, on_otherother10);
vc = ones(length(on_otherother10),1);
vc(ib) = 2;

%pm maxmingrad dot not good
%maxval minval tc pc space zerocrossing 1-ratio

[COEFF1,SCORE1] = princomp(isinormnorm');
% [COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF,SCORE] = princomp(tcnormminn');

X = [];
X(:,1) = TCParams.minval;
X(:,2) = TCParams.zerocrossing;
X(:,3) = SCORE1(:,1);
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, on_otherother10, tcnormnorm,0,vc);
 %ismember(ont3_10_init, on_otherother10)
on_otherother10(idx==2)
ont3_10 = on_otherother10(idx==2);
on_otherotherother10 = on_otherother10(idx==1);
%  ismember(ont3_10_init, ont3_10);
ont3_10 = [1759,2134,3213,4443,4534,6275,6752,7564];
on_otherotherother10=[242,558,631,677,813,1100,1189,1281,1531,1594,1697,1863,1939,2118,2464,2748,3422,3515,3618,4024,4126,4127,4488,4548,4877,5496,5506,5777,5868,7096,7189,7232,7610,7621];
 %% on t3 extra cells
 plot_rf_summaries(datarun000_10, [6425,7384, 468, 842, 1129,1278,1384,1414,1925, 1967,2374,2687,2943,3125,3721,3887,4819,5148,5493,6137], 'coordinates', 'monitor', 'label', true);
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3extrarf', gcf)
 
 
 plot_time_courses(datarun000_10, [6425,7384, 468, 842, 1129,1278,1384,1414,1925, 1967,2374,2687,2943,3125,3721,3887,4819,5148,5493,6137],'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3extratc', gcf)

 %%
 scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),49, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
hold on;
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),49, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('Minimum Value - Time Course')
ylabel('Zero Crossing - Time Course ')
zlabel('PC-1 - Interspike interval')

legend('Cluster 1 - ON T3 Cells','Cluster 2 - Other ON Cells', 'Location','NW')
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);

%%
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3classf', gcf)

%%
plot_rf_summaries(datarun000_10, ont3_10, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3rf', gcf)

%%
close all;
plot_time_courses(datarun000_10,ont3_10, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3tc', gcf)

%%
datarun000_10 = get_interspikeinterval(datarun000_10, ont3_10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont3_10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, ont3_10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
 maxisi(1,i) = max(isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
maxisi = repmat(maxisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;
isimax = isi./maxisi;

figure();
shadedErrorBar(x2(1:50),mean(isinormnorm(1:50, :)'),std(isinormnorm(1:50, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3isifull', gcf)
%%
wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); 
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,ont3_10), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:24; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont3_10)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end

b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
plot(binSize,psthnorm(5,:), 'b');
stairs([0 3 12 15 24],[0.4 0.35 0.3 0.35 0.35], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3psth2', gcf)

spikesbytrials{1,1} = get_raster(datarun001_10.spikes{get_cell_indices(datarun001_10, ont3_10(1,5)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 24 0 4]);
hold on;
stairs([0 3 12 15 24],[w g b g g], 'Color', 'r', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T3/', 'ont3psthcell4534', gcf);
close;



%% OFF T1
off_10 = [17,91,317,362,437,454,481,496,512,571,647,796,800,917,932,1126,1142,1231,1246,1263,1279,1368,1396,1459,1486,1518,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2146,2193,2461,2462,2507,2536,2597,2701,2750,2794,2836,3016,3017,3076,3093,3121,3182,3363,3421,3423,3452,3496,3617,3736,3766,3856,3874,3947,4038,4068,4128,4188,4261,4307,4321,4336,4367,4384,4471,4487,4519,4580,4606,4621,4728,4771,4953,4981,4983,4996,5011,5131,5328,5358,5403,5449,5508,5536,5614,5732,5746,5761,5836,5959,6031,6047,6170,6196,6226,6241,6378,6391,6438,6483,6511,6556,6588,6617,6676,6813,6842,6856,6932,6949,6992,7007,7052,7276,7336,7352,7353,7382,7561,7608];
offt1_10_init = [317,481,571,917,1126,1142,1518,2146,2461,2836,3017,3421,3736,3856,4128,4261,4336,4728,4771,5131,5403,5508,5761,6047,6196,6438,6676,6813,6992,7052,7382];
 
 [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_10, off_10, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 32),:)); %CHANGE ACCORDING TO CORRECT TEMPORAL PERIOD
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
close all;


 temp_tcs = get_time_courses_matrix(datarun000_10, off_10);
tc_fit = [];
final_params  =[];
for i = 1:length(off_10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)



[tc nontc] = get_time_courses_matrix(datarun000_10, off_10); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
tcnormauc = [];
for i = 1:length(off_10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
  auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 
normval = repmat(normval, size(tc, 1), 1);
auc = repmat(auc, size(tc, 1), 1);

tcnormnorm = tc./normval;
tcnormauc = tc./auc;



vc = [];
[C ia ib] = intersect(offt1_10_init, off_10);
vc = ones(length(off_10),1);
vc(ib) = 2;

%[COEFF,SCORE] = princomp(tcnormauc');


X = [];
X(:,1) = A32';
X(:,2) = A256';
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_10, tcnormnorm,0, vc);
off_10(idx==2)
offt1_10 = off_10(idx==2);
off_other10 = off_10(idx==1);
offt1_10 =[317,481,571,917,1126,1142,1518,2146,2461,2836,3017,3421,3736,3856,4128,4261,4336,4606,4728,4771,5131,5403,5508,5761,6047,6196,6438,6676,6813,6992,7052,7382];
%off_other10 =[17,91,362,437,454,496,512,647,796,800,932,1231,1246,1263,1279,1368,1396,1459,1486,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2193,2462,2507,2536,2597,2701,2750,2794,3016,3076,3093,3121,3182,3363,3423,3452,3496,3617,3766,3874,3947,4038,4068,4188,4307,4321,4367,4384,4471,4487,4519,4580,4621,4953,4981,4983,4996,5011,5328,5358,5449,5536,5614,5732,5746,5836,5959,6031,6170,6226,6241,6378,6391,6483,6511,6556,6588,6617,6842,6856,6932,6949,7007,7276,7336,7352,7353,7561,7608];
 %%
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),49, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
hold on;
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),49, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('Average Spikes Fired - Fast Speed Drifting Grating')
ylabel('Average Spikes Fired - Slow Speed Drifting Grating')
zlabel('Degree of Transience of Time Course')

legend('Cluster 1 - OFF T1 Cells','Cluster 2 - Other OFF Cells', 'Location','NW')
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);

%%
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1classf2', gcf)

%%
plot_rf_summaries(datarun000_10, offt1_10, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1rf', gcf);

close all;
plot_time_courses(datarun000_10,offt1_10, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1tc', gcf)

%% extra cells off t1
plot_rf_summaries(datarun000_10, [offt1_10 2162 ], 'coordinates', 'monitor');
plot_time_courses(datarun000_10,[offt1_10 2162 ], 'all', true, 'bw', true);

%%
plot_rf_summaries(datarun000_10, [2162], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1extrarf', gcf);

close all;
plot_time_courses(datarun000_10,[2162], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1extratc', gcf)


%%
datarun000_10 = get_interspikeinterval(datarun000_10, offt1_10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt1_10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, offt1_10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
 maxisi(1,i) = max(isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
maxisi = repmat(maxisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;
isimax = isi./maxisi;

figure();
shadedErrorBar(x2(1:50),mean(isinormnorm(1:50, :)'),std(isinormnorm(1:50, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1isifull', gcf)
%%
wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); 
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,offt1_10), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:24; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt1_10)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end

b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
plot(binSize,psthnorm(1,:), 'b');
stairs([0 3 12 15 24],[0.4 0.35 0.3 0.35 0.35], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1psth2', gcf)

spikesbytrials{1,1} = get_raster(datarun001_10.spikes{get_cell_indices(datarun001_10, offt1_10(1,1)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 24 0 4]);
hold on;
stairs([0 3 12 15 24],[w g b g g], 'Color', 'r', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T1/', 'offt1psthcell317', gcf);
close;



%% OFF T2
offt2_10_init = [91,362,437,647,796,932,1246,1263,1396,1726,1876,1894,1938,1983,2193,2536,2701,3121,3182,3452,3496,3766,3947,4038,4307,4471,4621,4953,4981,5328,5358,6031,6226,6391,6556,6842,6932,7007,7276,7352];
off_other10 =[17,91,362,437,454,496,512,647,796,800,932,1231,1246,1263,1279,1368,1396,1459,1486,1624,1701,1726,1801,1803,1876,1894,1895,1938,1983,2029,2057,2073,2117,2193,2462,2507,2536,2597,2701,2750,2794,3016,3076,3093,3121,3182,3363,3423,3452,3496,3617,3766,3874,3947,4038,4068,4188,4307,4321,4367,4384,4471,4487,4519,4580,4621,4953,4981,4983,4996,5011,5328,5358,5449,5536,5614,5732,5746,5836,5959,6031,6170,6226,6241,6378,6391,6483,6511,6556,6588,6617,6842,6856,6932,6949,7007,7276,7336,7352,7353,7561,7608];

[tc nontc] = get_time_courses_matrix(datarun000_10, off_other10); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
tcnormauc = [];
for i = 1:length(off_other10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
  auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 
normval = repmat(normval, size(tc, 1), 1);
auc = repmat(auc, size(tc, 1), 1);

tcnormnorm = tc./normval;
tcnormauc = tc./auc;


temp_tcs = get_time_courses_matrix(datarun000_10, off_other10);
tc_fit = [];
final_params  =[];
for i = 1:length(off_other10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_other10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_other10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0);
datarun000_10 = get_interspikeinterval(datarun000_10, off_other10);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_other10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, off_other10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

radius = [];
radius = get_rf_fit_radius(datarun000_10, off_other10);

vc = [];
[C ia ib] = intersect(offt2_10_init, off_other10);
vc = ones(length(off_other10),1);
vc(ib) = 2;

%[COEFF,SCORE] = princomp(tcnormnorm');

%minval maxval zc dot maxt mint
% 
X = [];
X(:,1) = TCParams.dot;
X(:,2) = TCParams.maxtim;
X(:,3) = TCParams.mintim;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_other10, tcnormnorm,0, vc);
%off_other10(idx==2)
offt2_10 = off_other10(idx==2);
 off_otherother10 = off_other10(idx==1);%%
offt2_10 = [91,362,437,647,796,932,1246,1263,1396,1726,1876,1894,1938,1983,2193,2536,2701,3121,3182,3452,3496,3766,3947,4038,4307,4471,4621,4953,4981,5328,5358,6031,6226,6391,6556,6842,6932,7007,7276,7352];
% off_otherother10 = [17,454,496,512,800,1231,1279,1368,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2117,2462,2507,2597,2750,2794,3016,3076,3093,3363,3423,3617,3874,4068,4188,4321,4367,4384,4487,4519,4580,4983,4996,5011,5449,5536,5614,5732,5746,5836,5959,6170,6241,6378,6483,6511,6588,6617,6856,6949,7336,7353,7561,7608];

%%
plot_rf_summaries(datarun000_10, [offt2_10], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T2/', 'offt2rf', gcf);

close all;
plot_time_courses(datarun000_10,[offt2_10], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T2/', 'offt2tc', gcf)
%% extra cells off t2
plot_rf_summaries(datarun000_10, [offt2_10 4728], 'coordinates', 'monitor');
plot_time_courses(datarun000_10,[offt2_10 4728], 'all', true, 'bw', true);

%%
plot_rf_summaries(datarun000_10, [4728], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T2/', 'offt2extrarf', gcf);

close all;
plot_time_courses(datarun000_10,[4728], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T2/', 'offt2extratc', gcf)

 
%%
datarun000_10 = get_interspikeinterval(datarun000_10, offt2_10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt2_10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, offt2_10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
 maxisi(1,i) = max(isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
maxisi = repmat(maxisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;
isimax = isi./maxisi;

figure();
shadedErrorBar(x2(1:50),mean(isinormnorm(1:50, :)'),std(isinormnorm(1:50, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T2/', 'offt2isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T2/', 'offt2isifull', gcf)
%%
wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); 
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,offt2_10), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:24; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt2_10)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end

b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
%plot(binSize,psthnorm(1,:), 'b');
stairs([0 3 12 15 24],[0.4 0.35 0.3 0.35 0.35], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T2/', 'offt2psth2', gcf)

spikesbytrials{1,1} = get_raster(datarun001_10.spikes{get_cell_indices(datarun001_10, offt2_10(1,1)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 24 0 4]);
hold on;
stairs([0 3 12 15 24],[w g b g g], 'Color', 'r', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T2/', 'offt2psthcell91', gcf);
close;



 %%
off_otherother10 = [17,454,496,512,800,1231,1279,1368,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2117,2462,2507,2597,2750,2794,3016,3076,3093,3363,3423,3617,3874,4068,4188,4321,4367,4384,4487,4519,4580,4983,4996,5011,5449,5536,5614,5732,5746,5836,5959,6170,6241,6378,6483,6511,6588,6617,6856,6949,7336,7353,7561,7608];
offt4_10_init = [2117,3076,4068,4983,5011,6241,7608,6617,2507,1279,3016];

datarun000_10 = get_interspikeinterval(datarun000_10, off_otherother10);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherother10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, off_otherother10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_10, off_otherother10); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(off_otherother10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

temp_tcs = get_time_courses_matrix(datarun000_10, off_otherother10);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherother10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherother10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherother10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);

% wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
% gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); %Might change with dataset
% [h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,off_otherother10), 0, '/0', wh, gr, 10, false);
% binSize = 0.1:0.1:10; %change depending on length of trial
% pulsePSTH = [];
% normvalpulsePSTH = [];
% pulsenormnormPSTH = [];
% maxpulse = [];
% maxpulsetime =[];
%  for a = 1:length(off_otherother10)
%  pulsePSTH(:,a) = sum(nhist{a,1})./50; %change depending on num of trials
% end
%  
% for i = 1:length(off_otherother10)
%  normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
% end
% normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
% pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;
% 
% [maxpulse maxpulsetime] = max(pulsePSTH);
% maxpulsetime = maxpulsetime*0.1;

 
[C ia ib] = intersect(offt4_10_init, off_otherother10);
vc = ones(length(off_otherother10),1);
vc(ib) = 2;


[COEFF1,SCORE1] = princomp(isinormnorm');
%[COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF2,SCORE2] = princomp(tcnormminn');

X = [];
X(:,1) = SCORE1(:,1);
X(:,2) = TCParams.minval;
X(:,3) = TCParams.mintim;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherother10, tcnormnorm,0,vc);
 offt4_10 = off_otherother10(idx==2);
 off_otherotherother10 = off_otherother10(idx==1);
 %off_otherother10(idx==2)
%ismember(offt4_10_init, off_otherother10(idx==2))
 %off_otherotherother10 = [17,454,496,512,800,1231,1368,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2462,2597,2750,2794,3016,3093,3363,3423,3617,3874,4188,4321,4367,4384,4487,4519,4580,4996,5449,5536,5614,5732,5746,5836,5959,6170,6378,6483,6511,6588,6856,6949,7336,7353,7561];
 offt4_10 = [1279,2117,2507,3076,4068,4983,5011,6241,6617,7608];

 %% extra cells off t4
plot_rf_summaries(datarun000_10, [offt4_10  3 3901 4115 5117 7667], 'coordinates', 'monitor');
plot_time_courses(datarun000_10,[offt4_10 3 3901 4115 5117 7667 ], 'all', true, 'bw', true);
%%
plot_rf_summaries(datarun000_10, [3 3901 4115 5117 7667], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T4/', 'offt4extrarf', gcf);

close all;
plot_time_courses(datarun000_10,[3 3901 4115 5117 7667], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T4/', 'offt4extratc', gcf)

%%
plot_rf_summaries(datarun000_10, [offt4_10], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T4/', 'offt4rf', gcf);

close all;
plot_time_courses(datarun000_10,[offt4_10], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T4/', 'offt4tc', gcf)

%%
datarun000_10 = get_interspikeinterval(datarun000_10, offt4_10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt4_10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, offt4_10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
 maxisi(1,i) = max(isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
maxisi = repmat(maxisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;
isimax = isi./maxisi;

figure();
shadedErrorBar(x2(1:50),mean(isinormnorm(1:50, :)'),std(isinormnorm(1:50, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T4/', 'offt4isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T4/', 'offt4isifull', gcf)
%%
wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); 
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,offt4_10), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:24; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt4_10)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end

b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
%plot(binSize,psthnorm(1,:), 'b');
stairs([0 3 12 15 24],[0.4 0.35 0.3 0.35 0.35], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T4/', 'offt4psth2', gcf)

spikesbytrials{1,1} = get_raster(datarun001_10.spikes{get_cell_indices(datarun001_10, offt4_10(1,5)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 24 0 4]);
hold on;
stairs([0 3 12 15 24],[w g b g g], 'Color', 'r', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T4/', 'offt4psthcell4068', gcf);
close;


%%
off_otherotherother10 = [17,454,496,512,800,1231,1368,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2462,2597,2750,2794,3016,3093,3363,3423,3617,3874,4188,4321,4367,4384,4487,4519,4580,4996,5449,5536,5614,5732,5746,5836,5959,6170,6378,6483,6511,6588,6856,6949,7336,7353,7561];
offt3_10_init =[512,1368,2750,4188,4487,5732,6170,6511,6588];

[tc nontc] = get_time_courses_matrix(datarun000_10, off_otherotherother10); %or cellids
x = 1:1:30;
 normval = [];
 mx = [];
tcnormnorm = [];
tcnormmx = [];
for i = 1:length(off_otherotherother10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);

tcnormnorm = tc./normval;
[mx mxt] = max(tc);
mx = repmat(mx, size(tc, 1), 1);
tcnormmx = tc./mx;

temp_tcs = get_time_courses_matrix(datarun000_10, off_otherotherother10);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherother10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherother10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherother10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0);


vc = [];
[C ia ib] = intersect(offt3_10_init, off_otherotherother10);
vc = ones(length(off_otherotherother10),1);
vc(ib) = 2;
%minval maxval zc dot maxt mint

%[COEFF,SCORE] = princomp(tcnormmx');

X = [];
X(:,1) = TCParams.maxmingrad;
X(:,2) = TCParams.maxtim;
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherotherother10, tcnormnorm,0, vc);
%off_otherotherother10(idx==2)
offt3_10 = off_otherotherother10(idx==2);
off_otherotherotherother10 = off_otherotherother10(idx==1);
offt3_10 = [512,1368,2750,4188,4487,5732,6170,6511,6588];
off_otherotherotherother10 = [17,454,496,800,1231,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2462,2597,2794,3016,3093,3363,3423,3617,3874,4321,4367,4384,4519,4580,4996,5449,5536,5614,5746,5836,5959,6378,6483,6856,6949,7336,7353,7561];
 %% extra cells off t3
plot_rf_summaries(datarun000_10, [offt3_10 1579 2012 3199 3710 4160 6079 6151 7324 1580], 'coordinates', 'monitor');
plot_time_courses(datarun000_10,[offt3_10 1579 2012 3199 3710 4160 6079 6151 7324 1580], 'all', true, 'bw', true);


%%
plot_rf_summaries(datarun000_10, [1579 2012 3199 3710 4160 6079 6151 7324 1580], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/', 'offt3extrarf', gcf)

%%
close all;
plot_time_courses(datarun000_10,[1579 2012 3199 3710 4160 6079 6151 7324 1580], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/', 'offt3tcextra', gcf)





%% OFF T5
%off_otherotherotherother10
%offt5_10_init = [800,2029,2597,3423,3617,4321,4519,5536,5614,6949,7353];
%plot_time_courses(datarun000_10, [800,2029,2597,3423,3617,4321,4519,5536,5614,6949,7353,3874], 'all', true, 'bw', true);
off_otherotherotherother10 = [17,454,496,800,1231,1459,1486,1624,1701,1801,1803,1895,2029,2057,2073,2462,2597,2794,3016,3093,3363,3423,3617,3874,4321,4367,4384,4519,4580,4996,5449,5536,5614,5746,5836,5959,6378,6483,6856,6949,7336,7353,7561];
offt5_10_init = [800,2029,2597,3423,3617,4321,4519,5536,5614,6949,7353,3874];

[tc nontc] = get_time_courses_matrix(datarun000_10, off_otherotherotherother10); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(off_otherotherotherother10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

temp_tcs = get_time_courses_matrix(datarun000_10, off_otherotherotherother10);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherotherother10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherotherother10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherotherother10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);


[C ia ib] = intersect(offt5_10_init, off_otherotherotherother10);
vc = ones(length(off_otherotherotherother10),1);
vc(ib) = 2;


X = [];
X(:,1) = TCParams.maxtim;
X(:,2) = TCParams.mintim;
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherotherotherother10, tcnormnorm,0, vc);
offt5_10 = off_otherotherotherother10(idx==2);
off_otherotherotherotherother10 = off_otherotherotherother10(idx==1);%%
offt5_10 = [800,1459,2029,2597,3423,3617,3874,4321,4519,5536,5614,6949,7353];
off_otherotherotherotherother10 = [17,454,496,1231,1486,1624,1701,1801,1803,1895,2057,2073,2462,2794,3016,3093,3363,4367,4384,4580,4996,5449,5746,5836,5959,6378,6483,6856,7336,7561];

 %% extra cells off t5
 figure();
plot_rf_summaries(datarun000_10, [offt5_10 574 1533 1698 2686 2971 4878 5446 6573 2390 2495], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_10,[offt5_10  574 1533 1698 2686 2971 4878 5446 6573 2495 ], 'all', true, 'bw', true);


%%
plot_rf_summaries(datarun000_10, [offt5_10], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/', 'offt5rf', gcf)

%%
close all;
plot_time_courses(datarun000_10,[offt5_10], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/', 'offt5tc', gcf)





%% off t6

off_otherotherotherotherother10 = [17,454,496,1231,1486,1624,1701,1801,1803,1895,2057,2073,2462,2794,3016,3093,3363,4367,4384,4580,4996,5449,5746,5836,5959,6378,6483,6856,7336,7561];
offt6_10_init = [5836 6856 1486 496];

datarun000_10 = get_interspikeinterval(datarun000_10, off_otherotherotherotherother10);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherotherotherotherother10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, off_otherotherotherotherother10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_10, off_otherotherotherotherother10); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(off_otherotherotherotherother10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

close all;
rstd = [];
meanpix = [];
cellind = get_cell_indices(datarun000_10, off_otherotherotherotherother10);
stamat = cell(length(cellind),1);
 
for i = 1:length(cellind)
    B = [];
    B = datarun000_10.stas.rfs{cellind(i), 1}(:)';
    meanpix(i) = mean(B);
    rstd(i) = robust_std(B, [1]);
    stamat{i,1} = zeros(size(datarun000_10.stas.rfs{cellind(i),1},1),size(datarun000_10.stas.rfs{cellind(i),1},2));
    for j = 1:size(datarun000_10.stas.rfs{cellind(i),1},1)
            for k = 1:size(datarun000_10.stas.rfs{cellind(i),1},2)
                if (datarun000_10.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
                    stamat{i,1}(j,k) = 1;
                end
            end
    end
 
    
end
    
 
 
 
close all;
v = [];
pm = [];
for i = 1:length(cellind)
    rm = [];
    rn = [];
    rmm = [];
    rnn = [];
    DT = [];
    kr = [];
    points = [];
    [rm,rn] = find(stamat{i,1}); %finds all the significant pixels
            if ~isempty(rm)
                [rmm,rnn] = pix_border(rm,rn); %adding everything by .5 and subtracting by 0.5 - i think it is calculating the 4 coordinate values for each pixel - edge correction
                DT = delaunayTriangulation(rmm,rnn); %prediction is still correct with inputs rm and rn, but edges need correcting
                if ~isempty(DT.ConnectivityList)
                    [kr v(1,i)] = convexHull(DT);
                end
                %spy(stamat{i,1}, 'LineSpec', 'r');
                [x,y] = find(stamat{i,1});
                clr = stamat{i,1}(stamat{i,1}~=0);
                %scatter(y,x,20,'MarkerFaceColor',[rand(1) rand(1) rand(1)],'LineWidth',0.05)
                %set(gca,'Xdir','reverse');%'Ydir','reverse')
                %plot(DT.Points(:,2),DT.Points(:,1), '.','markersize',3);
                %hold on;
                %plot(DT.Points(kr,2), DT.Points(kr,1), 'Color',[rand(1) rand(1) rand(1)]);
                %hold off;
                %pause;
                points(:,1) = DT.Points(kr,2);
                points(:,2) = DT.Points(kr,1);
                perimeter = 0;
                for j = 1:size(points, 1)-1
                    perimeter = perimeter + norm(points(j, :) - points(j+1, :));
                end
                perimeter = perimeter + norm(points(end, :) - points(1, :)); % Last point to first
                pm(i) = perimeter;
            else
            end
            v(2,i) = length(rm);
end



 num_rgcs = length(off_otherotherotherotherother10);

% initialize some variables for the look
rf_areas = zeros(num_rgcs,1);
abs_mean_pixel_val = zeros(num_rgcs,1);
snrs = zeros(num_rgcs, 1);
contrast_index = zeros(num_rgcs, 1);
[Y, X] = meshgrid(1:1:40, 1:1:80);


for rgc = 1:num_rgcs
    
    temp_index= get_cell_indices(datarun000_10, off_otherotherotherotherother10(rgc));
    
    %plot_rf(datarun000_10, datarun000_10.cell_ids(rgc), 'sig_stix', true)
    
    [I, J] = find(full(datarun000_10.stas.marks{temp_index}));
    if length(I) < 4
        rf_areas(rgc) = 1;
        warning('report')
    else
        % get the convex hull defined by the significant stixels
        hull_coords = convhull(I,J);
        
        % get the area of the convex hull and save it for every rgc
        rf_areas(rgc) = polyarea(J(hull_coords),I(hull_coords));
        
        % get pixels that are in polygon defined by convex hull
        IN = inpolygon(X, Y, J(hull_coords), I(hull_coords));
        
        % get the number of pixels for averaging
        temp_num_pixels = length(find(IN > 0));

        % get the RF for this rgc
        temp_rf = get_rf(datarun000_10, off_otherotherotherotherother10(rgc));
 
      
        % normalize the rf
        temp_rf = temp_rf ./ abs(ext(reshape(temp_rf, [], 1)));
        
        % mask the RF by the pixels in the polygon 
        masked_rf = temp_rf .* IN';
        
        % number of positive pixels
        num_positive = length(find(masked_rf > 0));
        num_negative = length(find(masked_rf < 0));
        
        % compute contrast index
        contrast_index(rgc) = 1 - abs(((num_positive - num_negative) ./ temp_num_pixels));
    end
end


[COEFF,SCORE] = princomp(isinormnorm');


[C ia ib] = intersect(offt6_10_init, off_otherotherotherotherother10);
vc = ones(length(off_otherotherotherotherother10),1);
vc(ib) = 2;


X = [];
X(:,1) = SCORE(:,1);
X(:,2) = 1 - v(2,:)./v(1,:);
X(:,3) = contrast_index;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherotherotherotherother10, tcnormnorm,0, vc);
offt6_10 = off_otherotherotherotherother10(idx==2);
off_otherotherotherotherotherother10 = off_otherotherotherotherother10(idx==1);%%
offt6_10 = [496,1486,5836,6856];
off_otherotherotherotherotherother10 = [17,454,1231,1624,1701,1801,1803,1895,2057,2073,2462,2794,3016,3093,3363,4367,4384,4580,4996,5449,5746,5959,6378,6483,7336,7561];

 %% extra cells off t6 - same
 figure();
plot_rf_summaries(datarun000_10, [offt6_10 ], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_10,[offt6_10  ], 'all', true, 'bw', true);
%%
plot_rf_summaries(datarun000_10, [offt6_10], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T6/', 'offt6rf', gcf);

close all;
plot_time_courses(datarun000_10,[offt6_10], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T6/', 'offt6tc', gcf)



%%
datarun000_10 = get_interspikeinterval(datarun000_10, offt6_10);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt6_10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, offt6_10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
 maxisi(1,i) = max(isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
maxisi = repmat(maxisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;
isimax = isi./maxisi;

figure();
shadedErrorBar(x2(1:50),mean(isinormnorm(1:50, :)'),std(isinormnorm(1:50, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T6/', 'offt6isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T6/', 'offt6isifull', gcf)
%%
wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); 
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,offt6_10), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:24; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt6_10)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end

b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
%plot(binSize,psthnorm(1,:), 'b');
stairs([0 3 12 15 24],[0.4 0.35 0.3 0.35 0.35], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T6/', 'offt6psth2', gcf)

spikesbytrials{1,1} = get_raster(datarun001_10.spikes{get_cell_indices(datarun001_10, offt6_10(1,1)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 24 0 4]);
hold on;
stairs([0 3 12 15 24],[w g b g g], 'Color', 'r', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/OFF/T6/', 'offt6psthcell496', gcf);
close;


%%
plot_rf_summaries(datarun000_10, [on_otherotherother10], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/', 'onunclassfrf', gcf)

%%
close all;
plot_time_courses(datarun000_10,[on_otherotherother10], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/', 'obunclassftc', gcf)



%%
offt4_31_init =[681,1772,2178,2899,3946,6499,2973,7444,591,4036];

datarun000_31 = get_interspikeinterval(datarun000_31, off_otherother31);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherother31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, off_otherother31(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_31, off_otherother31); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(off_otherother31) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

temp_tcs = get_time_courses_matrix(datarun000_31, off_otherother31);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherother31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherother31)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherother31) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);

 
[C ia ib] = intersect(offt4_31_init, off_otherother31);
vc = ones(length(off_otherother31),1);
vc(ib) = 2;


[COEFF1,SCORE1] = princomp(isinormnorm');

X = [];
X(:,1) = SCORE1(:,1);
X(:,2) = TCParams.minval;
X(:,3) = TCParams.mintim;
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_otherother31, tcnormnorm,0,vc);
 %offt4_31 = off_otherother31(idx==2);
 %off_otherotherotherother31 = off_otherother31(idx==1);
 off_otherother31(idx==2)
ismember(offt4_31_init, off_otherother31(idx==2))

 off_otherotherother31 = off_otherother31(idx==1);
 
 %%
 
 offt3_31_init =[1577,1954,2659,7324,6033,4383,2884];
offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];

offt3_10_init =[512,1368,2750,4188,4487,5732,6170,6511,6588];
offt5_10_init = [800,2597,3617,4321,4519,5536,5614,7353,3423,2029,6949];

offt3_15_init = [62,991,1156,4234,4278,4487,5733,6286,6931];
offt5_15_init = [1246,2253,3695,5116,6260,4771,6380,226];
 

plot_time_courses(datarun000_15, offt3_15_init, 'all', true, 'bw', true);
hold on;
plot_time_courses(datarun000_15, offt5_15_init, 'all', true, 'bw', true);

plot_time_courses(datarun000_15, off_otherotherother15, 'all', true, 'bw', true);
plot_time_courses(datarun000_10, off_otherotherother10, 'all', true, 'bw', true);
plot_time_courses(datarun000_31, off_otherotherother31, 'all', true, 'bw', true);


%%
[tc nontc] = get_time_courses_matrix(datarun000_15, offt3_15_init); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
for i = 1:length(offt3_15_init) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

plot(1:30, tcnormnorm(:,:), 'r');
hold on;

[tc nontc] = get_time_courses_matrix(datarun000_15, offt5_15_init); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
for i = 1:length(offt5_15_init) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

plot(1:30, tcnormnorm(:,:), 'b');

 legend('T3', 'T5');


%%
%%
[tc nontc] = get_time_courses_matrix(datarun000_31, offt5_31_init); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
for i = 1:size(tc,2) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

plot(1:30, tcnormnorm(:,:), 'r');
hold on;

[tc nontc] = get_time_courses_matrix(datarun000_31, off_otherotherother31(~ismember(off_otherotherother31,offt5_31_init))); %or cellids % %offt5_31_init
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
for i = 1:size(tc,2) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

plot(1:30, tcnormnorm(:,:), 'b');

 legend('T5');
 title('31')
 
 %%
 
 %close all;
[tc nontc] = get_time_courses_matrix(datarun000_31, [offt5_31_init]); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
for i = 1:size(tc,2) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;
for i = 1:size(tc,2)
    plot(1:30, tcnormnorm(:,i), 'b');
    pause;
    hold on;
end



%%

[tc nontc] = get_time_courses_matrix(datarun000_10, off_otherotherother10); %or cellids
x = 1:1:30;
 normval = [];
 mx = [];
tcnormnorm = [];
tcnormmx = [];
for i = 1:length(off_otherotherother10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);

tcnormnorm = tc./normval;
[mx mxt] = max(tc);
mx = repmat(mx, size(tc, 1), 1);
tcnormmx = tc./mx;

temp_tcs = get_time_courses_matrix(datarun000_10, off_otherotherother10);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherother10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherother10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherother10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0);

radius = [];
radius = get_rf_fit_radius(datarun000_10, off_otherotherother10);

vc = [];
[C ia ib] = intersect(offt3_10_init, off_otherotherother10);
vc = ones(length(off_otherotherother10),1);
vc(ib) = 2;
%minval maxval zc dot maxt mint

%[COEFF,SCORE] = princomp(tcnormmx');

X = [];
X(:,1) = TCParams.maxmingrad;
X(:,2) = TCParams.maxtim;
X(:,3) = TCParams.dot;
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherotherother10, tcnormnorm,0, vc);
off_otherotherother10(idx==2)
ismember(off_otherotherother10(idx==2),offt3_10_init)

%offt3_10 = off_otherotherother10(idx==2);
%off_otherotherother10 = off_otherotherother10(idx==1);



 
 

%%


datarun001_10.names.rrs_neurons_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-10-1\data003-map\data003-map.neurons';
 datarun001_10 = load_neurons(datarun001_10);
 datarun001_10.names.rrs_params_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-10-1\data003-map\data003-map.params';
datarun001_10 = load_params(datarun001_10);

 
 datarun002_10.names.rrs_neurons_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-10-1\data002-map\data002-map.neurons';
  datarun002_10 = load_neurons(datarun002_10);
  datarun002_10.names.rrs_params_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-10-1\data002-map\data002-map.params';
 datarun002_10 = load_params( datarun002_10);
  datarun002_10.names.stimulus_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-10-1\stimuli\s02';
    datarun002_10 = load_stim(datarun002_10);
    
    
 datarun000_10.names.rrs_neurons_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-10-1\data001-map\data001-map.neurons';
 datarun000_10 = load_neurons(datarun000_10);
 datarun000_10.names.rrs_sta_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-10-1\data001-map\data001-map.sta';
 datarun000_10 = load_sta(datarun000_10, 'load_sta', 'all');
    marks_params.thresh = 4.0;
    datarun000_10 = get_sta_summaries(datarun000_10, 'all', 'marks_params', marks_params);
    %params and get sta fits and ei dont work
 

offt5_10 = [800 2597 3617 4321 4519 5536 5614 7353 3423 2029 6949];

ont3_10 = [1759 3213 4443 4534 6275 6752 7564];


plot_time_courses(datarun000_10,[1759 3213 4443 4534 6275 6752 7564 ], 'all' , true, 'bw', true, 'normalize', true);

datarun000_10 = get_interspikeinterval(datarun000_10, on_otherother10);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(on_otherother10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, on_otherother10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

plot(x2, isinormnorm(:,~ismember(on_otherother10,ont3_10)), 'b')
hold on;
plot(x2, isinormnorm(:,ismember(on_otherother10,[ont3_10])), 'r')
title('ISI-norm - T3 in red')

[tc nontc] = get_time_courses_matrix(datarun000_10, on_otherother10); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(on_otherother10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];


plot(x, tcnormnorm(:,~ismember(on_otherother10,ont3_10)), 'b')
hold on;
plot(x, tcnormnorm(:,ismember(on_otherother10,ont3_10)), 'r')
title('TC-norm - T3 in red')


wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); %For 2012-10-10-1 dataset, that is how the triggers are arranges
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,on_otherother10), 0, '/0', wh, gr, 24, false);
binSize = 0.1:0.1:24; %change depending on length of trial
pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(on_otherother10)
 pulsePSTH(:,a) = sum(nhist{a,1})./15; %change depending on num of trials
end
 
plot(binSize, pulsePSTH(:,~ismember(on_otherother10,ont3_10)), 'b')
hold on;
plot(binSize, pulsePSTH(:,ismember(on_otherother10,ont3_10)), 'r')
title('Pulse PSTH - T4 in red')

for i = 1:length(on_otherother10)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
end
normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;

[maxpulse maxpulsetime] = max(pulsePSTH);
maxpulsetime = maxpulsetime*0.1;

plot(binSize, pulsenormnormPSTH(:,~ismember(on_otherother10,ont3_10)), 'b')
hold on;
plot(binSize, pulsenormnormPSTH(:,ismember(on_otherother10,ont3_10)), 'r')


 [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_10,on_otherother10,0);
 
 NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 32),:));
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
plot(A32(~ismember(on_otherother10,ont3_10)), A256(~ismember(on_otherother10,ont3_10)), 'ob')
hold on;
plot(A32(ismember(on_otherother10,ont3_10)), A256(ismember(on_otherother10,ont3_10)), 'or')






[C ia ib] = intersect(ont3_10, on_otherother10);
vc = ones(length(on_otherother10),1);
vc(ib) = 2;


[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF2,SCORE2] = princomp(tcnormmx');

X = [];
X(:,1) = SCORE2(:,1);
X(:,2) = SCORE2(:,2);
X(:,3) = SCORE2(:,3);
[idx] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, on_otherother10, tcnormnorm,0,vc);
 on_otherother10(idx==2)
 ismember(ont3_10, on_otherother10(idx==2))

 %% ict 8th
 % 

offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376];
%offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];
offt3_31_init =[1577,1954,2659,7324,6033,4383,2884];


cellids = offt5_31_init(7:11);
datarun002 = datarun002_31;
% % 
%[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, datarun002_31.cell_ids, 0);
%[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);


% 
% 
sp = 64;
temp = 256;
% 
spatPer = find(unique(StimComb(:,1))==sp);
tempPer = find(unique(StimComb(:,2))==temp);

    h3 = figure;


for aa = 1:length(cellids)
cell_id = cellids(1, aa);  
figure;
[T R Unew Vnew] = polar_plots_one(rho(:,spatPer),theta(:,spatPer),U(:,spatPer),V(:,spatPer),num,get_cell_indices(datarun002, cell_id));
close;
all = 1:length(datarun002.stimulus.combinations);
[an indic] = sortrows(StimComb, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = D(aa,:); %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(6,9,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 'axis_range', [0 8 0 10], 'stop', 12, 'foa', destaxes, 'tic_color', [0 0 0]);
    b1 = findobj(gcf);
    h4 = findall(b1,'Type','text') ;
    set(h4,'Visible', 'off');
    allaxes=findall(b1,'Type','axes');
    set(allaxes,'Visible', 'off');
end
subplot(6,9,B(aa));
%set(gca, 'FontName', 'Helvetica', 'FontSize', 11)
set(gca, 'FontSize', 9)
p = polar(T{tempPer,1}, R{tempPer,1});
h = findall(gca, 'type', 'line'); % find all lines in the current figure
set(h(h==p), 'LineWidth',3)
h(h==p) = []; % remove the line you ploted from the list.
set(h, 'LineStyle', '--');
b = findobj(gcf);
h2 = findall(b,'Type','text') ;
set(h2,'Visible', 'off');
hold on;
h1 = compass(Unew(tempPer,1),Vnew(tempPer,1), 'r'); %Vector average plot
set(h1,'linewidth',3) 
hold off;
end

save_figure_pdf('/Users/sneharavi/Desktop/', 'T5-2-256-2012-10-31', h3)

%% TRY LOOKING AT DG FOR ALL SPEEDS FOR T3 AND T5

 
 

A(1,:) = [12 3 2 1 10 19 20 21];

B(1) = [11];

D(1,:) = [12 3 2 1 10 19 20 21];
D(2,:) = D(1,:)+3;
D(3,:) = D(2,:)+3;
D(4,:) = D(1,:)+27;
D(5,:) = D(4,:)+3;
D(6,:) = D(5,:)+3;
D(7,:) = D(4,:)+27;
D(8,:) = D(7,:)+3;
D(9,:) = D(8,:)+3;
D(10,:) = D(7,:)+27;
D(11,:) = D(10,:)+3;
D(12,:) = D(11,:)+3;



plot_rf_portraits(datarun000_31,offt3_31_init)
plot_rf_portraits(datarun000_31,offt5_31_init)
 
 



%%
%%
off_otherotherother10 = [17,454,496,1231,1459,1486,1624,1701,1801,1803,1895,2057,2073,2462,2794,3016,3093,3363,4367,4384,4580,4996,5449,5746,5836,5959,6378,6483,6856,7336,7561];
offt6_10_init = [5836 6856 1486 496]; %17];


datarun000_10 = get_interspikeinterval(datarun000_10, off_otherotherother10);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherotherother10) %or nonds
 isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, off_otherotherother10(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_10, off_otherotherother10); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(off_otherotherother10) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

temp_tcs = get_time_courses_matrix(datarun000_10, off_otherotherother10);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherother10)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherother10)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherother10) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);

wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,off_otherotherother10), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:24; %change depending on length of trial
pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(off_otherotherother10)
 pulsePSTH(:,a) = sum(nhist{a,1})./15; %change depending on num of trials
end
 
for i = 1:length(off_otherotherother10)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
end
normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;

[maxpulse maxpulsetime] = max(pulsePSTH);
maxpulsetime = maxpulsetime*0.1;

 
[C ia ib] = intersect(offt6_10_init, off_otherotherother10);
vc = ones(length(off_otherotherother10),1);
vc(ib) = 2;


sumsumsp = [];
for i = 1: length(sumSpTrTrig)
    sumsumsp(i,:) = sum(sumSpTrTrig{i, 1});
end
  
%% Calculate sta parameter
% marks_params.thresh = 4.0;
% datarun000_10 = get_significant_stixels(datarun000_10, off_otherotherother10, 'thresh_params', marks_params);
% 
% for i = 1:length(offt6_10_init)
%     plot_rf(datarun000_10,offt6_10_init(i));
%     hold on;
%     cell_index = get_cell_indices(datarun000_10, offt6_10_init(i));
%     sig_stix = datarun000_10.stas.marks{cell_index};
%     temp_indices = find(full(sig_stix));
%     [temp_rows, temp_cols] = ind2sub(size(sig_stix), temp_indices);
%     plot(temp_cols, temp_rows, 'ro', 'MarkerSize', 7);
%     hold off;
%      pause;
%      close all
% end

%%
close all;
rstd = [];
meanpix = [];
cellind = get_cell_indices(datarun000_10, off_otherotherother10);
stamat = cell(length(cellind),1);
 
for i = 1:length(cellind)
    B = [];
    B = datarun000_10.stas.rfs{cellind(i), 1}(:)';
    meanpix(i) = mean(B);
    rstd(i) = robust_std(B, [1]);
    stamat{i,1} = zeros(size(datarun000_10.stas.rfs{cellind(i),1},1),size(datarun000_10.stas.rfs{cellind(i),1},2));
    for j = 1:size(datarun000_10.stas.rfs{cellind(i),1},1)
            for k = 1:size(datarun000_10.stas.rfs{cellind(i),1},2)
                if (datarun000_10.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 4.6*rstd(i))
                    stamat{i,1}(j,k) = 1;
                end
            end
    end
 
    
end
    
 
 
 
close all;
v = [];
pm = [];
for i = 1:length(cellind)
    rm = [];
    rn = [];
    rmm = [];
    rnn = [];
    DT = [];
    kr = [];
    points = [];
    [rm,rn] = find(stamat{i,1}); %finds all the significant pixels
            if ~isempty(rm)
                [rmm,rnn] = pix_border(rm,rn); %adding everything by .5 and subtracting by 0.5 - i think it is calculating the 4 coordinate values for each pixel - edge correction
                DT = delaunayTriangulation(rmm,rnn); %prediction is still correct with inputs rm and rn, but edges need correcting
                if ~isempty(DT.ConnectivityList)
                    [kr v(1,i)] = convexHull(DT);
                end
                %plot_rf(datarun000_10,off_otherotherother10(i));
                %spy(stamat{i,1})%'LineSpec', 'or');
                [x,y] = find(stamat{i,1});
                clr = stamat{i,1}(stamat{i,1}~=0);
                %scatter(y,x,20,'MarkerFaceColor',[1 0 0],'LineWidth',0.05)
                %set(gca,'Xdir','reverse');%'Ydir','reverse')
                %plot(DT.Points(:,2),DT.Points(:,1), '.','markersize',3);
                %hold on;
                %plot(DT.Points(kr,2), DT.Points(kr,1), 'Color',[1 0 0]);%[rand(1) rand(1) rand(1)]);
                %hold off;
                %pause;
                points(:,1) = DT.Points(kr,2);
                points(:,2) = DT.Points(kr,1);
                perimeter = 0;
                for j = 1:size(points, 1)-1
                    perimeter = perimeter + norm(points(j, :) - points(j+1, :));
                end
                perimeter = perimeter + norm(points(end, :) - points(1, :)); % Last point to first
                pm(i) = perimeter;
            else
            end
            v(2,i) = length(rm);
end
%%
% datarun000_10 = get_significant_stixels(datarun000_10, off_otherotherother10);
% stavar = [];
% for i = 1:length(off_otherotherother10)
%         ii = get_cell_indices(datarun000_10, off_otherotherother10(1,i));
%         stavar(i,1) = var(datarun000_10.stas.rfs{ii,1}(datarun000_10.stas.significant_stixels{ii, 1}));
%         if(isnan(stavar(i,1)))
%             stavar(i,1) = 0;
%         end
% end

%%
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_10, off_otherotherother10, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V spave] = dscellanalysis(NumSpikesCell, StimComb);

%%

minfiringslow = [];
maxfiringslow = [];
minfiringfast = [];
maxfiringfast = [];

for i = 1:length(off_otherotherother10)
minfiringslow(i) =  (spave{2,1}(i,1)+spave{2,1}(i,2))/2;
maxfiringslow(i) =  (spave{2,1}(i,5)+spave{2,1}(i,8))/2;
minfiringfast(i) = (spave{1,1}(i,3)+spave{1,1}(i,5))/2;
maxfiringfast(i) = (spave{1,1}(i,2)+spave{1,1}(i,7))/2;
end
%%
scatter(minfiringslow(1,ismember(off_otherotherother10,offt6_10_init)),maxfiringslow(1,ismember(off_otherotherother10,offt6_10_init)),'r');
hold on;
scatter(minfiringslow(1,~ismember(off_otherotherother10,offt6_10_init)),maxfiringslow(1,~ismember(off_otherotherother10,offt6_10_init)),'b');
xlabel('Min Firing - Slow Speed')
ylabel('Max Firing - Slow Speed')
pause;
hold off;


scatter(minfiringfast(1,ismember(off_otherotherother10,offt6_10_init)),maxfiringfast(1,ismember(off_otherotherother10,offt6_10_init)),'r');
hold on;
scatter(minfiringfast(1,~ismember(off_otherotherother10,offt6_10_init)),maxfiringfast(1,~ismember(off_otherotherother10,offt6_10_init)),'b');
xlabel('Min Firing - Fast Speed')
ylabel('Max Firing - Fast Speed')
pause;
hold off;

xx = maxfiringslow./minfiringslow;
yy = maxfiringfast./minfiringfast;


scatter(xx(1,ismember(off_otherotherother10,offt6_10_init)),yy(1,ismember(off_otherotherother10,offt6_10_init)),'r');
hold on;
scatter(xx(1,~ismember(off_otherotherother10,offt6_10_init)),yy(1,~ismember(off_otherotherother10,offt6_10_init)),'b');
xlabel('Max/Min Firing - Slow Speed')
ylabel('Max/Min Firing - Fast Speed')
pause;
hold off;


xx = 1-(minfiringslow./maxfiringslow);
xx(minfiringslow./maxfiringslow > 1) = 0;
yy = 1-(minfiringfast./maxfiringfast);
yy(minfiringfast./maxfiringfast > 1) = 0;

scatter(xx(1,ismember(off_otherotherother10,offt6_10_init)),yy(1,ismember(off_otherotherother10,offt6_10_init)),'r');
hold on;
scatter(xx(1,~ismember(off_otherotherother10,offt6_10_init)),yy(1,~ismember(off_otherotherother10,offt6_10_init)),'b');
xlabel('1-Min/Max Firing - Slow Speed')
ylabel('1-Min/Max Firing - Fast Speed')
pause;
hold off;

% minangle = [];
% maxangle = [];
% minAxis = [];
% maxAxis = [];
% minAxFir = [];
% maxAxFir = [];
% minFir = cell(2,1);
% maxFir = cell(2,1);
% [CC,minAxis] = min(spave{2,1}');
% for i = 1:length(off_otherotherother10)
%     minangle(1,i) = theta{2,1}(i,minAxis(1,i));
%     minangle(2,i) = minangle(1,i)+pi;
%     if (minangle(2,i) < 0 || minangle(2,i) > pi + 3*pi/4)
%         minangle(2,i) = minangle(1,i) - pi;
%     end
%     minAxis(2,i) = find(minangle(2,i) == theta{2,1}(1,:));
%     maxangle(1,i) = minangle(1,i)+pi./2;
%     if (maxangle(1,i) < 0 || maxangle(1,i) > pi + 3*pi/4)
%         maxangle(1,i) = minangle(1,i)-pi./2;
%     end
%     maxAxis(1,i) = find(maxangle(1,i) == theta{2,1}(1,:));
%     maxangle(2,i) = maxangle(1,i)+pi;
%         if (maxangle(2,i) < 0 || maxangle(2,i) > pi + 3*pi/4)
%                 maxangle(2,i) = maxangle(1,i)-pi;
%         end
%       maxAxis(2,i) = find(maxangle(2,i) == theta{2,1}(1,:));
%       minAxFir(i) = (spave{2,1}(i,minAxis(1,i))+spave{2,1}(i,minAxis(2,i)))/2;
%       maxAxFir(i) = (spave{2,1}(i,maxAxis(1,i))+spave{2,1}(i,maxAxis(2,i)))/2;
% end

% scatter(minAxFir,maxAxFir);
% hold on;
% xlabel('minimum firing');
% ylabel('maximum firing')
% title('S 64 T 256')
% %plot(0:1:250, 0:1:250);
% hold off;
% figure()
% hist(maxAxFir./minAxFir, 20)
% xlabel('maximum firing / minimum firing');
% title('S 64 T 256')

% minFir{2,1} = minAxFir;
% maxFir{2,1} = maxAxFir;
% 
% 
% minangle = [];
% maxangle = [];
% minAxis = [];
% maxAxis = [];
% minAxFir = [];
% maxAxFir = [];
% [CC,minAxis] = min(spave{1,1}')
% for i = 1:length(off_otherotherother10)
%     minangle(1,i) = theta{1,1}(i,minAxis(i));
%     minangle(2,i) = minangle(1,i)+pi;
%     if (minangle(2,i) < 0 || minangle(2,i) > pi + 3*pi/4)
%         minangle(2,i) = minangle(1,i) - pi;
%     end
%     minAxis(2,i) = find(minangle(2,i) == theta{1,1}(1,:));
%     maxangle(1,i) = minangle(1,i)+pi./2;
%     if (maxangle(1,i) < 0 || maxangle(1,i) > pi + 3*pi/4)
%         maxangle(1,i) = minangle(1,i)-pi./2;
%     end
%     maxAxis(1,i) = find(maxangle(1,i) == theta{1,1}(1,:));
%     maxangle(2,i) = maxangle(1,i)+pi;
%         if (maxangle(2,i) < 0 || maxangle(2,i) > pi + 3*pi/4)
%                 maxangle(2,i) = maxangle(1,i)-pi;
%         end
%       maxAxis(2,i) = find(maxangle(2,i) == theta{1,1}(1,:));
%       minAxFir(i) = (spave{1,1}(i,minAxis(1,i))+spave{1,1}(i,minAxis(2,i)))/2;
%       maxAxFir(i) = (spave{1,1}(i,maxAxis(1,i))+spave{1,1}(i,maxAxis(2,i)))/2;
% end
% 
% % scatter(minAxFir,maxAxFir);
% % hold on;
% % xlabel('minimum firing');
% % ylabel('maximum firing')
% % title('S 64 T 32')
% % %plot(0:1:250, 0:1:250);
% % hold off;
% % figure();
% % hist(maxAxFir./minAxFir,20)
% % xlabel('maximum firing / minimum firing');
% % title('S 64 T 32')
% 
% minFir{1,1} = minAxFir;
% maxFir{1,1} = maxAxFir;
% 
% xxx = 1 - (minFir{1,1}./maxFir{1,1});
% xxx(minFir{1,1}./maxFir{1,1} > 1) = 0;
% yyy = 1 - (minFir{2,1}./maxFir{2,1});
% yyy(minFir{2,1}./maxFir{2,1} > 1) = 0;
% scatter(xxx,yyy);
% xlabel('maximum firing / minimum firing - T 32');
% ylabel('maximum firing / minimum firing - T 256')
% title('S 64 T 32-256')
% figure()
% hist(xxx.*yyy,100);
% xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');
%%


       num_rgcs = length(off_otherotherother10);

% initialize some variables for the look
rf_areas = zeros(num_rgcs,1);
abs_mean_pixel_val = zeros(num_rgcs,1);
snrs = zeros(num_rgcs, 1);
contrast_index = zeros(num_rgcs, 1);
[Y, X] = meshgrid(1:1:40, 1:1:80);


for rgc = 1:num_rgcs
    
    temp_index= get_cell_indices(datarun000_10, off_otherotherother10(rgc));
    
    %plot_rf(datarun000_10, datarun000_10.cell_ids(rgc), 'sig_stix', true)
    
    [I, J] = find(full(datarun000_10.stas.marks{temp_index}));
    if length(I) < 4
        rf_areas(rgc) = 1;
        warning('report')
    else
        % get the convex hull defined by the significant stixels
        hull_coords = convhull(I,J);
        
        % get the area of the convex hull and save it for every rgc
        rf_areas(rgc) = polyarea(J(hull_coords),I(hull_coords));
        
        % get pixels that are in polygon defined by convex hull
        IN = inpolygon(X, Y, J(hull_coords), I(hull_coords));
        
        % get the number of pixels for averaging
        temp_num_pixels = length(find(IN > 0));

        % get the RF for this rgc
        temp_rf = get_rf(datarun000_10, off_otherotherother10(rgc));
 
      
        % normalize the rf
        temp_rf = temp_rf ./ abs(ext(reshape(temp_rf, [], 1)));
        
        % mask the RF by the pixels in the polygon 
        masked_rf = temp_rf .* IN';
        
        % number of positive pixels
        num_positive = length(find(masked_rf > 0));
        num_negative = length(find(masked_rf < 0));
        
        % compute contrast index
        contrast_index(rgc) = 1 - abs(((num_positive - num_negative) ./ temp_num_pixels));
    end
end
%%

[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsenormnormPSTH');
%[COEFF2,SCORE2] = princomp(tcnormminn');

%%
[COEFF1,SCORE1] = princomp(isinormnorm');

X = [];

X(:,1) = contrast_index;%(sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = 1- (v(2,:)./v(1,:));%SCORE1(:,1);
scatter(X(ismember(off_otherotherother10,offt6_10_init),1), X(ismember(off_otherotherother10,offt6_10_init),2), 'r');
hold on;
scatter(X(~ismember(off_otherotherother10,offt6_10_init),1), X(~ismember(off_otherotherother10,offt6_10_init),2), 'b');
xlabel('CI');
ylabel('1-ratio');
% zlabel('CI');
hold off;
pause;


X = [];
X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
scatter(X(ismember(off_otherotherother10,offt6_10_init),1), X(ismember(off_otherotherother10,offt6_10_init),2), 'r');
hold on;
scatter(X(~ismember(off_otherotherother10,offt6_10_init),1), X(~ismember(off_otherotherother10,offt6_10_init),2), 'b');
xlabel('Pulse');
ylabel('ISI');
hold off;
pause;
%%

% pulse isi 1-ratio CI 
% (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));SCORE1(:,1); 1- (v(2,:)./v(1,:)); contrast_index;

[COEFF1,SCORE1] = princomp(isinormnorm');

X = [];

X(:,1) = SCORE1(:,1);
X(:,2) = 1- (v(2,:)./v(1,:));
X(:,3) = contrast_index;
scatter3(X(ismember(off_otherotherother10,offt6_10_init),1), X(ismember(off_otherotherother10,offt6_10_init),2), X(ismember(off_otherotherother10,offt6_10_init),3),'r');
hold on;
scatter3(X(~ismember(off_otherotherother10,offt6_10_init),1), X(~ismember(off_otherotherother10,offt6_10_init),2),  X(~ismember(off_otherotherother10,offt6_10_init),3),'b');

%%
X = [];

X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
X(:,3) = contrast_index;
scatter3 (X(ismember(off_otherotherother10,offt6_10_init),1), X(ismember(off_otherotherother10,offt6_10_init),2), X(ismember(off_otherotherother10,offt6_10_init),3), 'r');
hold on;
scatter3 (X(~ismember(off_otherotherother10,offt6_10_init),1), X(~ismember(off_otherotherother10,offt6_10_init),2), X(~ismember(off_otherotherother10,offt6_10_init),3), 'b');
xlabel('pulse');
ylabel('isi');
zlabel('CI');

%%
vcc  =[];
vcc(1,:) = mean(X(ib,:));
vcc(2,:) = mean(X(vc==1,:));
[idx] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherotherother10, tcnormnorm,0,vc);
 %offt4_10 = off_otherotherother10(idx==2);
 %off_otherotherother10 = off_otherotherother10(idx==1);
 off_otherotherother10(idx==2)
ismember(offt6_10_init, off_otherotherother10(idx==2))
%%
plot(tcnormnorm(:,(vc==1)), 'b')
hold on;
pause
plot(tcnormnorm(:,(vc==2)), 'r')
hold off;
pause;
plot(isinormnorm(:,(vc==1)), 'b')
hold on;
pause
plot(isinormnorm(:,(vc==2)), 'r')
hold off;
pause
plot(pulsenormnormPSTH(:,(vc==1)), 'b')
hold on;
pause
plot(pulsenormnormPSTH(:,(vc==2)), 'r')

%%

%%
 off_otherotherother15 = [226,407,649,843,872,1130,1411,1427,1564,1817,1921,2206,2236,2326,2343,2373,2521,2555,2732,2794,3091,3286,3857,4097,4235,4697,4732,4771,5148,5672,6139,6257,6376,6380,6886,6992,7517];
offt6_15_init = [1817 2343 3286 6992];


datarun000_15 = get_interspikeinterval(datarun000_15, off_otherotherother15);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherotherother15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, off_otherotherother15(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherotherother15); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(off_otherotherother15) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

temp_tcs = get_time_courses_matrix(datarun000_15, off_otherotherother15);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherother15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherother15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherother15) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);

wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,off_otherotherother15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; %change depending on length of trial
pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(off_otherotherother15)
 pulsePSTH(:,a) = sum(nhist{a,1})./25; %change depending on num of trials
end
 
for i = 1:length(off_otherotherother15)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
end
normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;

[maxpulse maxpulsetime] = max(pulsePSTH);
maxpulsetime = maxpulsetime*0.1;

 
[C ia ib] = intersect(offt6_15_init, off_otherotherother15);
vc = ones(length(off_otherotherother15),1);
vc(ib) = 2;


sumsumsp = [];
for i = 1: length(sumSpTrTrig)
    sumsumsp(i,:) = sum(sumSpTrTrig{i, 1});
end
  
%% Calculate sta parameter
% marks_params.thresh = 4.0;
% datarun000_15 = get_significant_stixels(datarun000_15, off_otherotherother15, 'thresh_params', marks_params);
% 
% for i = 1:length(offt6_15_init)
%     plot_rf(datarun000_15,offt6_15_init(i));
%     hold on;
%     cell_index = get_cell_indices(datarun000_15, offt6_15_init(i));
%     sig_stix = datarun000_15.stas.marks{cell_index};
%     temp_indices = find(full(sig_stix));
%     [temp_rows, temp_cols] = ind2sub(size(sig_stix), temp_indices);
%     plot(temp_cols, temp_rows, 'ro', 'MarkerSize', 7);
%     hold off;
%      pause;
%      close all
% end

%%
close all;
rstd = [];
meanpix = [];
cellind = get_cell_indices(datarun000_15, off_otherotherother15);
stamat = cell(length(cellind),1);
 
for i = 1:length(cellind)
    B = [];
    B = datarun000_15.stas.rfs{cellind(i), 1}(:)';
    meanpix(i) = mean(B);
    rstd(i) = robust_std(B, [1]);
    stamat{i,1} = zeros(size(datarun000_15.stas.rfs{cellind(i),1},1),size(datarun000_15.stas.rfs{cellind(i),1},2));
    for j = 1:size(datarun000_15.stas.rfs{cellind(i),1},1)
            for k = 1:size(datarun000_15.stas.rfs{cellind(i),1},2)
                if (datarun000_15.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 4.6*rstd(i))
                    stamat{i,1}(j,k) = 1;
                end
            end
    end
 
    
end
    
 
 
 
close all;
v = [];
pm = [];
for i = 1:length(cellind)
    rm = [];
    rn = [];
    rmm = [];
    rnn = [];
    DT = [];
    kr = [];
    points = [];
    [rm,rn] = find(stamat{i,1}); %finds all the significant pixels
            if ~isempty(rm)
                [rmm,rnn] = pix_border(rm,rn); %adding everything by .5 and subtracting by 0.5 - i think it is calculating the 4 coordinate values for each pixel - edge correction
                DT = delaunayTriangulation(rmm,rnn); %prediction is still correct with inputs rm and rn, but edges need correcting
                if ~isempty(DT.ConnectivityList)
                    [kr v(1,i)] = convexHull(DT);
                end
                %plot_rf(datarun000_15,off_otherotherother15(i));
                %spy(stamat{i,1})%'LineSpec', 'or');
                [x,y] = find(stamat{i,1});
                clr = stamat{i,1}(stamat{i,1}~=0);
                %scatter(y,x,20,'MarkerFaceColor',[1 0 0],'LineWidth',0.05)
                %set(gca,'Xdir','reverse');%'Ydir','reverse')
                %plot(DT.Points(:,2),DT.Points(:,1), '.','markersize',3);
                %hold on;
                %plot(DT.Points(kr,2), DT.Points(kr,1), 'Color',[1 0 0]);%[rand(1) rand(1) rand(1)]);
                %hold off;
                %pause;
                points(:,1) = DT.Points(kr,2);
                points(:,2) = DT.Points(kr,1);
                perimeter = 0;
                for j = 1:size(points, 1)-1
                    perimeter = perimeter + norm(points(j, :) - points(j+1, :));
                end
                perimeter = perimeter + norm(points(end, :) - points(1, :)); % Last point to first
                pm(i) = perimeter;
            else
            end
            v(2,i) = length(rm);
end

% radius = [];
% radius = get_rf_fit_radius(datarun000_15, off_otherotherother15);
%%
% 
% datarun000_15 = get_significant_stixels(datarun000_15, off_otherotherother15);
% stavar = [];
% for i = 1:length(off_otherotherother15)
%         ii = get_cell_indices(datarun000_15, off_otherotherother15(1,i));
%         stavar(i,1) = var(datarun000_15.stas.rfs{ii,1}(datarun000_15.stas.significant_stixels{ii, 1}));
%         if(isnan(stavar(i,1)))
%             stavar(i,1) = 0;
%         end
% end

%%

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, off_otherotherother15, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V spave] = dscellanalysis(NumSpikesCell, StimComb);

%%

minfiringslow = [];
maxfiringslow = [];
minfiringfast = [];
maxfiringfast = [];

for i = 1:length(off_otherotherother15)
minfiringslow(i) =  (spave{2,1}(i,3)+spave{2,1}(i,7))/2;
maxfiringslow(i) =  (spave{2,1}(i,5)+spave{2,1}(i,1))/2;
minfiringfast(i) = (spave{1,1}(i,2)+spave{1,1}(i,3))/2;
maxfiringfast(i) = (spave{1,1}(i,5)+spave{1,1}(i,6))/2;
end
%%
scatter(minfiringslow(1,ismember(off_otherotherother15,offt6_15_init)),maxfiringslow(1,ismember(off_otherotherother15,offt6_15_init)),'r');
hold on;
scatter(minfiringslow(1,~ismember(off_otherotherother15,offt6_15_init)),maxfiringslow(1,~ismember(off_otherotherother15,offt6_15_init)),'b');
xlabel('Min Firing - Slow Speed')
ylabel('Max Firing - Slow Speed')
pause;
hold off;


scatter(minfiringfast(1,ismember(off_otherotherother15,offt6_15_init)),maxfiringfast(1,ismember(off_otherotherother15,offt6_15_init)),'r');
hold on;
scatter(minfiringfast(1,~ismember(off_otherotherother15,offt6_15_init)),maxfiringfast(1,~ismember(off_otherotherother15,offt6_15_init)),'b');
xlabel('Min Firing - Fast Speed')
ylabel('Max Firing - Fast Speed')
pause;
hold off;

xx = maxfiringslow./minfiringslow;
yy = maxfiringfast./minfiringfast;


scatter(xx(1,ismember(off_otherotherother15,offt6_15_init)),yy(1,ismember(off_otherotherother15,offt6_15_init)),'r');
hold on;
scatter(xx(1,~ismember(off_otherotherother15,offt6_15_init)),yy(1,~ismember(off_otherotherother15,offt6_15_init)),'b');
xlabel('Max/Min Firing - Slow Speed')
ylabel('Max/Min Firing - Fast Speed')
pause;
hold off;


xx = 1-(minfiringslow./maxfiringslow);
xx(minfiringslow./maxfiringslow > 1) = 0;
yy = 1-(minfiringfast./maxfiringfast);
yy(minfiringfast./maxfiringfast > 1) = 0;

scatter(xx(1,ismember(off_otherotherother15,offt6_15_init)),yy(1,ismember(off_otherotherother15,offt6_15_init)),'r');
hold on;
scatter(xx(1,~ismember(off_otherotherother15,offt6_15_init)),yy(1,~ismember(off_otherotherother15,offt6_15_init)),'b');
xlabel('1-Min/Max Firing - Slow Speed')
ylabel('1-Min/Max Firing - Fast Speed')
pause;
hold off;

% minangle = [];
% maxangle = [];
% minAxis = [];
% maxAxis = [];
% minAxFir = [];
% maxAxFir = [];
% minFir = cell(2,1);
% maxFir = cell(2,1);
% [CC,minAxis] = min(spave{2,1}');
% for i = 1:length(off_otherotherother15)
%     minangle(1,i) = theta{2,1}(i,minAxis(1,i));
%     minangle(2,i) = minangle(1,i)+pi;
%     if (minangle(2,i) < 0 || minangle(2,i) > pi + 3*pi/4)
%         minangle(2,i) = minangle(1,i) - pi;
%     end
%     minAxis(2,i) = find(minangle(2,i) == theta{2,1}(1,:));
%     maxangle(1,i) = minangle(1,i)+pi./2;
%     if (maxangle(1,i) < 0 || maxangle(1,i) > pi + 3*pi/4)
%         maxangle(1,i) = minangle(1,i)-pi./2;
%     end
%     maxAxis(1,i) = find(maxangle(1,i) == theta{2,1}(1,:));
%     maxangle(2,i) = maxangle(1,i)+pi;
%         if (maxangle(2,i) < 0 || maxangle(2,i) > pi + 3*pi/4)
%                 maxangle(2,i) = maxangle(1,i)-pi;
%         end
%       maxAxis(2,i) = find(maxangle(2,i) == theta{2,1}(1,:));
%       minAxFir(i) = (spave{2,1}(i,minAxis(1,i))+spave{2,1}(i,minAxis(2,i)))/2;
%       maxAxFir(i) = (spave{2,1}(i,maxAxis(1,i))+spave{2,1}(i,maxAxis(2,i)))/2;
% end
% 
% % scatter(minAxFir,maxAxFir);
% % hold on;
% % xlabel('minimum firing');
% % ylabel('maximum firing')
% % title('S 64 T 256')
% % %plot(0:1:250, 0:1:250);
% % hold off;
% % figure()
% % hist(maxAxFir./minAxFir, 20)
% % xlabel('maximum firing / minimum firing');
% % title('S 64 T 256')
% 
% minFir{2,1} = minAxFir;
% maxFir{2,1} = maxAxFir;
% 
% 
% minangle = [];
% maxangle = [];
% minAxis = [];
% maxAxis = [];
% minAxFir = [];
% maxAxFir = [];
% [CC,minAxis] = min(spave{1,1}')
% for i = 1:length(off_otherotherother15)
%     minangle(1,i) = theta{1,1}(i,minAxis(i));
%     minangle(2,i) = minangle(1,i)+pi;
%     if (minangle(2,i) < 0 || minangle(2,i) > pi + 3*pi/4)
%         minangle(2,i) = minangle(1,i) - pi;
%     end
%     minAxis(2,i) = find(minangle(2,i) == theta{1,1}(1,:));
%     maxangle(1,i) = minangle(1,i)+pi./2;
%     if (maxangle(1,i) < 0 || maxangle(1,i) > pi + 3*pi/4)
%         maxangle(1,i) = minangle(1,i)-pi./2;
%     end
%     maxAxis(1,i) = find(maxangle(1,i) == theta{1,1}(1,:));
%     maxangle(2,i) = maxangle(1,i)+pi;
%         if (maxangle(2,i) < 0 || maxangle(2,i) > pi + 3*pi/4)
%                 maxangle(2,i) = maxangle(1,i)-pi;
%         end
%       maxAxis(2,i) = find(maxangle(2,i) == theta{1,1}(1,:));
%       minAxFir(i) = (spave{1,1}(i,minAxis(1,i))+spave{1,1}(i,minAxis(2,i)))/2;
%       maxAxFir(i) = (spave{1,1}(i,maxAxis(1,i))+spave{1,1}(i,maxAxis(2,i)))/2;
% end
% 
% % scatter(minAxFir,maxAxFir);
% % hold on;
% % xlabel('minimum firing');
% % ylabel('maximum firing')
% % title('S 64 T 32')
% % %plot(0:1:250, 0:1:250);
% % hold off;
% % figure();
% % hist(maxAxFir./minAxFir,20)
% % xlabel('maximum firing / minimum firing');
% % title('S 64 T 32')
% 
% minFir{1,1} = minAxFir;
% maxFir{1,1} = maxAxFir;
% 
% xxx = 1 - (minFir{1,1}./maxFir{1,1});
% xxx(minFir{1,1}./maxFir{1,1} > 1) = 0;
% yyy = 1 - (minFir{2,1}./maxFir{2,1});
% yyy(minFir{2,1}./maxFir{2,1} > 1) = 0;
% scatter(xxx,yyy);
% xlabel('maximum firing / minimum firing - T 32');
% ylabel('maximum firing / minimum firing - T 256')
% title('S 64 T 32-256')
% figure()
% hist(xxx.*yyy,100);
% xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');

%%


       num_rgcs = length(off_otherotherother15);

% initialize some variables for the look
rf_areas = zeros(num_rgcs,1);
abs_mean_pixel_val = zeros(num_rgcs,1);
snrs = zeros(num_rgcs, 1);
contrast_index = zeros(num_rgcs, 1);
[Y, X] = meshgrid(1:1:40, 1:1:80);


for rgc = 1:num_rgcs
    
    temp_index= get_cell_indices(datarun000_15, off_otherotherother15(rgc));
    
    %plot_rf(datarun000_15, datarun000_15.cell_ids(rgc), 'sig_stix', true)
    
    [I, J] = find(full(datarun000_15.stas.marks{temp_index}));
    if length(I) < 4
        rf_areas(rgc) = 1;
        warning('report')
    else
        % get the convex hull defined by the significant stixels
        hull_coords = convhull(I,J);
        
        % get the area of the convex hull and save it for every rgc
        rf_areas(rgc) = polyarea(J(hull_coords),I(hull_coords));
        
        % get pixels that are in polygon defined by convex hull
        IN = inpolygon(X, Y, J(hull_coords), I(hull_coords));
        
        % get the number of pixels for averaging
        temp_num_pixels = length(find(IN > 0));

        % get the RF for this rgc
        temp_rf = get_rf(datarun000_15, off_otherotherother15(rgc));
 
      
        % normalize the rf
        temp_rf = temp_rf ./ abs(ext(reshape(temp_rf, [], 1)));
        
        % mask the RF by the pixels in the polygon 
        masked_rf = temp_rf .* IN';
        
        % number of positive pixels
        num_positive = length(find(masked_rf > 0));
        num_negative = length(find(masked_rf < 0));
        
        % compute contrast index
        contrast_index(rgc) = 1 - abs(((num_positive - num_negative) ./ temp_num_pixels));
    end
end
%%
[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsenormnormPSTH');
%[COEFF2,SCORE2] = princomp(tcnormminn');


%%

% pulse isi 1-ratio CI 
% (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));SCORE1(:,1); 1- (v(2,:)./v(1,:)); contrast_index;

[COEFF1,SCORE1] = princomp(isinormnorm');

X = [];

X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = 1- (v(2,:)./v(1,:));
X(:,3) = contrast_index;
scatter3(X(ismember(off_otherotherother15,offt6_15_init),1), X(ismember(off_otherotherother15,offt6_15_init),2), X(ismember(off_otherotherother15,offt6_15_init),3),'r');
hold on;
scatter3(X(~ismember(off_otherotherother15,offt6_15_init),1), X(~ismember(off_otherotherother15,offt6_15_init),2),  X(~ismember(off_otherotherother15,offt6_15_init),3),'b');


%%
[COEFF1,SCORE1] = princomp(isinormnorm');

X = [];

X(:,1) = contrast_index;%(sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = 1- (v(2,:)./v(1,:));%SCORE1(:,1);
scatter(X(ismember(off_otherotherother15,offt6_15_init),1), X(ismember(off_otherotherother15,offt6_15_init),2), 'r');
hold on;
scatter(X(~ismember(off_otherotherother15,offt6_15_init),1), X(~ismember(off_otherotherother15,offt6_15_init),2), 'b');
xlabel('CI');
ylabel('1-ratio');
% zlabel('CI');
hold off;
pause;


X = [];
X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
scatter(X(ismember(off_otherotherother15,offt6_15_init),1), X(ismember(off_otherotherother15,offt6_15_init),2), 'r');
hold on;
scatter(X(~ismember(off_otherotherother15,offt6_15_init),1), X(~ismember(off_otherotherother15,offt6_15_init),2), 'b');
xlabel('Pulse');
ylabel('ISI');
hold off;
pause;

%%

X = [];

X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
X(:,3) = contrast_index;
scatter3 (X(ismember(off_otherotherother15,offt6_15_init),1), X(ismember(off_otherotherother15,offt6_15_init),2), X(ismember(off_otherotherother15,offt6_15_init),3), 'r');
hold on;
scatter3 (X(~ismember(off_otherotherother15,offt6_15_init),1), X(~ismember(off_otherotherother15,offt6_15_init),2), X(~ismember(off_otherotherother15,offt6_15_init),3), 'b');
xlabel('pulse');
ylabel('isi');
zlabel('CI');


%%

[idx] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherotherother15, tcnormnorm,0,vc);
 %offt4_15 = off_otherotherother15(idx==2);
%  %off_otherotherother15 = off_otherotherother15(idx==1);
%  off_otherotherother15(idx==2)
% ismember(offt6_15_init, off_otherotherother15(idx==2))


%%
vcc  =[];
vcc(1,:) = mean(X(ib,:));
vcc(2,:) = mean(X(vc==1,:));
[idx] = clustering_analysis_plots(X, 1,0, 2, 1, 0, datarun000_15, off_otherotherother15, tcnormnorm,0,vcc);
%%
scatter3 (X(:,1), X(:,2), X(:,3))

%%
plot(tcnormnorm(:,(vc==1)), 'b')
hold on;
pause
plot(tcnormnorm(:,(vc==2)), 'r')
hold off;
pause;
plot(isinormnorm(:,(vc==1)), 'b')
hold on;
pause
plot(isinormnorm(:,(vc==2)), 'r')
hold off;
pause
plot(pulsenormnormPSTH(:,(vc==1)), 'b')
hold on;
pause
plot(pulsenormnormPSTH(:,(vc==2)), 'r')

%%
%% oct 28th

 off_otherotherother31 = [49,259,272,391,393,783,1368,1430,1741,1862,2177,2539,2597,2971,3616,3905,4066,4518,4681,4727,4730,4771,4786,5072,5089,5495,5506,5596,5841,5926,5943,6062,6091,6212,6213,6483,6541,6888,6901,6993,7157,7203,7502,7562,7593,7668];
offt6_31_init = [ 4786 5596 6888];%1501 2177];

%off_otherotherotherother31 = [724,1501,1731,1773,1816,2867,3661,4204,4578,5117,5431,5853,6106,6376,6542,6796];
%3616? 7502 7668
%offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];

datarun000_31 = get_interspikeinterval(datarun000_31, off_otherotherother31);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherotherother31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, off_otherotherother31(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_31, off_otherotherother31); %or cellids
x = 1:1:30;
normval = [];
auc = [];
minn = [];
tcnormauc = [];
tcnormnorm = [];
tcnormminn = [];
mx = [];
vr = [];
tcnormvr = [];
tcnormmx  = [];
for i = 1:length(off_otherotherother31) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
 auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 

vr = var(tc);
[mx mxt] = max(tc);
[minn minnt] = min(tc);
mxovmn = mx./minn;
mxtmintdiff = mxt - minnt;

auc = repmat(auc, size(tc, 1), 1);
normval = repmat(normval, size(tc, 1), 1);
vr = repmat(vr, size(tc, 1), 1);
mx = repmat(mx, size(tc, 1), 1);
minn = repmat(minn, size(tc, 1), 1);

tcnormauc = tc./auc;
tcnormnorm = tc./normval;
tcnormvr = tc./vr;
tcnormmx = tc./mx;
tcnormminn = tc./minn;

times = [mxt; minnt; mxtmintdiff;];
amp = [minn(1,:); mx(1,:); mxovmn(1,:);];

temp_tcs = get_time_courses_matrix(datarun000_31, off_otherotherother31);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherother31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherother31)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherother31) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);

wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,off_otherotherother31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; %change depending on length of trial
pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(off_otherotherother31)
 pulsePSTH(:,a) = sum(nhist{a,1})./50; %change depending on num of trials
end
 
for i = 1:length(off_otherotherother31)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
end
normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;

[maxpulse maxpulsetime] = max(pulsePSTH);
maxpulsetime = maxpulsetime*0.1;

 
[C ia ib] = intersect(offt6_31_init, off_otherotherother31);
vc = ones(length(off_otherotherother31),1);
vc(ib) = 2;


sumsumsp = [];
for i = 1: length(sumSpTrTrig)
    sumsumsp(i,:) = sum(sumSpTrTrig{i, 1});
end
  
%% Calculate sta parameter
% marks_params.thresh = 4.0;
% datarun000_31 = get_significant_stixels(datarun000_31, off_otherotherother31, 'thresh_params', marks_params);
% 
% for i = 1:length(offt6_31_init)
%     plot_rf(datarun000_31,offt6_31_init(i));
%     hold on;
%     cell_index = get_cell_indices(datarun000_31, offt6_31_init(i));
%     sig_stix = datarun000_31.stas.marks{cell_index};
%     temp_indices = find(full(sig_stix));
%     [temp_rows, temp_cols] = ind2sub(size(sig_stix), temp_indices);
%     plot(temp_cols, temp_rows, 'ro', 'MarkerSize', 7);
%     hold off;
%      pause;
%      close all
% end

%%
close all;
rstd = [];
meanpix = [];
cellind = get_cell_indices(datarun000_31, off_otherotherother31);
stamat = cell(length(cellind),1);
 
for i = 1:length(cellind)
    B = [];
    B = datarun000_31.stas.rfs{cellind(i), 1}(:)';
    meanpix(i) = mean(B);
    rstd(i) = robust_std(B, [1]);
    stamat{i,1} = zeros(size(datarun000_31.stas.rfs{cellind(i),1},1),size(datarun000_31.stas.rfs{cellind(i),1},2));
    for j = 1:size(datarun000_31.stas.rfs{cellind(i),1},1)
            for k = 1:size(datarun000_31.stas.rfs{cellind(i),1},2)
                if (datarun000_31.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 4.6*rstd(i))
                    stamat{i,1}(j,k) = 1;
                end
            end
    end
 
    
end
    
 
 
 
close all;
v = [];
pm = [];
for i = 1:length(cellind)
    rm = [];
    rn = [];
    rmm = [];
    rnn = [];
    DT = [];
    kr = [];
    points = [];
    [rm,rn] = find(stamat{i,1}); %finds all the significant pixels
            if ~isempty(rm)
                [rmm,rnn] = pix_border(rm,rn); %adding everything by .5 and subtracting by 0.5 - i think it is calculating the 4 coordinate values for each pixel - edge correction
                DT = delaunayTriangulation(rmm,rnn); %prediction is still correct with inputs rm and rn, but edges need correcting
                if ~isempty(DT.ConnectivityList)
                    [kr v(1,i)] = convexHull(DT);
                end
                %plot_rf(datarun000_31,off_otherotherother31(i));
                %spy(stamat{i,1})%'LineSpec', 'or');
                [x,y] = find(stamat{i,1});
                clr = stamat{i,1}(stamat{i,1}~=0);
                %scatter(y,x,20,'MarkerFaceColor',[1 0 0],'LineWidth',0.05)
                %set(gca,'Xdir','reverse');%'Ydir','reverse')
                %plot(DT.Points(:,2),DT.Points(:,1), '.','markersize',3);
                %hold on;
                %plot(DT.Points(kr,2), DT.Points(kr,1), 'Color',[1 0 0]);%[rand(1) rand(1) rand(1)]);
                %hold off;
                %pause;
                points(:,1) = DT.Points(kr,2);
                points(:,2) = DT.Points(kr,1);
                perimeter = 0;
                for j = 1:size(points, 1)-1
                    perimeter = perimeter + norm(points(j, :) - points(j+1, :));
                end
                perimeter = perimeter + norm(points(end, :) - points(1, :)); % Last point to first
                pm(i) = perimeter;
            else
            end
            v(2,i) = length(rm);
end

% radius = [];
% radius = get_rf_fit_radius(datarun000_31, off_otherotherother31);
%%
% 
% stavar = [];
% for i = 1:length(off_otherotherother31)
%         ii = get_cell_indices(datarun000_31, off_otherotherother31(1,i));
%   
%         stavar(i,1) = var(datarun000_31.stas.rfs{ii,1}(datarun000_31.stas.significant_stixels{ii, 1}));
%         if(isnan(stavar(i,1)))
%             stavar(i,1) = 0;
%         end
% end

%%
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, off_otherotherother31, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V spave] = dscellanalysis(NumSpikesCell, StimComb);

minfiringslow = [];
maxfiringslow = [];
minfiringfast = [];
maxfiringfast = [];

for i = 1:length(off_otherotherother31)
minfiringslow(i) =  (spave{2,1}(i,3)+spave{2,1}(i,4))/2;
maxfiringslow(i) =  (spave{2,1}(i,1)+spave{2,1}(i,2))/2;
minfiringfast(i) = (spave{1,1}(i,4)+spave{1,1}(i,5))/2;
maxfiringfast(i) = (spave{1,1}(i,3)+spave{1,1}(i,8))/2;
end
%%
scatter(minfiringslow(1,ismember(off_otherotherother31,offt6_31_init)),maxfiringslow(1,ismember(off_otherotherother31,offt6_31_init)),'r');
hold on;
scatter(minfiringslow(1,~ismember(off_otherotherother31,offt6_31_init)),maxfiringslow(1,~ismember(off_otherotherother31,offt6_31_init)),'b');
xlabel('Min Firing - Slow Speed')
ylabel('Max Firing - Slow Speed')
pause;
hold off;


scatter(minfiringfast(1,ismember(off_otherotherother31,offt6_31_init)),maxfiringfast(1,ismember(off_otherotherother31,offt6_31_init)),'r');
hold on;
scatter(minfiringfast(1,~ismember(off_otherotherother31,offt6_31_init)),maxfiringfast(1,~ismember(off_otherotherother31,offt6_31_init)),'b');
xlabel('Min Firing - Fast Speed')
ylabel('Max Firing - Fast Speed')
pause;
hold off;

xx = maxfiringslow./minfiringslow;
yy = maxfiringfast./minfiringfast;


scatter(xx(1,ismember(off_otherotherother31,offt6_31_init)),yy(1,ismember(off_otherotherother31,offt6_31_init)),'r');
hold on;
scatter(xx(1,~ismember(off_otherotherother31,offt6_31_init)),yy(1,~ismember(off_otherotherother31,offt6_31_init)),'b');
xlabel('Max/Min Firing - Slow Speed')
ylabel('Max/Min Firing - Fast Speed')
pause;
hold off;


xx = 1-(minfiringslow./maxfiringslow);
xx(minfiringslow./maxfiringslow > 1) = 0;
yy = 1-(minfiringfast./maxfiringfast);
yy(minfiringfast./maxfiringfast > 1) = 0;

scatter(xx(1,ismember(off_otherotherother31,offt6_31_init)),yy(1,ismember(off_otherotherother31,offt6_31_init)),'r');
hold on;
scatter(xx(1,~ismember(off_otherotherother31,offt6_31_init)),yy(1,~ismember(off_otherotherother31,offt6_31_init)),'b');
xlabel('1-Min/Max Firing - Slow Speed')
ylabel('1-Min/Max Firing - Fast Speed')
pause;
hold off;



% xxx = 1 - (minFir{1,1}./maxFir{1,1});
% xxx(minFir{1,1}./maxFir{1,1} > 1) = 0;
% yyy = 1 - (minFir{2,1}./maxFir{2,1});
% yyy(minFir{2,1}./maxFir{2,1} > 1) = 0;
% scatter(xxx,yyy);

% minangle = [];
% maxangle = [];
% minAxis = [];
% maxAxis = [];
% minAxFir = [];
% maxAxFir = [];
% minFir = cell(2,1);
% maxFir = cell(2,1);
% [CC,minAxis] = min(spave{2,1}');
% for i = 1:length(off_otherotherother31)
%     minangle(1,i) = theta{2,1}(i,minAxis(1,i));
%     minangle(2,i) = minangle(1,i)+pi;
%     if (minangle(2,i) < 0 || minangle(2,i) > pi + 3*pi/4)
%         minangle(2,i) = minangle(1,i) - pi;
%     end
%     minAxis(2,i) = find(minangle(2,i) == theta{2,1}(1,:));
%     maxangle(1,i) = minangle(1,i)+pi./2;
%     if (maxangle(1,i) < 0 || maxangle(1,i) > pi + 3*pi/4)
%         maxangle(1,i) = minangle(1,i)-pi./2;
%     end
%     maxAxis(1,i) = find(maxangle(1,i) == theta{2,1}(1,:));
%     maxangle(2,i) = maxangle(1,i)+pi;
%         if (maxangle(2,i) < 0 || maxangle(2,i) > pi + 3*pi/4)
%                 maxangle(2,i) = maxangle(1,i)-pi;
%         end
%       maxAxis(2,i) = find(maxangle(2,i) == theta{2,1}(1,:));
%       minAxFir(i) = (spave{2,1}(i,minAxis(1,i))+spave{2,1}(i,minAxis(2,i)))/2;
%       maxAxFir(i) = (spave{2,1}(i,maxAxis(1,i))+spave{2,1}(i,maxAxis(2,i)))/2;
% end
% 
% % scatter(minAxFir,maxAxFir);
% % hold on;
% % xlabel('minimum firing');
% % ylabel('maximum firing')
% % title('S 64 T 256')
% % %plot(0:1:250, 0:1:250);
% % hold off;
% % figure()
% % hist(maxAxFir./minAxFir, 20)
% % xlabel('maximum firing / minimum firing');
% % title('S 64 T 256')
% 
% minFir{2,1} = minAxFir;
% maxFir{2,1} = maxAxFir;
% 
% 
% minangle = [];
% maxangle = [];
% minAxis = [];
% maxAxis = [];
% minAxFir = [];
% maxAxFir = [];
% [CC,minAxis] = min(spave{1,1}')
% for i = 1:length(off_otherotherother31)
%     minangle(1,i) = theta{1,1}(i,minAxis(i));
%     minangle(2,i) = minangle(1,i)+pi;
%     if (minangle(2,i) < 0 || minangle(2,i) > pi + 3*pi/4)
%         minangle(2,i) = minangle(1,i) - pi;
%     end
%     minAxis(2,i) = find(minangle(2,i) == theta{1,1}(1,:));
%     maxangle(1,i) = minangle(1,i)+pi./2;
%     if (maxangle(1,i) < 0 || maxangle(1,i) > pi + 3*pi/4)
%         maxangle(1,i) = minangle(1,i)-pi./2;
%     end
%     maxAxis(1,i) = find(maxangle(1,i) == theta{1,1}(1,:));
%     maxangle(2,i) = maxangle(1,i)+pi;
%         if (maxangle(2,i) < 0 || maxangle(2,i) > pi + 3*pi/4)
%                 maxangle(2,i) = maxangle(1,i)-pi;
%         end
%       maxAxis(2,i) = find(maxangle(2,i) == theta{1,1}(1,:));
%       minAxFir(i) = (spave{1,1}(i,minAxis(1,i))+spave{1,1}(i,minAxis(2,i)))/2;
%       maxAxFir(i) = (spave{1,1}(i,maxAxis(1,i))+spave{1,1}(i,maxAxis(2,i)))/2;
% end
% 
% % scatter(minAxFir,maxAxFir);
% % hold on;
% % xlabel('minimum firing');
% % ylabel('maximum firing')
% % title('S 64 T 32')
% % %plot(0:1:250, 0:1:250);
% % hold off;
% % figure();
% % hist(maxAxFir./minAxFir,20)
% % xlabel('maximum firing / minimum firing');
% % title('S 64 T 32')
% 
% minFir{1,1} = minAxFir;
% maxFir{1,1} = maxAxFir;

% xxx = 1 - (minFir{1,1}./maxFir{1,1});
% xxx(minFir{1,1}./maxFir{1,1} > 1) = 0;
% yyy = 1 - (minFir{2,1}./maxFir{2,1});
% yyy(minFir{2,1}./maxFir{2,1} > 1) = 0;
% scatter(xxx,yyy);
% xlabel('maximum firing / minimum firing - T 32');
% ylabel('maximum firing / minimum firing - T 256')
% title('S 64 T 32-256')
% figure()
% hist(xxx.*yyy,100);
% xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');

%%
      num_rgcs = length(off_otherotherother31);

% initialize some variables for the look
rf_areas = zeros(num_rgcs,1);
abs_mean_pixel_val = zeros(num_rgcs,1);
snrs = zeros(num_rgcs, 1);
contrast_index = zeros(num_rgcs, 1);
[Y, X] = meshgrid(1:1:40, 1:1:80);


for rgc = 1:num_rgcs
    
    temp_index= get_cell_indices(datarun000_31, off_otherotherother31(rgc));
    
    %plot_rf(datarun000_31, datarun000_31.cell_ids(rgc), 'sig_stix', true)
    
    [I, J] = find(full(datarun000_31.stas.marks{temp_index}));
    if length(I) < 4
        rf_areas(rgc) = 1;
        warning('report')
    else
        % get the convex hull defined by the significant stixels
        hull_coords = convhull(I,J);
        
        % get the area of the convex hull and save it for every rgc
        rf_areas(rgc) = polyarea(J(hull_coords),I(hull_coords));
        
        % get pixels that are in polygon defined by convex hull
        IN = inpolygon(X, Y, J(hull_coords), I(hull_coords));
        
        % get the number of pixels for averaging
        temp_num_pixels = length(find(IN > 0));

        % get the RF for this rgc
        temp_rf = get_rf(datarun000_31, off_otherotherother31(rgc));
 
      
        % normalize the rf
        temp_rf = temp_rf ./ abs(ext(reshape(temp_rf, [], 1)));
        
        % mask the RF by the pixels in the polygon 
        masked_rf = temp_rf .* IN';
        
        % number of positive pixels
        num_positive = length(find(masked_rf > 0));
        num_negative = length(find(masked_rf < 0));
        
        % compute contrast index
        contrast_index(rgc) = 1 - abs(((num_positive - num_negative) ./ temp_num_pixels));
    end
end
%%

% pulse isi 1-ratio CI 
% (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));SCORE1(:,1); 1- (v(2,:)./v(1,:)); contrast_index;

[COEFF1,SCORE1] = princomp(isinormnorm');

X = [];

X(:,1) = SCORE1(:,1);
X(:,2) = 1- (v(2,:)./v(1,:));
X(:,3) = contrast_index;
scatter3(X(ismember(off_otherotherother31,offt6_31_init),1), X(ismember(off_otherotherother31,offt6_31_init),2), X(ismember(off_otherotherother31,offt6_31_init),3),'r');
hold on;
scatter3(X(~ismember(off_otherotherother31,offt6_31_init),1), X(~ismember(off_otherotherother31,offt6_31_init),2),  X(~ismember(off_otherotherother31,offt6_31_init),3),'b');


%%



[COEFF1,SCORE1] = princomp(isinormnorm');

X = [];

X(:,1) = contrast_index;%(sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = 1- (v(2,:)./v(1,:));%SCORE1(:,1);
scatter(X(ismember(off_otherotherother31,offt6_31_init),1), X(ismember(off_otherotherother31,offt6_31_init),2), 'r');
hold on;
scatter(X(~ismember(off_otherotherother31,offt6_31_init),1), X(~ismember(off_otherotherother31,offt6_31_init),2), 'b');
xlabel('CI');
ylabel('1-ratio');
% zlabel('CI');
hold off;
pause;


X = [];
X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
scatter(X(ismember(off_otherotherother31,offt6_31_init),1), X(ismember(off_otherotherother31,offt6_31_init),2), 'r');
hold on;
scatter(X(~ismember(off_otherotherother31,offt6_31_init),1), X(~ismember(off_otherotherother31,offt6_31_init),2), 'b');
xlabel('Pulse');
ylabel('ISI');
hold off;
pause;



%% ISI, TC PARAMS, TC PC, PULSE PC, PULSE PARAMS, STA
%isi 1 2 3  pulsesum pulsenorm 1 pulsenorm 2 pulsenorm 3 pulse 1 pulse 2 pulse 3  stavar
%X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));

[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsePSTH');
[COEFF2,SCORE2] = princomp(tcnormminn');

X = [];

X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
X(:,3) = (1-(v(2,:)./v(1,:)));%contrast_index';
scatter3 (X(ismember(off_otherotherother31,offt6_31_init),1), X(ismember(off_otherotherother31,offt6_31_init),2), X(ismember(off_otherotherother31,offt6_31_init),3), 'r');
hold on;
scatter3 (X(~ismember(off_otherotherother31,offt6_31_init),1), X(~ismember(off_otherotherother31,offt6_31_init),2), X(~ismember(off_otherotherother31,offt6_31_init),3), 'b');
xlabel('pulse');
ylabel('isi');
zlabel('CI');

%%
[idx] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_otherotherother31, tcnormnorm,0,vc);
%%
vcc  =[];
vcc(1,:) = mean(X(ib,:));
vcc(2,:) = mean(X(vc==1,:));
[idx] = clustering_analysis_plots(X, 1,0, 2, 1, 0, datarun000_31, off_otherotherother31, tcnormnorm,0,vcc);

%%
scatter3 (X(:,1), X(:,2), X(:,3))
 %offt4_31 = off_otherotherother31(idx==2);
 %off_otherotherother31 = off_otherotherother31(idx==1);
%  off_otherotherother31(idx==2)
% ismember(offt6_31_init, off_otherotherother31(idx==2))
plot(tcnormnorm(:,(vc==1)), 'b')
hold on;
pause
plot(tcnormnorm(:,(vc==2)), 'r')
hold off;
pause;
plot(isinormnorm(:,(vc==1)), 'b')
hold on;
pause
plot(isinormnorm(:,(vc==2)), 'r')
hold off;
pause
plot(pulsenormnormPSTH(:,(vc==1)), 'b')
hold on;
pause
plot(pulsenormnormPSTH(:,(vc==2)), 'r')

%% on t3 - 10-10

plot_rf_summaries(datarun000_10, [ont3_10,6425,7384, 468, 842, 1129,1278,1384,1414,1925, 1967,2374,2687,2943,3125,3721,3887,4819,5148,5493,6137], 'coordinates', 'monitor');
plot_time_courses(datarun000_10,  [ont3_10,6425,7384, 468, 842, 1129,1278,1384,1414,1925, 1967,2374,2687,2943,3125,3721,3887,4819,5148,5493,6137],'all', true, 'bw', true);
ismember([6425,7384, 468, 842, 1129,1278,1384,1414,1925, 1967,2374,2687,2943,3125,3721,3887,4819,5148,5493,6137], cellids_10)

%% ont2; 10-10

plot_rf_summaries(datarun000_10, [ont2_10 241 1291 1802 4681 5222 5237 5867], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_10,  [ont2_10 241 1291 1802 4681 5222 5237 5867],'all', true, 'bw', true);

%% ont1 10-10

plot_rf_summaries(datarun000_10, [ont1_10 4216 7501 7576 6571 2626 ], 'coordinates', 'monitor', 'label', true);
figure();
plot_time_courses(datarun000_10,  [ont1_10 4216 7501 7576 6571 2626 ],'all', true, 'bw', true);

%% nov 26th finding extra cells for off t4, t3, t5, t6
temp_tcs = get_time_courses_matrix(datarun000_10, datarun000_10.cell_ids);
tc_fit = [];
final_params  =[];
for i = 1:length(datarun000_10.cell_ids)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(datarun000_10.cell_ids)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(datarun000_10.cell_ids) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all fitted time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[minn minnt] = min(tcfittednormnorm);
[eval extt] = max(abs(tcfittednormnorm));
extrval = [];
for i = 1:length(extt)
extrval(i) = tcfittednormnorm(extt(i),i);
end


[tc nontc] = get_time_courses_matrix(datarun000_10, datarun000_10.cell_ids); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(datarun000_10.cell_ids) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

on_init10 = [31,77,211,242,244,272,558,631,646,677,783,813,1036,1081,1100,1144,1189,1281,1397,1442,1517,1531,1594,1696,1697,1759,1863,1939,2118,2134,2177,2401,2402,2435,2464,2491,2703,2716,2748,3002,3031,3213,3287,3316,3331,3422,3512,3515,3587,3618,3811,4024,4055,4126,4127,4174,4262,4443,4488,4534,4548,4711,4712,4877,5087,5146,5176,5496,5506,5642,5671,5702,5777,5821,5868,5957,6136,6197,6272,6275,6319,6572,6722,6751,6752,6767,6796,6977,7096,7189,7232,7354,7426,7472,7564,7610,7621,7668];
[C ia ib] = intersect(on_init10, datarun000_10.cell_ids);
vc = ones(length(datarun000_10.cell_ids),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_10, datarun000_10.cell_ids, tcnormnorm,0, vc);
on_allsnr10_all = datarun000_10.cell_ids(idx ==2);
off_allsnr10_all = datarun000_10.cell_ids(idx ==1);

%%

c = get_cell_indices(datarun000_10, on_allsnr10_all);
snronall = [];
for i = 1:length(c)
    r1 = sort(datarun000_10.stas.rfs{c(1,i),1}(:)', 'descend');
    snronall(1,i) = mean(r1(1:4))./std(r1);
end
on_10_all = on_allsnr10_all(snronall > (mean(snronall) - 2.5*std(snronall)));

% hax=axes; 
% hold on;
% hist(snronall)
% SP= mean(snronall) - 2.5*std(snronall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('On cutoff')

c = get_cell_indices(datarun000_10, off_allsnr10_all);
snroffall = [];
for i = 1:length(c)
    r1 = sort(datarun000_10.stas.rfs{c(1,i),1}(:)', 'descend');
    snroffall(1,i) = mean(r1(1:4))./std(r1);
end
off_10_all = off_allsnr10_all(snroffall > (mean(snroffall) - 2.5*std(snroffall)));
 %% end of nov 26th finding extra cells for off t4, t3, t5, t6

offt1_10_all = [offt1_10 2162 ];
offt2_10_all =  [offt2_10 4728];
offt4_10_all = [offt4_10  3 3901 4115 5117 7667];
offt3_10_all =  [offt3_10 1579 2012 3199 3710 4160 6079 6151 7324 1580];
offt5_10_all = [offt5_10 574 1533 1698 2686 2971 4878 5446 6573  2495];
off_otherotherotherotherother10_all = off_10_all(~ismember(off_10_all, [offt1_10_all offt2_10_all offt4_10_all offt3_10_all offt5_10_all]));
offt6_10_all = [5836 6856 1486 496];

%% Drifting Grating PSTH

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_10, ont1_10, 0);

ind = find(StimComb(:,2)==256);
zrind = find(StimComb(:,3)==0);
zeroind = intersect(ind, zrind);
psthall = [];
for i = 1:length(ont1_10)
    psthzero = [];
    psth = [];
    [T, psthzero, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont1_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,zeroind)), 'stop', 8, 'bin_size', 0.1);
    psthzerocut = psthzero(21:61);
    for k = 1:length(ind)
        psth01 = [];
        [T, psth01, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont1_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,ind(k))), 'stop', 8, 'bin_size', 0.1);
        y = [];yy = [];
        y = psth01(21:61);
        dp = [];
        for j = 1:length(y)
            yy = circshift(y,[0,j]);
            dp(j) = dot(psthzerocut,yy);
            yy = [];
        end
        ymx = []; imx = [];
        [ymx,imx]=max(dp);
        yy = circshift(y,[0,imx]);
        psth(k,:) = yy;
    end
    
    x = [];xmx=[];imx = []; xx = [];
    x = mean(psth);
    [xmx,imx]=max(x);
    xx = circshift(x,[0,length(x)-imx+1]);
    
    psthall(i,:) = xx;
    %plot(2:0.1:6,xx);
% hold on;
% plot(2:0.1:6,yy, 'r');
% title('Dot');
% hold off; 


%pause;

end

plot(mean(psthall));
pause;
shadedErrorBar(2:0.1:6,mean(psthall),std(psthall),'k');

%%
dgresp_10(1,:) = mean(psthall);
dgresp_10(2,:) = mean(psthall);
dgresp_10(3,:) = mean(psthall);

plot(dgresp_10','DisplayName','dgresp')
legend('t1', 't2', 't3');
title('2012-10-10');

%%

spamp = [];
spamp(1,:) = mean(dgresp_10');
spamp(2,:) = mean(dgresp_15');
spamp(3,:) = mean(dgresp_31');


y(2,1) = mean(A256(ismember(on_10,ont1_10)));
y(2,2) = mean(A256(~ismember(on_10,ont1_10)));


b = bar(spamp, 0.5);


set(b(1), 'LineWidth', 2, 'EdgeColor' , [1 0 0], 'FaceColor' , [1 0.7 0.8])
set(b(2), 'LineWidth', 2, 'EdgeColor' , [0 0 1], 'FaceColor' , [.7 0.8 1])
set(b(3), 'LineWidth', 2, 'EdgeColor' , [0 1 0], 'FaceColor' , [.7 1 0.8])
set(gca,'XTickLabel',{'10-10', '10-15', '10-31'})
legend('t1', 't2', 't3');
title('Average Spike Amplitude')

%% Plot Drifitng grating PSTH and Pulse PSTH Side by Side
%Check if On-Off cells have F2 and On cells have F1 - for DS cell
%classification into On and On-Off
% Normalized Pulse PSTH plotted
% DG: Plot PSTH (averaged over 8 trials) for preferred and null directions.
% Null direction is shifted wrt preferred using dot product overlap so both

close all;

wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,ds_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
b = 1; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5;

ind = find(StimComb(:,2)==256);
zrind = find(StimComb(:,3)==0);
zeroind = intersect(ind, zrind);
psthall = [];

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, ds_15, 0);

ind = find(StimComb(:,2)==256);
zrind = find(StimComb(:,3)==0);
zeroind = intersect(ind, zrind);
psthall = [];

for i = 1:length(ds_15)
    [maxf prf] = max(NumSpikesCell(i, ind));
    [minf nll] = min(NumSpikesCell(i, ind));
    psthzero = [];
    psth = [];
    [T, psthzero, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ds_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,ind(prf))), 'stop', 8, 'bin_size', 0.25);
    psthzerocut = psthzero(9:25);
    psth01 = [];
    [T, psth01, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ds_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,ind(nll))), 'stop', 8, 'bin_size', 0.25);
    y = [];yy = [];
    y = psth01(9:25);
    dp = [];
    for j = 1:length(y)
        yy = circshift(y,[0,j]);
        dp(j) = dot(psthzerocut,yy);
        yy = [];
    end
    ymx = []; imx = [];
    [ymx,imx]=max(dp);
    yy = circshift(y,[0,imx]);
    psth02 = yy;
    subplot(1,2,1);
    plot(2:0.25:6, psth02, 'b');
    hold on;
    plot(2:0.25:6, psthzerocut, 'r');
    legend('null', 'preferred');
    hold off;
    
     psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
    subplot(1,2,2);
    stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
    hold on;
    plot(binSize,psthind, 'b');
    plot(binSize,psthindnorm, 'r');
    set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
    set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
    set(gcf, 'Color', [1 1 1]);
    hold off;
    pause;
end

%% Space Time Receptive Field - offt4
% Plot x vs t or y vs t
%X vs t - collapse all y values onto x points - take average
%y vs t - collapse all x values onto y - take average
%Is average the right measure to collapse all pixels onto 1 d

for i = 1: length(ds_10)
    stayt = [];
    staxt = [];
    for j = 1:40
        stayt(j,:) = mean(datarun000_10.stas.stas{get_cell_indices(datarun000_10, ds_10(i)),1}(j,:,:,:));
    end
    for k = 1:80
        staxt(k,:) = mean(datarun000_10.stas.stas{get_cell_indices(datarun000_10, ds_10(i)),1}(:,k,:,:));
    end
    imagesc(norm_image(staxt));
    pause;
end

    




