%% load dataruns - drifting grating, white noise, pulses

addpath('/Users/sneharavi/Documents/MATLAB/Classification/');
addpath('/Users/sneharavi/Documents/MATLAB/DS cell analysis/');

[datarun000_15] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-15-0/data000-3600-7200s/datamaps-sr-model/', 'data000-map/data000-map', 0, 0, 1);
[datarun001_15] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-15-0/data000-3600-7200s/datamaps-sr-model/', 'data001-map/data001-map', 0, 0, 0);
[datarun002_15] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-15-0/data000-3600-7200s/datamaps-sr-model/', 'data002-map/data002-map', 1, '/stimuli/s02', 0);

%% DS

cellids_15 = intersect((intersect(datarun000_15.cell_ids, datarun001_15.cell_ids)), datarun002_15.cell_ids);
%cellids_15 = [4,31,46,62,76,94,154,182,226,257,272,301,333,347,407,424,438,454,467,496,514,528,586,649,692,708,751,753,766,768,782,783,843,857,860,872,888,889,901,979,991,1051,1067,1081,1084,1097,1098,1111,1127,1130,1156,1191,1246,1310,1339,1381,1384,1385,1398,1411,1427,1486,1532,1549,1564,1578,1581,1591,1595,1637,1652,1683,1684,1685,1786,1817,1876,1877,1892,1893,1895,1908,1921,1966,1969,1999,2011,2026,2042,2136,2161,2177,2192,2206,2208,2236,2253,2311,2326,2328,2343,2357,2373,2401,2417,2449,2461,2462,2478,2521,2522,2555,2597,2686,2716,2732,2746,2794,2809,2851,2868,2881,2896,2898,2929,3002,3061,3091,3139,3152,3198,3200,3226,3241,3244,3258,3286,3287,3303,3319,3422,3452,3512,3530,3559,3586,3589,3601,3634,3635,3636,3648,3679,3691,3692,3695,3721,3736,3767,3812,3813,3815,3842,3857,3859,3889,3917,3933,3934,3935,3946,3991,3994,4006,4021,4022,4069,4096,4097,4098,4130,4145,4157,4173,4188,4231,4234,4235,4246,4278,4279,4294,4324,4353,4427,4442,4459,4486,4487,4501,4503,4562,4591,4668,4697,4731,4732,4771,4774,4788,4789,4846,4864,4892,4941,4985,4997,4998,4999,5071,5073,5088,5104,5116,5148,5150,5179,5223,5359,5405,5433,5446,5464,5567,5569,5632,5641,5645,5657,5658,5672,5702,5703,5705,5719,5733,5791,5836,5851,5853,5896,5898,5927,5941,6034,6064,6093,6106,6125,6139,6140,6152,6155,6170,6196,6229,6257,6260,6286,6304,6321,6332,6361,6363,6376,6380,6391,6392,6422,6439,6451,6455,6512,6542,6589,6721,6722,6737,6751,6752,6797,6811,6826,6828,6886,6903,6931,6976,6980,6992,7021,7040,7067,7069,7096,7157,7186,7203,7234,7261,7278,7306,7308,7354,7442,7471,7475,7487,7503,7517,7520,7532,7562,7667];
[tc nontc] = get_time_courses_matrix(datarun000_15, cellids_15); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(cellids_15) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

%DS cells
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, cellids_15, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
ds_init15 = [257,301,333,438,467,708,1097,1595,1683,1685,1895,2042,2449,2898,3512,3636,3736,3815,3842,3934,4069,4157,4173,4353,4731,4789,4846,4985,5150,5632,5702,5719,6229,6321,6332,6751,6797,7308, 766,5148, 7475]

[C ia ib] = intersect(ds_init15, cellids_15);
vc = ones(length(cellids_15),1);
vc(ib) = 2;

close all;
X = [];
N = [];
p = [];
X(:,1) = log(mag{1,1})';
X(:,2) = log(mag{2,1})';
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_15, cellids_15, tcnormnorm,0, vc);

ds_15 = [];
ds_15 = cellids_15(idx==2);
nonds_15 = cellids_15(idx==1);
%lr = sum(p(idx==2,2))/length(ds_15)
%lr = sum(p(idx==1,1))/length(nonds_15)
%ds_15 = [257,301,333,438,467,708,1097,1595,1683,1685,1895,2042,2449,2898,3512,3636,3736,3815,3842,3934,4069,4157,4173,4353,4731,4789,4846,4985,5150,5632,5702,5719,6229,6321,6332,6751,6797,7308];
%nonds_15 = [4,31,46,62,76,94,154,182,226,272,347,407,424,454,496,514,528,586,649,692,751,753,766,768,782,783,843,857,860,872,888,889,901,979,991,1051,1067,1081,1084,1098,1111,1127,1130,1156,1191,1246,1310,1339,1381,1384,1385,1398,1411,1427,1486,1532,1549,1564,1578,1581,1591,1637,1652,1684,1786,1817,1876,1877,1892,1893,1908,1921,1966,1969,1999,2011,2026,2136,2161,2177,2192,2206,2208,2236,2253,2311,2326,2328,2343,2357,2373,2401,2417,2461,2462,2478,2521,2522,2555,2597,2686,2716,2732,2746,2794,2809,2851,2868,2881,2896,2929,3002,3061,3091,3139,3152,3198,3200,3226,3241,3244,3258,3286,3287,3303,3319,3422,3452,3530,3559,3586,3589,3601,3634,3635,3648,3679,3691,3692,3695,3721,3767,3812,3813,3857,3859,3889,3917,3933,3935,3946,3991,3994,4006,4021,4022,4096,4097,4098,4130,4145,4188,4231,4234,4235,4246,4278,4279,4294,4324,4427,4442,4459,4486,4487,4501,4503,4562,4591,4668,4697,4732,4771,4774,4788,4864,4892,4941,4997,4998,4999,5071,5073,5088,5104,5116,5148,5179,5223,5359,5405,5433,5446,5464,5567,5569,5641,5645,5657,5658,5672,5703,5705,5733,5791,5836,5851,5853,5896,5898,5927,5941,6034,6064,6093,6106,6125,6139,6140,6152,6155,6170,6196,6257,6260,6286,6304,6361,6363,6376,6380,6391,6392,6422,6439,6451,6455,6512,6542,6589,6721,6722,6737,6752,6811,6826,6828,6886,6903,6931,6976,6980,6992,7021,7040,7067,7069,7096,7157,7186,7203,7234,7261,7278,7306,7354,7442,7471,7475,7487,7503,7517,7520,7532,7562,7667];

%% ON - OFF Cells
temp_tcs = get_time_courses_matrix(datarun000_15, nonds_15);
tc_fit = [];
final_params  =[];
for i = 1:length(nonds_15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(nonds_15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(nonds_15) %or nonds
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


[tc nontc] = get_time_courses_matrix(datarun000_15, nonds_15); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(nonds_15) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

on_init15 = [4,31,154,182,272,424,454,496,528,692,766,783,857,860,888,979,1067,1111,1339,1385,1398,1486,1532,1549,1578,1591,1786,1876,1877,1892,1969,2011,2026,2136,2161,2208,2311,2328,2417,2461,2597,2716,2851,2929,3002,3152,3198,3241,3244,3303,3319,3422,3452,3530,3559,3586,3634,3691,3721,3812,3859,3917,3946,3991,3994,4022,4096,4130,4145,4231,4246,4279,4427,4501,4562,4668,4774,4892,4941,4997,4998,5073,5088,5104,5179,5359,5405,5567,5657,5791,5853,5898,5927,5941,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6542,6722,6737,6752,6826,6903,6980,7067,7096,7157,7186,7203,7261,7354,7442,7487,7520,7532,7562];
[C ia ib] = intersect(on_init15, nonds_15);
vc = ones(length(nonds_15),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_15, nonds_15, tcnormnorm,0, vc);
on_allsnr15 = nonds_15(idx ==2)
off_allsnr15 = nonds_15(idx ==1)
%on_allsnr15 = [4,31,154,182,272,424,454,496,528,692,766,783,857,860,888,979,1067,1111,1339,1385,1398,1486,1532,1549,1578,1591,1786,1876,1877,1892,1969,2011,2026,2136,2161,2208,2311,2328,2417,2461,2597,2716,2851,2929,3002,3152,3198,3241,3244,3303,3319,3422,3452,3530,3559,3586,3634,3691,3721,3812,3859,3917,3946,3991,3994,4022,4096,4130,4145,4231,4246,4279,4427,4501,4562,4668,4774,4892,4941,4997,4998,5073,5088,5104,5179,5359,5405,5567,5657,5791,5853,5898,5927,5941,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6542,6722,6737,6752,6826,6903,6980,7067,7096,7157,7186,7203,7261,7354,7442,7487,7520,7532,7562];
%off_allsnr15 = [46,62,76,94,226,347,407,514,586,649,751,753,768,782,843,872,889,901,991,1051,1081,1084,1098,1127,1130,1156,1191,1246,1310,1381,1384,1411,1427,1564,1581,1637,1652,1684,1817,1893,1908,1921,1966,1999,2177,2192,2206,2236,2253,2326,2343,2357,2373,2401,2462,2478,2521,2522,2555,2686,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3287,3589,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3933,3935,4006,4021,4097,4098,4188,4234,4235,4278,4294,4324,4442,4459,4486,4487,4503,4591,4697,4732,4771,4788,4864,4999,5071,5116,5148,5223,5433,5446,5464,5569,5641,5645,5658,5672,5703,5705,5733,5836,5851,5896,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6391,6422,6451,6589,6721,6811,6828,6886,6931,6976,6992,7021,7040,7069,7234,7278,7306,7471,7475,7503,7517,7667];
%% SNR CUTOFF
c = get_cell_indices(datarun000_15, on_allsnr15);
snronall = [];
for i = 1:length(c)
    r1 = sort(datarun000_15.stas.rfs{c(1,i),1}(:)', 'descend');
    snronall(1,i) = mean(r1(1:4))./std(r1);
end
on_15 = on_allsnr15(snronall > (mean(snronall) - 2.5*std(snronall)));

% hax=axes; 
% hold on;
% hist(snronall)
% SP= mean(snronall) - 2.5*std(snronall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('On cutoff')

c = get_cell_indices(datarun000_15, off_allsnr15);
snroffall = [];
for i = 1:length(c)
    r1 = sort(datarun000_15.stas.rfs{c(1,i),1}(:)', 'descend');
    snroffall(1,i) = mean(r1(1:4))./std(r1);
end
off_15 = off_allsnr15(snroffall > (mean(snroffall) - 2.5*std(snroffall)));

% hax=axes; 
% hold on;
% hist(snroffall)
% SP= mean(snroffall) - 2.5*std(snroffall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('Off cutoff')

snroff_15 = off_allsnr15(snroffall < (mean(snroffall) - 2.5*std(snroffall)));
snron_15 = on_allsnr15(snronall < (mean(snronall) - 2.5*std(snronall)));

% snron_15 = [424 4279];
% snroff_15 = [7475];
%on_15 = [4,31,154,182,272,454,496,528,692,766,783,857,860,888,979,1067,1111,1339,1385,1398,1486,1532,1549,1578,1591,1786,1876,1877,1892,1969,2011,2026,2136,2161,2208,2311,2328,2417,2461,2597,2716,2851,2929,3002,3152,3198,3241,3244,3303,3319,3422,3452,3530,3559,3586,3634,3691,3721,3812,3859,3917,3946,3991,3994,4022,4096,4130,4145,4231,4246,4427,4501,4562,4668,4774,4892,4941,4997,4998,5073,5088,5104,5179,5359,5405,5567,5657,5791,5853,5898,5927,5941,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6542,6722,6737,6752,6826,6903,6980,7067,7096,7157,7186,7203,7261,7354,7442,7487,7520,7532,7562];
% off_15 = [46,62,76,94,226,347,407,514,586,649,751,753,768,782,843,872,889,901,991,1051,1081,1084,1098,1127,1130,1156,1191,1246,1310,1381,1384,1411,1427,1564,1581,1637,1652,1684,1817,1893,1908,1921,1966,1999,2177,2192,2206,2236,2253,2326,2343,2357,2373,2401,2462,2478,2521,2522,2555,2686,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3287,3589,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3933,3935,4006,4021,4097,4098,4188,4234,4235,4278,4294,4324,4442,4459,4486,4487,4503,4591,4697,4732,4771,4788,4864,4999,5071,5116,5148,5223,5433,5446,5464,5569,5641,5645,5658,5672,5703,5705,5733,5836,5851,5896,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6391,6422,6451,6589,6721,6811,6828,6886,6931,6976,6992,7021,7040,7069,7234,7278,7306,7471,7503,7517,7667];
%% ON - OFF after SNR check

nondssnr_15 = [on_15 off_15];

temp_tcs = get_time_courses_matrix(datarun000_15, nondssnr_15);
tc_fit = [];
final_params  =[];
for i = 1:length(nondssnr_15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(nondssnr_15) %fit time course
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(nondssnr_15) %or nonds
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


[tc nontc] = get_time_courses_matrix(datarun000_15, nondssnr_15); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(nondssnr_15) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

on_initsnr15 = on_15;
[C ia ib] = intersect(on_initsnr15, nondssnr_15);
vc = ones(length(nondssnr_15),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_15, nondssnr_15, tcnormnorm,0, vc);
onon_15 = nondssnr_15(idx ==2); %idx might change so be careful - on might be idx 2 and off idx 1
offoff_snr15 = nondssnr_15(idx ==1);


    plot(X(idx==2,1),X(idx==2,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
    hold on;
        plot(X(idx==1,1),X(idx==1,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]); 
            legend('Cluster 1 - On Cells','Cluster 2 - Off Cells', 'Location','NW')

set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON-OFF/', 'classfaftersnr', gcf)



%% ON T1 ------
%2012-10-15-0: 123on, 38 cells 1 t1 missing, 1 added
on_15 = [4,31,154,182,272,454,496,528,692,766,783,857,860,888,979,1067,1111,1339,1385,1398,1486,1532,1549,1578,1591,1786,1876,1877,1892,1969,2011,2026,2136,2161,2208,2311,2328,2417,2461,2597,2716,2851,2929,3002,3152,3198,3241,3244,3303,3319,3422,3452,3530,3559,3586,3634,3691,3721,3812,3859,3917,3946,3991,3994,4022,4096,4130,4145,4231,4246,4427,4501,4562,4668,4774,4892,4941,4997,4998,5073,5088,5104,5179,5359,5405,5567,5657,5791,5853,5898,5927,5941,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6542,6722,6737,6752,6826,6903,6980,7067,7096,7157,7186,7203,7261,7354,7442,7487,7520,7532,7562];
ont1_init15 = [4,154,692,860,1111,1339,1786,1877,1892,2011,2161,2208,2461,2851,3002,3319,3422,3691,3994,4774,5073,5088,5359,5567,5853,6542,6826,7261,7442,7487,7532];% 1969 3721 4096 4427 4501 5104 5941];
% clustering adds 7 more cells: 1969 3721 4096 4427 4501 5104 5941 - don't need 1969, but it is also fine - tc is only a little different -this was previously - new clustering only adds 5 more cell
% unsure t1s 3152 is also okay - but adding it to IC makes clusters add many more or less cells
%plot_time_courses(datarun000_15,[4,154,692,860,1111,1339,1786,1877,1892,2011,2161,2208,2461,2851,3002,3319,3422,3691,3994,4774,5073,5088,5359,5567,5853,6542,6826,7261,7442,7487,7532 3721 4096 4427 4501 5104 5941 1969], 'all', true, 'bw', true);
% dg

% 2nd plot uses min value and maxvalue - DOT vs ZC
% pulses?
%Pm vs 1 - ratio 

temp_tcs = get_time_courses_matrix(datarun000_15, on_15);
tc_fit = [];
final_params  =[];
for i = 1:length(on_15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_15) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)


[tc nontc] = get_time_courses_matrix(datarun000_15, on_15); %or cellids
x = 1:1:30;
auc = [];
mx = [];
normval = [];
tcnormnorm = [];
tcnormauc = [];
tcnormmx = [];
for i = 1:length(on_15) %or nonds
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

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, on_15, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 64),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
close all;

vc = [];
[C ia ib] = intersect(ont1_init15, on_15);
vc = ones(length(on_15),1);
vc(ib) = 2;

[COEFF,SCORE] = princomp(tcnormmx');


X = [];
X(:,1) = A32';
X(:,2) = A256';
%X(:,3) = SCORE(:,3);
X(:,3) = TCParams.dot';
%X(:,1) = TCParams.mintim';
%X(:,2) = TCParams.maxtim';
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, on_15, tcnormnorm,0, vc);
on_15(idx==2)
%sum(log(p(idx==2, 2)./p(idx==2, 1)))
ont1_15 = on_15(idx==2);
on_other15 = on_15(idx==1);
ont1_15 = [4,154,692,860,1111,1339,1786,1877,1892,2011,2161,2208,2461,2851,3002,3319,3422,3691,3721,3994,4096,4427,4774,5073,5088,5104,5359,5567,5853,5941,6542,6826,7261,7442,7487,7532];
on_other15 = [31,182,272,454,496,528,766,783,857,888,979,1067,1385,1398,1486,1532,1549,1578,1591,1876,1969,2026,2136,2311,2328,2417,2597,2716,2929,3152,3198,3241,3244,3303,3452,3530,3559,3586,3634,3812,3859,3917,3946,3991,4022,4130,4145,4231,4246,4501,4562,4668,4892,4941,4997,4998,5179,5405,5657,5791,5898,5927,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6722,6737,6752,6903,6980,7067,7096,7157,7186,7203,7354,7520,7562];
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1classf2', gcf)

%%
plot_rf_summaries(datarun000_15, ont1_15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1rf', gcf)

%%
plot_time_courses(datarun000_15,ont1_15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1tc', gcf)
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1tccomp', gcf)

%% Plots ISIs of all cells - ave and std deviation
datarun000_15 = get_interspikeinterval(datarun000_15, ont1_15);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont1_15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, ont1_15(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1isifull', gcf)

%% Plot Pulse response PSTH - ave and std deviation
wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,ont1_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont1_15)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end
b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

spikesbytrials{1,1} = get_raster(datarun001_15.spikes{get_cell_indices(datarun001_15, ont1_15(1,1)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1psthcell4', gcf);
close;
figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
plot(binSize,psthnorm(1,:), 'b');
stairs([0 3 5 8 10],[0.4 0.35 0.3 0.35 0.35], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1psth2', gcf)

%% Plot Drifting Grating Response

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, ont1_15, 0);
NS2 = [];
A32 = [];
A256 = [];
B = [];
B2  =[];
IX = [];
IX2 = [];
NS2 = NumSpikesCell';
A32 = NS2(find(StimComb(:,2) == 64),:); 
A32 = A32./8;
A256 = NS2(find(StimComb(:,2) == 256),:);
A256 = A256./8;
[B, IX] = sort(StimComb(find(StimComb(:,2) == 64),3));
A32 = A32(IX,:);
[B2, IX2] = sort(StimComb(find(StimComb(:,2) == 256),3));
A256 = A256(IX2,:);

A256norm = [];
A32norm = [];
for i = 1:length(A256)
    A256norm(:,i) = A256(:,i)./norm(A256(:,i));
    A32norm(:,i) = A32(:,i)./norm(A32(:,i));
end

figure();
shadedErrorBar(0:45:315,mean(A32'), std(A32'), 'k');
hold on;
shadedErrorBar(0:45:315,mean(A256'), std(A256'), 'k');
xlabel('angle');
ylabel('Average Spikes per s');
ylim([0 50]);

figure();
shadedErrorBar(0:45:315,mean(A32norm'), std(A32norm'), 'k');
hold on;
shadedErrorBar(0:45:315,mean(A256norm'), std(A256norm'), 'k');
xlabel('angle');
ylabel('Average Spikes per s');
ylim([0.1 0.5]);


%%
datarun000_15 = get_autocorrelations(datarun000_15, ont1_15);

%%
x2 = 0:0.0005:0.1; 
%nonds - cells not ds, tc - their time courses
acf = [];
normvalacf = [];
acfnormnorm = [];
maxacf = [];
acfmax = [];
for i = 1:length(ont1_15) %or nonds
 acf(:,i) = datarun000_15.autocorrelation{get_cell_indices(datarun000_15, ont1_15(1,i)), 1}.probabilities;
 normvalacf(1, i) = norm( acf(:,i));
 maxacf(1,i) = max(acf(:,i));
end 
normvalacf = repmat(normvalacf, size(acf, 1), 1);
maxacf = repmat(maxacf, size(acf, 1), 1);
acfnormnorm = acf./normvalacf;
acfmax = acf./maxacf;


figure();
shadedErrorBar(x2(1:100),mean(acfnormnorm(1:100, :)'),std(acfnormnorm(1:100, :)'),'k');
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1acfhalf', gcf)
figure();
shadedErrorBar(x2,mean(acfnormnorm(:, :)'),std(acfnormnorm(:, :)'),'k');
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T1/', 'ont1acffull', gcf)


%% ON T2
%dg, rf, time course params dot area zc time value, tc pc
%later: isi, pulses
on_other15 = [31,182,272,454,496,528,766,783,857,888,979,1067,1385,1398,1486,1532,1549,1578,1591,1876,1969,2026,2136,2311,2328,2417,2597,2716,2929,3152,3198,3241,3244,3303,3452,3530,3559,3586,3634,3812,3859,3917,3946,3991,4022,4130,4145,4231,4246,4501,4562,4668,4892,4941,4997,4998,5179,5405,5657,5791,5898,5927,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6722,6737,6752,6903,6980,7067,7096,7157,7186,7203,7354,7520,7562];
 ont2_15_init = [1398,2929,3559,3991,4998,6106,6722,7354,979,2026,3198,3634,4246,5657,6170,6903,7562,2417,3244,3917,454,5898,6455,7067,783,1067,272,3303,3946,4562,5927,6512,7203,888,1486,1876,2597,4668];

 temp_tcs = get_time_courses_matrix(datarun000_15, on_other15);
tc_fit = [];
final_params  =[];
for i = 1:length(on_other15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_other15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_other15) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)
 
 
 
 [tc nontc] = get_time_courses_matrix(datarun000_15,  on_other15); %or cellids
x = 1:1:30;
auc = [];
mx = [];
normval = [];
tcnormnorm = [];
tcnormauc = [];
tcnormmx = [];
for i = 1:length(on_other15) %or nonds
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
% [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, on_other15, 0);
% [mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
% NS2 = [];
% A32 = [];
% A256 = [];
% NS2 = NumSpikesCell';
% A32 = sum(NS2(find(StimComb(:,2) == 64),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
% A256 = sum(NS2(find(StimComb(:,2) == 256),:));
% close all;


datarun000_15 = get_interspikeinterval(datarun000_15, on_other15);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
for i = 1:length(on_other15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, on_other15(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

% close all;
% rstd = [];
% meanpix = [];
% cellind = get_cell_indices(datarun000_15, on_other15);
% stamat = cell(length(cellind),1);
%  
% for i = 1:length(cellind)
%     B = [];
%     B = datarun000_15.stas.rfs{cellind(i), 1}(:)';
%     meanpix(i) = mean(B);
%     rstd(i) = robust_std(B, [1]);
%     stamat{i,1} = zeros(size(datarun000_15.stas.rfs{cellind(i),1},1),size(datarun000_15.stas.rfs{cellind(i),1},2));
%     for j = 1:size(datarun000_15.stas.rfs{cellind(i),1},1)
%             for k = 1:size(datarun000_15.stas.rfs{cellind(i),1},2)
%                 if (datarun000_15.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
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
% radius = get_rf_fit_radius(datarun000_15, on_other15);

vc = [];
[C ia ib] = intersect(ont2_15_init, on_other15);
vc = ones(length(on_other15),1);
vc(ib) = 2;

 
[COEFF,SCORE] = princomp(tcnormnorm');
[COEFF1,SCORE1] = princomp(isinormnorm');
 
 X = [];
X(:,1) = TCParams.minval;
X(:,2) = TCParams.zerocrossing;
X(:,3) = SCORE1(:,1);
% X(:,1) = SCORE(:,2);
% X(:,2) = SCORE(:,3);
%X(:,3) = TCParams.mintim;
%X(:,4) = pm;
%X(:,5) = 1-(v(2,:)./v(1,:));
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, on_other15, tcnormnorm,0, vc);
on_other15(idx==2)
ont2_15 = on_other15(idx==2);
on_otherother15 = on_other15(idx==1);
ont2_15 = [272,454,783,888,979,1067,1398,1486,1876,2026,2417,2597,2929,3198,3244,3303,3559,3634,3917,3946,3991,4246,4562,4998,5657,5898,5927,6106,6170,6455,6512,6722,6903,7067,7203,7354,7562];
%on_otherother15 = [31,182,496,528,766,857,1385,1532,1549,1578,1591,1969,2136,2311,2328,2716,3152,3241,3452,3530,3586,3812,3859,4022,4130,4145,4231,4501,4668,4892,4941,4997,5179,5405,5791,6093,6125,6152,6155,6304,6363,6392,6439,6737,6752,6980,7096,7157,7186,7520];
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2classf', gcf)

%%
plot_rf_summaries(datarun000_15, ont2_15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2rf', gcf)

%%
plot_time_courses(datarun000_15,ont2_15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2tc', gcf)
%%
datarun000_15 = get_interspikeinterval(datarun000_15, ont2_15);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont2_15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, ont2_15(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2isifull', gcf)

%%
wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,ont2_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont2_15)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end

b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

spikesbytrials{1,1} = get_raster(datarun001_15.spikes{get_cell_indices(datarun001_15, ont2_15(1,1)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2psthcell272', gcf);
close;



figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
plot(binSize,psthnorm(1,:), 'b');
stairs([0 3 5 8 10],[0.6 0.55 0.5 0.55 0.55], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2psth2', gcf)











%%
datarun000_15 = get_autocorrelations(datarun000_15, ont2_15);

%%
x2 = 0:0.0005:0.1; 
%nonds - cells not ds, tc - their time courses
acf = [];
normvalacf = [];
acfnormnorm = [];
maxacf = [];
acfmax = [];
for i = 1:length(ont2_15) %or nonds
 acf(:,i) = datarun000_15.autocorrelation{get_cell_indices(datarun000_15, ont2_15(1,i)), 1}.probabilities;
 normvalacf(1, i) = norm( acf(:,i));
 maxacf(1,i) = max(acf(:,i));
end 
normvalacf = repmat(normvalacf, size(acf, 1), 1);
maxacf = repmat(maxacf, size(acf, 1), 1);
acfnormnorm = acf./normvalacf;
acfmax = acf./maxacf;


figure();
shadedErrorBar(x2(1:100),mean(acfnormnorm(1:100, :)'),std(acfnormnorm(1:100, :)'),'k');
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2acfhalf', gcf)
figure();
shadedErrorBar(x2,mean(acfnormnorm(:, :)'),std(acfnormnorm(:, :)'),'k');
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2acffull', gcf)



%%  ON T3
on_otherother15 = [31,182,496,528,766,857,1385,1532,1549,1578,1591,1969,2136,2311,2328,2716,3152,3241,3452,3530,3586,3812,3859,4022,4130,4145,4231,4501,4668,4892,4941,4997,5179,5405,5791,6093,6125,6152,6155,6304,6363,6392,6439,6737,6752,6980,7096,7157,7186,7520];
ont3_15_init = [1385 1549 2136 3452 4892 4941 5179 5405 6093 6980 7157 7520];

datarun000_15 = get_interspikeinterval(datarun000_15, on_otherother15);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(on_otherother15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, on_otherother15(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


  temp_tcs = get_time_courses_matrix(datarun000_15, on_otherother15);
tc_fit = [];
final_params  =[];
for i = 1:length(on_otherother15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_otherother15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_otherother15) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)


[tc nontc] = get_time_courses_matrix(datarun000_15, on_otherother15); %or cellids
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
for i = 1:length(on_otherother15) %or nonds
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

% [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, on_otherother15, 0);
% [mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
% NS2 = [];
% A32 = [];
% A256 = [];
% NS2 = NumSpikesCell';
% A32 = sum(NS2(find(StimComb(:,2) == 64),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
% A256 = sum(NS2(find(StimComb(:,2) == 256),:));
% close all;

% close all;
% rstd = [];
% meanpix = [];
% cellind = get_cell_indices(datarun000_15, on_otherother15);
% stamat = cell(length(cellind),1);
%  
% for i = 1:length(cellind)
%     B = [];
%     B = datarun000_15.stas.rfs{cellind(i), 1}(:)';
%     meanpix(i) = mean(B);
%     rstd(i) = robust_std(B, [1]);
%     stamat{i,1} = zeros(size(datarun000_15.stas.rfs{cellind(i),1},1),size(datarun000_15.stas.rfs{cellind(i),1},2));
%     for j = 1:size(datarun000_15.stas.rfs{cellind(i),1},1)
%             for k = 1:size(datarun000_15.stas.rfs{cellind(i),1},2)
%                 if (datarun000_15.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
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
% radius = get_rf_fit_radius(datarun000_15, on_otherother15);


[C ia ib] = intersect(ont3_15_init, on_otherother15);
vc = ones(length(on_otherother15),1);
vc(ib) = 2;


[COEFF1,SCORE1] = princomp(isinormnorm');
% [COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF,SCORE] = princomp(tcnormminn');

X = [];
X(:,1) = TCParams.minval;
X(:,2) = TCParams.zerocrossing;
X(:,3) = SCORE1(:,1);
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, on_otherother15, tcnormnorm,0,vc);
 ismember(ont3_15_init, on_otherother15)
on_otherother15(idx==2)
ont3_15 = on_otherother15(idx==2);
on_otherotherother15 = on_otherother15(idx==1);
%  ismember(ont3_15_init, ont3_15);
 ont3_15 = [1385,1549,2136,3452,4892,4941,5179,5405,6093,6980,7157,7520];
 on_otherotherother15 = [31,182,496,528,766,857,1532,1578,1591,1969,2311,2328,2716,3152,3241,3530,3586,3812,3859,4022,4130,4145,4231,4501,4668,4997,5791,6125,6152,6155,6304,6363,6392,6439,6737,6752,7096,7186];
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3classf', gcf)

%%
plot_rf_summaries(datarun000_15, ont3_15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3rf', gcf)

%%
close all;
plot_time_courses(datarun000_15,ont3_15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3tc', gcf)
%%
datarun000_15 = get_interspikeinterval(datarun000_15, ont3_15);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont3_15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, ont3_15(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3isifull', gcf)

%%
wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,ont3_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont3_15)
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
stairs([0 3 5 8 10],[0.6 0.55 0.5 0.55 0.55], 'Color', 'k', 'LineWidth',1);

set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3psth2', gcf)


spikesbytrials{1,1} = get_raster(datarun001_15.spikes{get_cell_indices(datarun001_15, ont3_15(1,1)), 1}, wh(1:8), 'tic_color', [0 0 0], 'axis_range', [0 10 0 8]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3psthcell1385', gcf);
close;






%%
datarun000_15 = get_autocorrelations(datarun000_15, ont3_15);

%%
x2 = 0:0.0005:0.1; 
%nonds - cells not ds, tc - their time courses
acf = [];
normvalacf = [];
acfnormnorm = [];
maxacf = [];
acfmax = [];
for i = 1:length(ont3_15) %or nonds
 acf(:,i) = datarun000_15.autocorrelation{get_cell_indices(datarun000_15, ont3_15(1,i)), 1}.probabilities;
 normvalacf(1, i) = norm( acf(:,i));
 maxacf(1,i) = max(acf(:,i));
end 
normvalacf = repmat(normvalacf, size(acf, 1), 1);
maxacf = repmat(maxacf, size(acf, 1), 1);
acfnormnorm = acf./normvalacf;
acfmax = acf./maxacf;


figure();
shadedErrorBar(x2(1:100),mean(acfnormnorm(1:100, :)'),std(acfnormnorm(1:100, :)'),'k');
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3acfhalf', gcf)
figure();
shadedErrorBar(x2,mean(acfnormnorm(:, :)'),std(acfnormnorm(:, :)'),'k');
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3acffull', gcf)



%%
offt1_15_init = [46,751,768,782,1051,1081,1381,1637,1652,1966,2177,2401,2686,3287,3589,3933,4006,4021,4324,4503,4591,5223,5446,5641,5658,5836,5896,6391,6451,6811,6976,7021,7306,7471];
off_15 =[46,62,76,94,226,347,407,514,586,649,751,753,768,782,843,872,889,901,991,1051,1081,1084,1098,1127,1130,1156,1191,1246,1310,1381,1384,1411,1427,1564,1581,1637,1652,1684,1817,1893,1908,1921,1966,1999,2177,2192,2206,2236,2253,2326,2343,2357,2373,2401,2462,2478,2521,2522,2555,2686,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3287,3589,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3933,3935,4006,4021,4097,4098,4188,4234,4235,4278,4294,4324,4442,4459,4486,4487,4503,4591,4697,4732,4771,4788,4864,4999,5071,5116,5148,5223,5433,5446,5464,5569,5641,5645,5658,5672,5703,5705,5733,5836,5851,5896,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6391,6422,6451,6589,6721,6811,6828,6886,6931,6976,6992,7021,7040,7069,7234,7278,7306,7471,7503,7517,7667];

 [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, off_15, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 64),:)); %CHANGE ACCORDING TO CORRECT TEMPORAL PERIOD
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
close all;


 temp_tcs = get_time_courses_matrix(datarun000_15, off_15);
tc_fit = [];
final_params  =[];
for i = 1:length(off_15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_15) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)



[tc nontc] = get_time_courses_matrix(datarun000_15, off_15); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
tcnormauc = [];
for i = 1:length(off_15) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
  auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 
normval = repmat(normval, size(tc, 1), 1);
auc = repmat(auc, size(tc, 1), 1);

tcnormnorm = tc./normval;
tcnormauc = tc./auc;



vc = [];
[C ia ib] = intersect(offt1_15_init, off_15);
vc = ones(length(off_15),1);
vc(ib) = 2;

%[COEFF,SCORE] = princomp(tcnormauc');


X = [];
X(:,1) = A32';
X(:,2) = A256';
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_15, tcnormnorm,0, vc);
off_15(idx==2)
offt1_15 = off_15(idx==2);
off_other15 = off_15(idx==1);
offt1_15 =[46,751,768,782,1051,1081,1381,1637,1652,1966,2177,2401,2686,3287,3589,3933,4006,4021,4324,4503,4591,5223,5446,5641,5658,5836,5896,6391,6451,6811,6976,7021,7306,7471];
%off_other15 = [62,76,94,226,347,407,514,586,649,753,843,872,889,901,991,1084,1098,1127,1130,1156,1191,1246,1310,1384,1411,1427,1564,1581,1684,1817,1893,1908,1921,1999,2192,2206,2236,2253,2326,2343,2357,2373,2462,2478,2521,2522,2555,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3935,4097,4098,4188,4234,4235,4278,4294,4442,4459,4486,4487,4697,4732,4771,4788,4864,4999,5071,5116,5148,5433,5464,5569,5645,5672,5703,5705,5733,5851,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6422,6589,6721,6828,6886,6931,6992,7040,7069,7234,7278,7503,7517,7667];
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T1/', 'offt1classf2', gcf)

%%
plot_rf_summaries(datarun000_15, offt1_15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T1/', 'offt1rf', gcf)

%%
close all;
plot_time_courses(datarun000_15,offt1_15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T1/', 'offt1tc', gcf)
%% cell 1652 for individual DG and pulse response

datarun000_15 = get_interspikeinterval(datarun000_15, offt1_15);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt1_15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, offt1_15(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T1/', 'offt1isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T1/', 'offt1isifull', gcf)

%%
wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,offt1_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt1_15)
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
plot(binSize,psthnorm(9,:), 'b');
stairs([0 3 5 8 10],[0.6 0.55 0.5 0.55 0.55], 'Color', 'k', 'LineWidth',1);

set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T1/', 'offt1psth2', gcf)


spikesbytrials{1,1} = get_raster(datarun001_15.spikes{get_cell_indices(datarun001_15, offt1_15(1,9)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T1/', 'offt1psthcell1652', gcf);
close;



%% off t1 extra cells
plot_rf_summaries(datarun000_15, [offt1_15], 'coordinates', 'monitor');
plot_time_courses(datarun000_15,[offt1_15], 'all', true, 'bw', true);

%% OFF T2
offt2_15_init = [76,753,1127,1921,5703,4442,4459,4486,1098,4864,1310,4999,1384,5071,1581,514,1893,5433,1908,5464,1999,5569,2192,5851,2357,586,2462,6034,2522,6140,2746,6196,2809,6422,2868,6589,2881,6721,3061,6828,3139,7234,3226,7503,3258,7667,347,872,3601,94,3692,3813,4098];
off_other15 = [62,76,94,226,347,407,514,586,649,753,843,872,889,901,991,1084,1098,1127,1130,1156,1191,1246,1310,1384,1411,1427,1564,1581,1684,1817,1893,1908,1921,1999,2192,2206,2236,2253,2326,2343,2357,2373,2462,2478,2521,2522,2555,2732,2746,2794,2809,2868,2881,2896,3061,3091,3139,3200,3226,3258,3286,3601,3635,3648,3679,3692,3695,3767,3813,3857,3889,3935,4097,4098,4188,4234,4235,4278,4294,4442,4459,4486,4487,4697,4732,4771,4788,4864,4999,5071,5116,5148,5433,5464,5569,5645,5672,5703,5705,5733,5851,6034,6064,6139,6140,6196,6257,6260,6286,6361,6376,6380,6422,6589,6721,6828,6886,6931,6992,7040,7069,7234,7278,7503,7517,7667];

[tc nontc] = get_time_courses_matrix(datarun000_15, off_other15); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
tcnormauc = [];
for i = 1:length(off_other15) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
  auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 
normval = repmat(normval, size(tc, 1), 1);
auc = repmat(auc, size(tc, 1), 1);

tcnormnorm = tc./normval;
tcnormauc = tc./auc;


temp_tcs = get_time_courses_matrix(datarun000_15, off_other15);
tc_fit = [];
final_params  =[];
for i = 1:length(off_other15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_other15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_other15) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0);
datarun000_15 = get_interspikeinterval(datarun000_15, off_other15);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_other15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, off_other15(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

radius = [];
radius = get_rf_fit_radius(datarun000_15, off_other15);

vc = [];
[C ia ib] = intersect(offt2_15_init, off_other15);
vc = ones(length(off_other15),1);
vc(ib) = 2;

%[COEFF,SCORE] = princomp(tcnormnorm');

%minval maxval zc dot maxt mint

X = [];
X(:,1) = TCParams.dot;
X(:,2) = TCParams.maxtim;
X(:,3) = TCParams.mintim;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_other15, tcnormnorm,0, vc);
off_other15(idx==2)
 %ismember(off_other15(idx==2),offt2_15_init)
 offt2_15 = off_other15(idx==2);
 off_otherother15 = off_other15(idx==1);%%
offt2_15 = [76,94,347,514,586,753,901,1098,1127,1310,1384,1581,1893,1908,1999,2192,2357,2462,2522,2746,2809,2868,2881,3061,3139,3226,3258,3601,3692,3813,4098,4442,4459,4486,4864,4999,5071,5433,5464,5569,5703,5851,6034,6140,6196,6422,6589,6721,6828,7234,7503,7667];
off_otherother15 = [62,226,407,649,843,872,889,991,1084,1130,1156,1191,1246,1411,1427,1564,1684,1817,1921,2206,2236,2253,2326,2343,2373,2478,2521,2555,2732,2794,2896,3091,3200,3286,3635,3648,3679,3695,3767,3857,3889,3935,4097,4188,4234,4235,4278,4294,4487,4697,4732,4771,4788,5116,5148,5645,5672,5705,5733,6064,6139,6257,6260,6286,6361,6376,6380,6886,6931,6992,7040,7069,7278,7517];

%% off t2 extra cells
plot_rf_summaries(datarun000_15, [offt2_15 3816 7400], 'coordinates', 'monitor');
plot_time_courses(datarun000_15,[offt2_15 3816 7400 ], 'all', true, 'bw', true);
%%
plot_rf_summaries(datarun000_15, [ 3816 7400], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2extrarf', gcf)

close all;
plot_time_courses(datarun000_15,[3816 7400], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2extratc', gcf)
%%
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),49, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
hold on;
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),49, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('Degree of Transience')
ylabel('Time of peak')
zlabel('Time of Trough')

legend('Cluster 1 - OFF T2 Cells','Cluster 2 - Other OFF Cells', 'Location','NW')
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
%%
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2classf2', gcf)

%%
plot_rf_summaries(datarun000_15, offt2_15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2rf', gcf)

%%
close all;
plot_time_courses(datarun000_15,offt2_15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2tc', gcf)

%% 
datarun000_15 = get_interspikeinterval(datarun000_15, offt2_15);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt2_15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, offt2_15(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2isifull', gcf)

%%
wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,offt2_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt2_15)
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
stairs([0 3 5 8 10],[0.6 0.55 0.5 0.55 0.55], 'Color', 'k', 'LineWidth',1);

set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2psth2', gcf)


spikesbytrials{1,1} = get_raster(datarun001_15.spikes{get_cell_indices(datarun001_15, offt2_15(1,3)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T2/', 'offt2psthcel347', gcf);
close;



%%
off_otherother15 = [62,226,407,649,843,872,889,991,1084,1130,1156,1191,1246,1411,1427,1564,1684,1817,1921,2206,2236,2253,2326,2343,2373,2478,2521,2555,2732,2794,2896,3091,3200,3286,3635,3648,3679,3695,3767,3857,3889,3935,4097,4188,4234,4235,4278,4294,4487,4697,4732,4771,4788,5116,5148,5645,5672,5705,5733,6064,6139,6257,6260,6286,6361,6376,6380,6886,6931,6992,7040,7069,7278,7517];
 offt4_15_init = [889,1191, 1684,2478,2896,3200,3648,3679,3767,3935,4188,4788, 4294,5645,6064,6361,7278, 1084, 3635,7040,7069,5705];


datarun000_15 = get_interspikeinterval(datarun000_15, off_otherother15);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherother15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, off_otherother15(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherother15); %or cellids
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
for i = 1:length(off_otherother15) %or nonds
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

temp_tcs = get_time_courses_matrix(datarun000_15, off_otherother15);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherother15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherother15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherother15) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);

% wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); %For 2012-15-31-1 dataset, that is how the triggers are arranges - need to change with dataset
% gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); %Might change with dataset
% [h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,off_otherother15), 0, '/0', wh, gr, 10, false);
% binSize = 0.1:0.1:10; %change depending on length of trial
% pulsePSTH = [];
% normvalpulsePSTH = [];
% pulsenormnormPSTH = [];
% maxpulse = [];
% maxpulsetime =[];
%  for a = 1:length(off_otherother15)
%  pulsePSTH(:,a) = sum(nhist{a,1})./50; %change depending on num of trials
% end
%  
% for i = 1:length(off_otherother15)
%  normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
% end
% normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
% pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;
% 
% [maxpulse maxpulsetime] = max(pulsePSTH);
% maxpulsetime = maxpulsetime*0.1;

 
[C ia ib] = intersect(offt4_15_init, off_otherother15);
vc = ones(length(off_otherother15),1);
vc(ib) = 2;


[COEFF1,SCORE1] = princomp(isinormnorm');
%[COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF2,SCORE2] = princomp(tcnormminn');

X = [];
X(:,1) = SCORE1(:,1);
X(:,2) = TCParams.minval;
X(:,3) = TCParams.mintim;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherother15, tcnormnorm,0,vc);
 offt4_15 = off_otherother15(idx==2);
 off_otherotherother15 = off_otherother15(idx==1);
 %off_otherother15(idx==2)
%ismember(offt4_15_init, off_otherother15(idx==2))
offt4_15 = [889,1084,1191,1684,2478,2896,3200,3635,3648,3679,3767,3889,3935,4188,4294,4788,5645,5705,6064,6361,7040,7069,7278];
%off_otherotherother15 = [62,226,407,649,843,872,991,1130,1156,1246,1411,1427,1564,1817,1921,2206,2236,2253,2326,2343,2373,2521,2555,2732,2794,3091,3286,3695,3857,4097,4234,4235,4278,4487,4697,4732,4771,5116,5148,5672,5733,6139,6257,6260,6286,6376,6380,6886,6931,6992,7517];
%% off t4 extra cells
figure();
plot_rf_summaries(datarun000_15, [offt4_15 4743 5042 4955 5131 ], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_15,[offt4_15 4743 5042 4955 5131], 'all', true, 'bw', true);

%%
plot_rf_summaries(datarun000_15, [4743 5042 4955 5131], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4extrarf', gcf)

close all;
plot_time_courses(datarun000_15,[4743 5042 4955 5131], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4extratc', gcf)
%%
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),49, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
hold on;
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),49, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('ISI - 1');
ylabel('Value of trough');
zlabel('Time of trough');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
%%
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4classf', gcf)

%%
plot_rf_summaries(datarun000_15, offt4_15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4rf', gcf)

%%
close all;
plot_time_courses(datarun000_15,offt4_15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4tc', gcf)

%% 
datarun000_15 = get_interspikeinterval(datarun000_15, offt4_15);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt4_15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, offt4_15(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4isifull', gcf)

%%
wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,offt4_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt4_15)
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
stairs([0 3 5 8 10],[0.6 0.55 0.5 0.55 0.55], 'Color', 'k', 'LineWidth',1);

set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4psth2', gcf)


spikesbytrials{1,1} = get_raster(datarun001_15.spikes{get_cell_indices(datarun001_15, offt4_15(1,5)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T4/', 'offt4psthcell2478', gcf);
close;




%%
 off_otherotherother15 = [62,226,407,649,843,872,991,1130,1156,1246,1411,1427,1564,1817,1921,2206,2236,2253,2326,2343,2373,2521,2555,2732,2794,3091,3286,3695,3857,4097,4234,4235,4278,4487,4697,4732,4771,5116,5148,5672,5733,6139,6257,6260,6286,6376,6380,6886,6931,6992,7517];
 offt3_15_init = [62,991,1156,4234,4278,4487,5733,6286,6931];

[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherotherother15); %or cellids
x = 1:1:30;
 normval = [];
 mx = [];
tcnormnorm = [];
tcnormmx = [];
for i = 1:length(off_otherotherother15) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);

tcnormnorm = tc./normval;
[mx mxt] = max(tc);
mx = repmat(mx, size(tc, 1), 1);
tcnormmx = tc./mx;

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

[TCParams] = time_course_parameters(tcfittednormnorm, 0);


vc = [];
[C ia ib] = intersect(offt3_15_init, off_otherotherother15);
vc = ones(length(off_otherotherother15),1);
vc(ib) = 2;
%minval maxval zc dot maxt mint

%[COEFF,SCORE] = princomp(tcnormmx');

X = [];
X(:,1) = TCParams.maxmingrad;
X(:,2) = TCParams.maxtim;
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherotherother15, tcnormnorm,0, vc);
%off_otherotherother15(idx==2)
offt3_15 = off_otherotherother15(idx==2);
off_otherotherotherother15 = off_otherotherother15(idx==1);%%
offt3_15 = [62,991,1130,1156,4234,4278,4487,4771,5733,6286,6886,6931];
off_otherotherotherother15 = [226,407,649,843,872,1246,1411,1427,1564,1817,1921,2206,2236,2253,2326,2343,2373,2521,2555,2732,2794,3091,3286,3695,3857,4097,4235,4697,4732,5116,5148,5672,6139,6257,6260,6376,6380,6992,7517];
%% off t3 extra cells
figure();
plot_rf_summaries(datarun000_15, [offt3_15 1069 2147 2313 6198 6875 7429 7637], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_15,[offt3_15 1069 2147 2313 6198 6875 7429 7637], 'all', true, 'bw', true);


%%
plot_rf_summaries(datarun000_15, [1069 2147 2313 6198 6875 7429 7637], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/', 'offt3rfextra', gcf)

%%
close all;
plot_time_courses(datarun000_15,[1069 2147 2313 6198 6875 7429 7637], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/', 'offt3extratc', gcf)



%% OFF T5
%offt5_15_init = [1246,2253,3695,5116,6260,4771,6380,226];
%plot_time_courses(datarun000_15, [1246,2253,3695,5116,6260,4771,6380,226], 'all', true, 'bw', true);
off_otherotherotherother15 = [226,407,649,843,872,1246,1411,1427,1564,1817,1921,2206,2236,2253,2326,2343,2373,2521,2555,2732,2794,3091,3286,3695,3857,4097,4235,4697,4732,5116,5148,5672,6139,6257,6260,6376,6380,6992,7517];
offt5_15_init = [1246,2253,3695,5116,6260]%4771,6380,226];


[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherotherotherother15); %or cellids
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
for i = 1:length(off_otherotherotherother15) %or nonds
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

temp_tcs = get_time_courses_matrix(datarun000_15, off_otherotherotherother15);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherotherother15)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherotherother15)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherotherother15) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);


[C ia ib] = intersect(offt5_15_init, off_otherotherotherother15);
vc = ones(length(off_otherotherotherother15),1);
vc(ib) = 2;


X = [];
X(:,1) = TCParams.maxtim;
X(:,2) = TCParams.mintim;
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherotherotherother15, tcnormnorm,0, vc);
offt5_15 = off_otherotherotherother15(idx==2);
off_otherotherotherotherother15 = off_otherotherotherother15(idx==1);%%
offt5_15 = [1246,2253,2373,3695,5116,6260];
off_otherotherotherotherother15 = [226,407,649,843,872,1411,1427,1564,1817,1921,2206,2236,2326,2343,2521,2555,2732,2794,3091,3286,3857,4097,4235,4697,4732,5148,5672,6139,6257,6376,6380,6992,7517];
%% off t5 extra cells
figure();
plot_rf_summaries(datarun000_15, [offt5_15 ], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_15,[offt5_15], 'all', true, 'bw', true);


%%
plot_rf_summaries(datarun000_15, offt5_15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/', 'offt5rf', gcf)

%%
close all;
plot_time_courses(datarun000_15,offt5_15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/', 'offt5tc', gcf)



%% off t6

off_otherotherotherotherother15 = [226,407,649,843,872,1411,1427,1564,1817,1921,2206,2236,2326,2343,2521,2555,2732,2794,3091,3286,3857,4097,4235,4697,4732,5148,5672,6139,6257,6376,6380,6992,7517];
offt6_15_init = [1817 2343 3286 6992];

datarun000_15 = get_interspikeinterval(datarun000_15, off_otherotherotherotherother15);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherotherotherotherother15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, off_otherotherotherotherother15(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherotherotherotherother15); %or cellids
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
for i = 1:length(off_otherotherotherotherother15) %or nonds
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
cellind = get_cell_indices(datarun000_15, off_otherotherotherotherother15);
stamat = cell(length(cellind),1);
 
for i = 1:length(cellind)
    B = [];
    B = datarun000_15.stas.rfs{cellind(i), 1}(:)';
    meanpix(i) = mean(B);
    rstd(i) = robust_std(B, [1]);
    stamat{i,1} = zeros(size(datarun000_15.stas.rfs{cellind(i),1},1),size(datarun000_15.stas.rfs{cellind(i),1},2));
    for j = 1:size(datarun000_15.stas.rfs{cellind(i),1},1)
            for k = 1:size(datarun000_15.stas.rfs{cellind(i),1},2)
                if (datarun000_15.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
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



 num_rgcs = length(off_otherotherotherotherother15);

% initialize some variables for the look
rf_areas = zeros(num_rgcs,1);
abs_mean_pixel_val = zeros(num_rgcs,1);
snrs = zeros(num_rgcs, 1);
contrast_index = zeros(num_rgcs, 1);
[Y, X] = meshgrid(1:1:40, 1:1:80);


for rgc = 1:num_rgcs
    
    temp_index= get_cell_indices(datarun000_15, off_otherotherotherotherother15(rgc));
    
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
        temp_rf = get_rf(datarun000_15, off_otherotherotherotherother15(rgc));
 
      
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


[C ia ib] = intersect(offt6_15_init, off_otherotherotherotherother15);
vc = ones(length(off_otherotherotherotherother15),1);
vc(ib) = 2;


X = [];
X(:,1) = SCORE(:,1);
X(:,2) = 1 - v(2,:)./v(1,:);
X(:,3) = contrast_index;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherotherotherotherother15, tcnormnorm,0, vc);
offt6_15 = off_otherotherotherotherother15(idx==2);
off_otherotherotherotherotherother15 = off_otherotherotherotherother15(idx==1);%%
offt6_15 = [1817,2343,3286,6992];
off_otherotherotherotherotherother15 = [226,407,649,843,872,1411,1427,1564,1921,2206,2236,2326,2521,2555,2732,2794,3091,3857,4097,4235,4697,4732,5148,5672,6139,6257,6376,6380,7517];



%% off t6 extra cells
figure();
plot_rf_summaries(datarun000_15, [offt6_15 4606 6485], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_15,[offt6_15 4606 6485], 'all', true, 'bw', true);


%%
plot_rf_summaries(datarun000_15, [4606 6485], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T6/', 'offt6extrarf', gcf)

close all;
plot_time_courses(datarun000_15,[4606 6485], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T6/', 'offt6extratc', gcf)
%%
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),49, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
hold on;
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),49, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('ISI - PC 1');
ylabel('1 - ratio')
zlabel('contrast index');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
%%
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T6/', 'offt6classf2', gcf)

%%
close all;
plot_rf_summaries(datarun000_15, on_otherotherother15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/', 'unclassfon', gcf)

%%
close all;
plot_time_courses(datarun000_15,off_otherotherotherotherotherother15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/', 'unclassofftc', gcf)

%% 
datarun000_15 = get_interspikeinterval(datarun000_15, offt6_15);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt6_15) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, offt6_15(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T6/', 'offt6isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T6/', 'offt6isifull', gcf)

%%
wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,offt6_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt6_15)
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
stairs([0 3 5 8 10],[0.6 0.55 0.5 0.55 0.55], 'Color', 'k', 'LineWidth',1);

set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T6/', 'offt6psth2', gcf)


spikesbytrials{1,1} = get_raster(datarun001_15.spikes{get_cell_indices(datarun001_15, offt6_15(1,2)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/OFF/T6/', 'offt6psthcell2343', gcf);
close;

%%
plot_rf_summaries(datarun000_15, ont3_15, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3rf', gcf)

%%
close all;
plot_time_courses(datarun000_15,ont3_15, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3tc', gcf)
%%

%%
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002, datarun002.cell_ids, 0);
NumSpikesCellIn = NumSpikesCell(get_cell_indices(datarun002, cellids), :);
[magin dsindexin magmaxin magavein anglein rhoin thetain numin Uin Vin] = dscellanalysis(NumSpikesCellIn, StimComb);

plot(magin{1,1}, magin{3,1}, 'o')
addpath('/Users/sravi/matlab/DS cell analysis/2012-10-15-0/Data001PulsePlots');
pathname = ['/Users/sravi/matlab/DS cell analysis/2012-10-15-0/Data001PulsePlots/'  num2str(chos(1,i)) '-' num2str(datarun.cell_ids(chos(1,i)))];
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(gcf,'-dpdf', pathname);

%Only extract information from 324 cells instead of 1336

get_cell_indices(datarun002, cellids)
NumSpikesCellIn = NumSpikesCell(get_cell_indices(datarun002, cellids), :);
[magin dsindexin magmaxin magavein anglein] = dscellanalysis(NumSpikesCellIn, StimComb);

plot(mag{1,1}, mag{2,1},'o');
plot(mag{1,1}, mag{3,1},'o');
plot(mag{2,1}, mag{3,1},'o');

hist(dsindex{1,1},100);
hist(dsindex{2,1},100);
hist(dsindex{3,1},100);

plot(angle{1,1},angle{2,1}, '+');
plot(angle{1,1},angle{3,1}, '+');
plot(angle{2,1},angle{3,1}, '+');

hist(magmax, 50);
hist(magave, 50)




plot(magin{1,1}, magin{2,1},'o');
plot(magin{1,1}, magin{3,1},'o');
plot(magin{2,1}, magin{3,1},'o');

hist(dsindexin{1,1},100);
hist(dsindexin{2,1},100);
hist(dsindexin{3,1},100);

plot(dsindexin{1,1}, dsindexin{2,1}, 'o')
plot(dsindexin{1,1}, dsindexin{3,1}, 'o')
plot(dsindexin{2,1}, dsindexin{3,1}, 'o')

plot(anglein{1,1},anglein{2,1}, '+');
plot(anglein{1,1},anglein{3,1}, '+');
plot(anglein{2,1},anglein{3,1}, '+');

hist(magmaxin, 50);
hist(magavein, 50) 


ginput(1);
Y = get(get(gca,'Children'),'YData');
X = get(get(gca,'Children'),'XData');
h = impoly;
accepted_pos = wait(h);
hold on;
drawPolygon(accepted_pos(:,1), accepted_pos(:,2));
in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));
plot(X(in),Y(in),'r+');
cid = 1:length(X);
chos = cid(in);

%DS Cells
%chos = ds cells index wrt 324 cells
ds = cellids(1, chos); %cell ids of ds cells

NumSpikesCellDS = NumSpikesCell(get_cell_indices(datarun002, ds), :);
[magds dsindexds magmaxds magaveds angleds rhods thetads numds] = dscellanalysis(NumSpikesCellDS, StimComb);

plot(magds{1,1}, magds{2,1},'o');
plot(magds{1,1}, magds{3,1},'o');
plot(magds{2,1}, magds{3,1},'o');

hist(dsindexds{1,1},100);
hist(dsindexds{2,1},100);
hist(dsindexds{3,1},100);

plot(dsindexds{1,1}, dsindexds{2,1}, 'o')
plot(dsindexds{1,1}, dsindexds{3,1}, 'o')
plot(dsindexds{2,1}, dsindexds{3,1}, 'o')

plot(angleds{1,1},angleds{2,1}, '+');
plot(angleds{1,1},angleds{3,1}, '+');
plot(angleds{2,1},angleds{3,1}, '+');

hist(magmaxds, 50);
hist(magaveds, 50) 






plot(angle{1,1}(1, chos(1,:)), angle{2,1}(1, chos(1,:)), '+')
plot(angle{1,1}(1, chos(1,:)), angle{3,1}(1, chos(1,:)), '+')
plot(angle{2,1}(1, chos(1,:)), angle{3,1}(1, chos(1,:)), '+')
plot(mag{2,1}(1, chos(1,:)), mag{3,1}(1, chos(1,:)), '+')

A = C12(1,chos);

plot_rf_summaries(datarun000, ds,'label', true)

plot_rf_portraits(datarun000,ds)

plot_rfs(datarun000, ds)

plot_time_courses(datarun000, ds, 1, 1)

datarun000 = get_autocorrelations(datarun000, ds);

plot_autocorrelograms(datarun000, ds, 'foa', 0);

 Locb(1, chos)

chos = get_cell_indices(datarun002,A);

pulse_analysis(datarun001, chos)

 
%ON type 1: 
A = [62 991 1156 4234 4278 4487 5733 6286];
% A = [ 2555  4732 6257 ];
% A = [226 3091 6380];
A = [889 1191 1684 2478 2896 3200 3635 3648 3679 3767 3889 3935 4188 4294 4788 5645 6064 6361 7040 7278];
plot_rf_summaries(datarun000, A,'label', true)
plot_rf_portraits(datarun000,A)
plot_rfs(datarun000, A)
plot_time_courses(datarun000, A, 1, 1)
datarun000 = get_autocorrelations(datarun000, A);
plot_autocorrelograms(datarun000, A, 'foa', 0);


plot_rf_summaries(datarun000, [1817 2343 3286], 'label', true)
plot_time_courses(datarun000, [1817 2343 3286], 1, 1)
plot_autocorrelograms(datarun000, [1817 2343 3286], 'foa', 0);


A = [4 154 692 860 1111 1339 1786 1877 1892 1969 2011 2161 2461 2851 3002 3152 3319 3422 3691 3994 4096 4427 4774 5073 5088 5359 5567 5791 5853 5941 6542 6826 7096 7261 7442 7487 7532];
plot_rf_summaries(datarun000, A,'label', true)
plot_rf_portraits(datarun000,A)
plot_rfs(datarun000, A)
plot_time_courses(datarun000, A, 1, 1)
datarun000 = get_autocorrelations(datarun000, A);
plot_autocorrelograms(datarun000, A, 'foa', 0);



[Lia2,Locb2] = ismember(C12, cellids2);



[datarun000] = load_dsdata('/Analysis/sravi/Rat/WildType/2013-04-01-0/data000-3600-7200s/', 'data000-map/data000-map', 0, 0, 1);
[datarun001] = load_dsdata('/Analysis/sravi/Rat/WildType/2013-04-01-0/data000-3600-7200s/', 'data001-map/data001-map', 0, 0, 0);
[datarun002] = load_dsdata('/Analysis/sravi/Rat/WildType/2013-04-01-0/data000-3600-7200s/', 'data002-map/data002-map', 1, '/stimuli/s02', 0);


[datarun004] = load_dsdata('/Analysis/sravi/Rat/WildType/2013-04-01-0/data000-3600-7200s/', 'data004-map/data004-map', 1, '/stimuli/s04', 0);


cellids = intersect((intersect((intersect(datarun000.cell_ids, datarun001.cell_ids)), datarun004.cell_ids)), datarun002.cell_ids);


[NumSpikesCell, StimComb] = get_spikescellstim(datarun004, datarun004.cell_ids, 0);
NumSpikesCellIn = NumSpikesCell(get_cell_indices(datarun004, cellids), :);
[mag dsindex magmax magave angle rho theta num U V] = dscellanalysis(NumSpikesCell, StimComb);
[magin dsindexin magmaxin magavein anglein rhoin thetain numin Uin Vin] = dscellanalysis(NumSpikesCellIn, StimComb);
plot(magin{2,1}, magin{5,1}, 'o');
plot(magin{3,1}, magin{5,1}, 'o');

(magin{5,1}(1,ds03))./(magin{3,1}(1,ds03))

all_rasterplots_pdf(NumSpikesCell,StimComb,datarun004, rho,theta,U,V,num, get_cell_indices(datarun004, cellids), 1,2)

%% 2012-10-15 extra cells
%%on t3
plot_rf_summaries(datarun000_15, [ont3_15,  800,1368,1996,2237,2941,3017,3185,3410,3755, 6242,4966,3811,3979,4143,4337,7592], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_15, [ont3_15,  800,1368,1996,2237,2941,3017,3185,3410,3755,6242,4966,3811,3979,4143,4337,7592], 'all', true, 'bw', true);


plot_rf_summaries(datarun000_15, [800,1368,1996,2237,2941,3017,3185,3410,3755, 6242,4966,3811,3979,4143,4337,7592], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3extrarf', gcf)

plot_time_courses(datarun000_15,[800,1368,1996,2237,2941,3017,3185,3410,3755, 6242,4966,3811,3979,4143,4337,7592], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T3/', 'ont3extratc', gcf)


%% ont2 - 10-15
plot_rf_summaries(datarun000_15, [ont2_15], 'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_15, [ont2_15 227 678 1501 2101 2266 2356 2671 2821 2867 4322 4534 5390 5492 5583 6001 7171 7399 7638], 'all', true, 'bw', true);
plot_rf_summaries(datarun000_15, [ont2_15 227 678 1501 2101 2266 2356 2671 2821 2867 4322 4534 5390 5492 5583 6001 7171 7399 7638], 'coordinates', 'monitor');

%%
plot_rf_summaries(datarun000_15, [227 678 1501 2101 2266 2356 2671 2821 2867 4322 4534 5390 5492 5583 6001 7171 7399 7638], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2extrarf', gcf)

plot_time_courses(datarun000_15,[227 678 1501 2101 2266 2356 2671 2821 2867 4322 4534 5390 5492 5583 6001 7171 7399 7638], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont2extratc', gcf)




%% ont1 - 2012-10-15

plot_time_courses(datarun000_15, [ont1_15 4638 4893 5751 7307 4248 ], 'all', true, 'bw', true);
figure();
plot_rf_summaries(datarun000_15, [ont1_15 4638 4893 5751 7307 4248 ], 'coordinates', 'monitor');

% 4638 4893 5751 7307 4248 6181


%% nov 26th finding extra cells for off t4, t3, t5, t6
temp_tcs = get_time_courses_matrix(datarun000_15, datarun000_15.cell_ids);
tc_fit = [];
final_params  =[];
for i = 1:length(datarun000_15.cell_ids)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(datarun000_15.cell_ids)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(datarun000_15.cell_ids) %or nonds
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

%%
[tc nontc] = get_time_courses_matrix(datarun000_15, datarun000_15.cell_ids); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(datarun000_15.cell_ids) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;
%%
on_init15 = [4,31,154,182,272,424,454,496,528,692,766,783,857,860,888,979,1067,1111,1339,1385,1398,1486,1532,1549,1578,1591,1786,1876,1877,1892,1969,2011,2026,2136,2161,2208,2311,2328,2417,2461,2597,2716,2851,2929,3002,3152,3198,3241,3244,3303,3319,3422,3452,3530,3559,3586,3634,3691,3721,3812,3859,3917,3946,3991,3994,4022,4096,4130,4145,4231,4246,4279,4427,4501,4562,4668,4774,4892,4941,4997,4998,5073,5088,5104,5179,5359,5405,5567,5657,5791,5853,5898,5927,5941,6093,6106,6125,6152,6155,6170,6304,6363,6392,6439,6455,6512,6542,6722,6737,6752,6826,6903,6980,7067,7096,7157,7186,7203,7261,7354,7442,7487,7520,7532,7562];
[C ia ib] = intersect(on_init15, datarun000_15.cell_ids);
vc = ones(length(datarun000_15.cell_ids),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_15, datarun000_15.cell_ids, tcnormnorm,0, vc);
on_allsnr15_all = datarun000_15.cell_ids(idx ==2);
off_allsnr15_all = datarun000_15.cell_ids(idx ==1);

%%

c = get_cell_indices(datarun000_15, on_allsnr15_all);
snronall = [];
for i = 1:length(c)
    r1 = sort(datarun000_15.stas.rfs{c(1,i),1}(:)', 'descend');
    snronall(1,i) = mean(r1(1:4))./std(r1);
end
on_15_all = on_allsnr15_all(snronall > (mean(snronall) - 2.5*std(snronall)));

% hax=axes; 
% hold on;
% hist(snronall)
% SP= mean(snronall) - 2.5*std(snronall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('On cutoff')

c = get_cell_indices(datarun000_15, off_allsnr15_all);
snroffall = [];
for i = 1:length(c)
    r1 = sort(datarun000_15.stas.rfs{c(1,i),1}(:)', 'descend');
    snroffall(1,i) = mean(r1(1:4))./std(r1);
end
off_15_all = off_allsnr15_all(snroffall > (mean(snroffall) - 2.5*std(snroffall)));

%%

off_otherother15_all = off_15_all(~ismember(off_15_all, [offt1_15_all offt2_15_all]));
 offt4_15_init = [889,1191, 1684,2478,2896,3200,3648,3679,3767,3935,4188,4788, 4294,5645,6064,6361,7278, 1084, 3635,7040,7069,5705];


datarun000_15 = get_interspikeinterval(datarun000_15, off_otherother15_all);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherother15_all) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, off_otherother15_all(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherother15_all); %or cellids
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
for i = 1:length(off_otherother15_all) %or nonds
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

temp_tcs = get_time_courses_matrix(datarun000_15, off_otherother15_all);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherother15_all)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherother15_all)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherother15_all) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);

% wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); %For 2012-15-31-1 dataset, that is how the triggers are arranges - need to change with dataset
% gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); %Might change with dataset
% [h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,off_otherother15_all), 0, '/0', wh, gr, 10, false);
% binSize = 0.1:0.1:10; %change depending on length of trial
% pulsePSTH = [];
% normvalpulsePSTH = [];
% pulsenormnormPSTH = [];
% maxpulse = [];
% maxpulsetime =[];
%  for a = 1:length(off_otherother15_all)
%  pulsePSTH(:,a) = sum(nhist{a,1})./50; %change depending on num of trials
% end
%  
% for i = 1:length(off_otherother15_all)
%  normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
% end
% normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
% pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;
% 
% [maxpulse maxpulsetime] = max(pulsePSTH);
% maxpulsetime = maxpulsetime*0.1;

 
[C ia ib] = intersect(offt4_15_init, off_otherother15_all);
vc = ones(length(off_otherother15_all),1);
vc(ib) = 2;


[COEFF1,SCORE1] = princomp(isinormnorm');
%[COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF2,SCORE2] = princomp(tcnormminn');

X = [];
X(:,1) = SCORE1(:,1);
X(:,2) = TCParams.minval;
X(:,3) = TCParams.mintim;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherother15_all, tcnormnorm,0,vc);
 %offt4_15 = off_otherother15_all(idx==2);
 off_otherotherother15_all = off_otherother15_all(idx==1);
 offt4_15_all = off_otherother15_all(idx==2);
ismember(offt4_15_all, cellids_15)

%%
offt4_15_all = [offt4_15 4743 5042 4955 5131 ];

%%

off_otherotherother15_all = [62,226,407,649,708,796,843,872,991,1008,1069,1130,1156,1246,1276,1411,1415,1427,1564,1595,1816,1817,1879,1921,2147,2206,2222,2236,2253,2299,2313,2326,2343,2373,2521,2555,2732,2794,2898,3091,3286,3483,3527,3558,3695,3842,3857,3889,4097,4173,4234,4235,4278,4487,4550,4606,4697,4732,4771,4940,4985,5116,5148,5341,5643,5672,5705,5733,5793,6139,6198,6257,6260,6286,6376,6380,6485,6875,6886,6931,6992,7052,7054,7429,7517,7547,7637];
offt3_15_init = [62,991,1156,4234,4278,4487,5733,6286,6931];

[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherotherother15_all); %or cellids
x = 1:1:30;
 normval = [];
 mx = [];
tcnormnorm = [];
tcnormmx = [];
for i = 1:length(off_otherotherother15_all) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);

tcnormnorm = tc./normval;
[mx mxt] = max(tc);
mx = repmat(mx, size(tc, 1), 1);
tcnormmx = tc./mx;

temp_tcs = get_time_courses_matrix(datarun000_15, off_otherotherother15_all);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherother15_all)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherother15_all)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherother15_all) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0);


vc = [];
[C ia ib] = intersect(offt3_15_init, off_otherotherother15_all);
vc = ones(length(off_otherotherother15_all),1);
vc(ib) = 2;
%minval maxval zc dot maxt mint

%[COEFF,SCORE] = princomp(tcnormmx');

X = [];
X(:,1) = TCParams.maxmingrad;
X(:,2) = TCParams.maxtim;
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherotherother15_all, tcnormnorm,0, vc);
offt3_15_all = off_otherotherother15_all(idx==2);
off_otherotherotherother15_all = off_otherotherother15_all(idx==1);%

%%
offt3_15_all = [offt3_15 1069 2147 2313 6198 6875 7429 7637];

%%
off_otherotherotherother15_all = [226,407,649,708,796,843,872,1008,1246,1276,1411,1415,1427,1564,1595,1816,1817,1879,1921,2206,2222,2236,2253,2299,2326,2343,2373,2521,2555,2732,2794,2898,3091,3286,3483,3527,3558,3695,3842,3857,3889,4097,4173,4235,4550,4606,4697,4732,4940,4985,5116,5148,5341,5643,5672,5705,5793,6139,6257,6260,6376,6380,6485,6992,7052,7054,7517,7547];
offt5_15_init = [1246,2253,3695,5116,6260];%4771,6380,226];

[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherotherotherother15_all); %or cellids
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
for i = 1:length(off_otherotherotherother15_all) %or nonds
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

temp_tcs = get_time_courses_matrix(datarun000_15, off_otherotherotherother15_all);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherotherother15_all)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherotherother15_all)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherotherother15_all) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);


[C ia ib] = intersect(offt5_15_init, off_otherotherotherother15_all);
vc = ones(length(off_otherotherotherother15_all),1);
vc(ib) = 2;


X = [];
X(:,1) = TCParams.maxtim;
X(:,2) = TCParams.mintim;
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherotherotherother15_all, tcnormnorm,0, vc);
offt5_15_all = off_otherotherotherother15_all(idx==2);
off_otherotherotherotherother15_all = off_otherotherotherother15_all(idx==1);%%

%%
offt5_15_all = [offt5_15 ];

%%

off_otherotherotherotherother15_all = [226,407,649,708,796,843,872,1008,1276,1411,1415,1427,1564,1595,1816,1817,1879,1921,2206,2222,2236,2299,2326,2343,2373,2521,2555,2732,2794,2898,3091,3286,3483,3527,3558,3842,3857,3889,4097,4173,4235,4550,4606,4697,4732,4940,4985,5148,5341,5672,5705,5793,6139,6257,6376,6380,6485,6992,7052,7054,7517,7547];
offt6_15_init = [1817 2343 3286 6992];

datarun000_15 = get_interspikeinterval(datarun000_15, off_otherotherotherotherother15_all);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherotherotherotherother15_all) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, off_otherotherotherotherother15_all(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_15, off_otherotherotherotherother15_all); %or cellids
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
for i = 1:length(off_otherotherotherotherother15_all) %or nonds
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
cellind = get_cell_indices(datarun000_15, off_otherotherotherotherother15_all);
stamat = cell(length(cellind),1);
 
for i = 1:length(cellind)
    B = [];
    B = datarun000_15.stas.rfs{cellind(i), 1}(:)';
    meanpix(i) = mean(B);
    rstd(i) = robust_std(B, [1]);
    stamat{i,1} = zeros(size(datarun000_15.stas.rfs{cellind(i),1},1),size(datarun000_15.stas.rfs{cellind(i),1},2));
    for j = 1:size(datarun000_15.stas.rfs{cellind(i),1},1)
            for k = 1:size(datarun000_15.stas.rfs{cellind(i),1},2)
                if (datarun000_15.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
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



 num_rgcs = length(off_otherotherotherotherother15_all);

% initialize some variables for the look
rf_areas = zeros(num_rgcs,1);
abs_mean_pixel_val = zeros(num_rgcs,1);
snrs = zeros(num_rgcs, 1);
contrast_index = zeros(num_rgcs, 1);
[Y, X] = meshgrid(1:1:40, 1:1:80);


for rgc = 1:num_rgcs
    
    temp_index= get_cell_indices(datarun000_15, off_otherotherotherotherother15_all(rgc));
    
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
        temp_rf = get_rf(datarun000_15, off_otherotherotherotherother15_all(rgc));
 
      
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


[C ia ib] = intersect(offt6_15_init, off_otherotherotherotherother15_all);
vc = ones(length(off_otherotherotherotherother15_all),1);
vc(ib) = 2;


X = [];
X(:,1) = SCORE(:,1);
X(:,2) = 1 - v(2,:)./v(1,:);
X(:,3) = contrast_index;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_15, off_otherotherotherotherother15_all, tcnormnorm,0, vc);
offt6_15_all = off_otherotherotherotherother15_all(idx==2);
off_otherotherotherotherotherother15_all = off_otherotherotherotherother15_all(idx==1);%%

%% end of nov 26th finding extra cells for off t4, t3, t5, t6

offt6_15_all = [offt6_15 4606 6485]

%% jan 21st 2015
wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); 
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,ont1_15), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psth = [];
psthind = [];
for i = 1:length(ont1_15)
    psthind = sum(nhist{i,1})/25;
    psthind = psthind./max(psthind);
    psth(i,:) = psthind;
end

plot(binSize,mean(psth), 'k');
hold on;
plot(binSize,psth(3,:), 'b');
plot(binSize,mean(psth)+std(psth), 'Color', [0.5 0.5 0.5]);
plot(binSize,mean(psth)-std(psth), 'Color', [0.5 0.5 0.5]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont1pulse', gcf)



%%
plot(binSize,mean(psth), 'k');
hold on;
plot(binSize,psth(3,:), 'b');
plot(binSize,mean(psth)+std(psth), 'Color', [0.5 0.5 0.5]);
plot(binSize,mean(psth)-std(psth), 'Color', [0.5 0.5 0.5]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-15-0/ON/T2/', 'ont1pulse', gcf)

%% Drifitng grating PSTH code
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, ont3_15, 0);

ind = find(StimComb(:,2)==256);
zrind = find(StimComb(:,3)==0);
zeroind = intersect(ind, zrind);
psthall = [];
for i = 1:length(ont3_15)
    psthzero = [];
    psth = [];
    [T, psthzero, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont3_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,zeroind)), 'stop', 8, 'bin_size', 0.1);
    psthzerocut = psthzero(21:61);
    for k = 1:length(ind)
        psth01 = [];
        [T, psth01, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont3_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,ind(k))), 'stop', 8, 'bin_size', 0.1);
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
dgresp_15(1,:) = mean(psthall);
dgresp_15(2,:) = mean(psthall);
dgresp_15(3,:) = mean(psthall);

plot(dgresp_15','DisplayName','dgresp')
legend('t1', 't2', 't3');
title('2012-10-15');

%%
Fs = 10;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 41;                     % Length of signal
t = (0:L-1)*T;                % Time vector
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(dgresp_31(3,:),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
spec = 2*abs(Y(1:NFFT/2+1));
% Plot single-sided amplitude spectrum.
plot(f,spec) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

frameRate = 120; %Our frame rate is 60 %Daniel FR = 120;
f1 = frameRate / 256; %initialize f1 and f2 given temporal period of stimulus
f2 = f1*2;

fourprop_31(3,1) = spec(ismember(f, f1));
fourprop_31(3,2) = spec(ismember(f, f2));
fourprop_31(3,3) = spec(ismember(f, 0));

set(gca,'XTickLabel',{'t1', 't2', 't3'})
legend('f1', 'f2', 'dc');
title('10-10')



%%


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

%% Orientation Selective Cells
% Calculating Orientation Selective Index = Pref - Null / Pref + Null
% Null direction defined first, followed by Preferred (either Orthogonal to
% Null or Direction with maximum firing rate)
% Plotting OSI - Histogram or 2D Plot of Speeds

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, nonds_31, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V spave] = dscellanalysis(NumSpikesCell, StimComb);

%%
minangle = [];
maxangle = [];
minAxis = [];
maxAxis = [];
minAxFir = [];
maxAxFir = [];
minFir = cell(2,1);
maxFir = cell(2,1);
[CC,minAxis] = min(spave{2,1}');
for i = 1:length(nonds_31)
    minangle(1,i) = theta{2,1}(i,minAxis(1,i));
    minangle(2,i) = minangle(1,i)+pi;
    if (minangle(2,i) < 0 || minangle(2,i) > pi + 3*pi/4)
        minangle(2,i) = minangle(1,i) - pi;
    end
    minAxis(2,i) = find(minangle(2,i) == theta{2,1}(1,:));
    maxangle(1,i) = minangle(1,i)+pi./2;
    if (maxangle(1,i) < 0 || maxangle(1,i) > pi + 3*pi/4)
        maxangle(1,i) = minangle(1,i)-pi./2;
    end
    maxAxis(1,i) = find(maxangle(1,i) == theta{2,1}(1,:));
    maxangle(2,i) = maxangle(1,i)+pi;
        if (maxangle(2,i) < 0 || maxangle(2,i) > pi + 3*pi/4)
                maxangle(2,i) = maxangle(1,i)-pi;
        end
      maxAxis(2,i) = find(maxangle(2,i) == theta{2,1}(1,:));
      minAxFir(i) = (spave{2,1}(i,minAxis(1,i))+spave{2,1}(i,minAxis(2,i)))/2;
      maxAxFir(i) = (spave{2,1}(i,maxAxis(1,i))+spave{2,1}(i,maxAxis(2,i)))/2;
end

scatter(minAxFir,maxAxFir);
hold on;
xlabel('minimum firing');
ylabel('maximum firing')
title('S 64 T 256')
%plot(0:1:250, 0:1:250);
hold off;

%%
figure()
hist(maxAxFir./minAxFir, 50)
xlabel('maximum firing / minimum firing');
title('S 64 T 256')
%%
minFir{2,1} = minAxFir;
maxFir{2,1} = maxAxFir;

%%

minangle = [];
maxangle = [];
minAxis = [];
maxAxis = [];
minAxFir = [];
maxAxFir = [];
[CC,minAxis] = min(spave{1,1}')
for i = 1:length(nonds_31)
    minangle(1,i) = theta{1,1}(i,minAxis(i));
    minangle(2,i) = minangle(1,i)+pi;
    if (minangle(2,i) < 0 || minangle(2,i) > pi + 3*pi/4)
        minangle(2,i) = minangle(1,i) - pi;
    end
    minAxis(2,i) = find(minangle(2,i) == theta{1,1}(1,:));
    maxangle(1,i) = minangle(1,i)+pi./2;
    if (maxangle(1,i) < 0 || maxangle(1,i) > pi + 3*pi/4)
        maxangle(1,i) = minangle(1,i)-pi./2;
    end
    maxAxis(1,i) = find(maxangle(1,i) == theta{1,1}(1,:));
    maxangle(2,i) = maxangle(1,i)+pi;
        if (maxangle(2,i) < 0 || maxangle(2,i) > pi + 3*pi/4)
                maxangle(2,i) = maxangle(1,i)-pi;
        end
      maxAxis(2,i) = find(maxangle(2,i) == theta{1,1}(1,:));
      minAxFir(i) = (spave{1,1}(i,minAxis(1,i))+spave{1,1}(i,minAxis(2,i)))/2;
      maxAxFir(i) = (spave{1,1}(i,maxAxis(1,i))+spave{1,1}(i,maxAxis(2,i)))/2;
end

scatter(minAxFir,maxAxFir);
hold on;
xlabel('minimum firing');
ylabel('maximum firing')
title('S 64 T 32')
%plot(0:1:250, 0:1:250);
hold off;
figure();
hist(maxAxFir./minAxFir,50)
xlabel('maximum firing / minimum firing');
title('S 64 T 32')

minFir{1,1} = minAxFir;
maxFir{1,1} = maxAxFir;

%%
xxx = maxFir{1,1}./minFir{1,1};
yyy = maxFir{2,1}./minFir{2,1};
scatter(xxx,yyy);
xlabel('maximum firing / minimum firing - T 32');
ylabel('maximum firing / minimum firing - T 256')
title('S 64 T 32-256')
figure()
hist(xxx.*yyy,50);
xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');

%%
xxx = [];
yyy = [];
xxx = 1-(minFir{1,1}./maxFir{1,1});
xxx(xxx<0) = 0;
yyy = 1-(minFir{2,1}./maxFir{2,1});
yyy(yyy<0) = 0;
scatter(xxx,yyy);
xlabel('1 - (null/pref) - T 32');
ylabel('1- (null/pref) - T 256')
title('S 64 T 32-256')
figure()
hist(xxx, 50);
figure()
hist(yyy,50);
figure()
hist(xxx.*yyy,50);
xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');




