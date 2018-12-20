%% load dataruns - drifting grating, white noise, pulses

addpath('/Users/sneharavi/Documents/MATLAB/Classification/');
addpath('/Users/sneharavi/Documents/MATLAB/DS cell analysis/');


[datarun002_31] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-31-0/data000-1800-7200/', 'data002-map/data002-map', 1,'stimuli/s02',0);
[datarun000_31] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-31-0/data000-1800-7200/', 'data000-map/data000-map', 0,0,1);
[datarun001_31] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-31-0/data000-1800-7200/', 'data001-map/data001-map', 0,0,0);

%% DS

cellids_31 = intersect((intersect(datarun000_31.cell_ids, datarun001_31.cell_ids)), datarun002_31.cell_ids); %intersecting cell ids through all 3 stimuli
%cellids_31 = [4,49,78,92,108,259,272,286,316,362,378,391,393,421,469,484,543,556,559,572,591,618,619,647,681,724,781,783,785,903,917,962,976,991,1006,1037,1067,1081,1083,1126,1172,1203,1248,1277,1306,1368,1381,1430,1442,1471,1487,1501,1576,1577,1595,1606,1670,1683,1712,1726,1731,1741,1756,1772,1773,1816,1862,1878,1922,1952,1954,1996,2028,2074,2086,2087,2101,2118,2146,2176,2177,2178,2356,2371,2419,2433,2494,2506,2536,2539,2581,2582,2597,2613,2656,2659,2687,2719,2747,2825,2851,2856,2867,2868,2884,2896,2897,2899,2942,2971,2973,3019,3046,3049,3076,3077,3121,3125,3215,3244,3260,3274,3289,3317,3482,3571,3616,3661,3843,3871,3887,3905,3931,3946,4036,4038,4066,4112,4113,4126,4142,4156,4159,4171,4172,4186,4204,4248,4277,4351,4366,4383,4384,4413,4518,4548,4578,4625,4681,4685,4697,4712,4713,4714,4727,4730,4771,4786,4846,4876,4892,4983,4999,5000,5026,5028,5042,5058,5072,5089,5117,5132,5148,5210,5240,5267,5281,5297,5328,5386,5401,5431,5462,5495,5506,5552,5596,5629,5642,5746,5748,5765,5791,5792,5794,5841,5853,5866,5897,5926,5943,5956,5961,6002,6003,6031,6033,6062,6091,6106,6122,6211,6212,6213,6242,6271,6289,6290,6376,6422,6483,6497,6499,6541,6542,6587,6590,6602,6617,6631,6646,6662,6692,6708,6752,6753,6781,6784,6796,6812,6827,6888,6901,6933,6962,6964,6993,7051,7067,7083,7098,7114,7156,7157,7159,7203,7291,7324,7336,7442,7444,7472,7502,7503,7518,7562,7593,7640,7668,7669];
[tc nontc] = get_time_courses_matrix(datarun000_31, cellids_31); 
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(cellids_31) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval; %normalize time courses by norm

%DS cells
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, cellids_31, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
ds_init31 = [92,316,484,559,619,1037,1442,1487,1683,1996,2074,2118,2433,2613,2942,3019,3077,3125,3260,4159,4548,4625,5000,5148,5210,6242,6290,6602,6708,6781,7083,7159,7503,7518]; %initial conditions for ds cells

[C ia ib] = intersect(ds_init31, cellids_31);
vc = ones(length(cellids_31),1);
vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

close all;
X = [];
N = [];
p = [];
X(:,1) = log(mag{1,1})';
X(:,2) = log(mag{2,1})';
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_31, cellids_31, tcnormnorm,0, vc);

ds_31 = [];
ds_31 = cellids_31(idx==2);
nonds_31 = cellids_31(idx==1);
%lr = sum(p(idx==2,2))/length(ds_31)
%lr = sum(p(idx==1,1))/length(nonds_31)

%ds_31 = [92,316,484,559,619,1037,1442,1487,1683,1996,2074,2118,2433,2613,2942,3019,3077,3125,3260,4159,4548,4625,5000,5148,5210,6242,6290,6602,6708,6781,7083,7159,7503,7518];
%nonds_31 = [4,49,78,108,259,272,286,362,378,391,393,421,469,543,556,572,591,618,647,681,724,781,783,785,903,917,962,976,991,1006,1067,1081,1083,1126,1172,1203,1248,1277,1306,1368,1381,1430,1471,1501,1576,1577,1595,1606,1670,1712,1726,1731,1741,1756,1772,1773,1816,1862,1878,1922,1952,1954,2028,2086,2087,2101,2146,2176,2177,2178,2356,2371,2419,2494,2506,2536,2539,2581,2582,2597,2656,2659,2687,2719,2747,2825,2851,2856,2867,2868,2884,2896,2897,2899,2971,2973,3046,3049,3076,3121,3215,3244,3274,3289,3317,3482,3571,3616,3661,3843,3871,3887,3905,3931,3946,4036,4038,4066,4112,4113,4126,4142,4156,4171,4172,4186,4204,4248,4277,4351,4366,4383,4384,4413,4518,4578,4681,4685,4697,4712,4713,4714,4727,4730,4771,4786,4846,4876,4892,4983,4999,5026,5028,5042,5058,5072,5089,5117,5132,5240,5267,5281,5297,5328,5386,5401,5431,5462,5495,5506,5552,5596,5629,5642,5746,5748,5765,5791,5792,5794,5841,5853,5866,5897,5926,5943,5956,5961,6002,6003,6031,6033,6062,6091,6106,6122,6211,6212,6213,6271,6289,6376,6422,6483,6497,6499,6541,6542,6587,6590,6617,6631,6646,6662,6692,6752,6753,6784,6796,6812,6827,6888,6901,6933,6962,6964,6993,7051,7067,7098,7114,7156,7157,7203,7291,7324,7336,7442,7444,7472,7502,7562,7593,7640,7668,7669];

%% ON - OFF Cells
temp_tcs = get_time_courses_matrix(datarun000_31, nonds_31);
tc_fit = [];
final_params  =[];
for i = 1:length(nonds_31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(nonds_31) %fit time course
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(nonds_31) %or nonds
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


[tc nontc] = get_time_courses_matrix(datarun000_31, nonds_31); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(nonds_31) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

on_init31 = [4,78,108,286,362,543,618,647,781,785,903,917,962,991,1067,1081,1203,1248,1306,1381,1471,1595,1606,1712,1756,1878,1952,2086,2087,2146,2371,2419,2494,2536,2582,2687,2719,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4038,4113,4156,4172,4186,4248,4277,4384,4685,4697,4713,4714,4876,4892,4983,5028,5042,5132,5240,5267,5297,5328,5401,5552,5642,5746,5765,5792,5794,5866,5897,5956,5961,6002,6031,6211,6289,6497,6590,6646,6662,6692,6752,6753,6784,6812,6962,7067,7098,7114,7336,7442,7472,7640,7669];
[C ia ib] = intersect(on_init31, nonds_31);
vc = ones(length(nonds_31),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_31, nonds_31, tcnormnorm,0, vc);
on_allsnr31 = nonds_31(idx ==2); %idx might change so be careful - on might be idx 2 and off idx 1
off_allsnr31 = nonds_31(idx ==1);
%on_allsnr31 = [4,78,108,286,362,543,618,647,781,785,903,917,962,991,1067,1081,1203,1248,1306,1381,1471,1595,1606,1712,1756,1878,1952,2086,2087,2146,2371,2419,2494,2536,2582,2687,2719,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4038,4113,4156,4172,4186,4248,4277,4384,4685,4697,4713,4714,4876,4892,4983,5028,5042,5132,5240,5267,5297,5328,5401,5552,5642,5746,5765,5792,5794,5866,5897,5956,5961,6002,6031,6211,6289,6497,6590,6646,6662,6692,6752,6753,6784,6812,6962,7067,7098,7114,7336,7442,7472,7640,7669];
%off_allsnr31 = [49,259,272,378,391,393,421,469,556,572,591,681,724,783,976,1006,1083,1126,1172,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2101,2176,2177,2178,2356,2506,2539,2581,2597,2656,2659,2825,2867,2868,2884,2899,2971,2973,3046,3049,3121,3274,3317,3482,3571,3616,3661,3871,3905,3946,4036,4066,4112,4126,4142,4171,4204,4351,4366,4383,4413,4518,4578,4681,4712,4727,4730,4771,4786,4846,4999,5026,5058,5072,5089,5117,5281,5386,5431,5462,5495,5506,5596,5629,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6122,6212,6213,6271,6376,6422,6483,6499,6541,6542,6587,6617,6631,6796,6827,6888,6901,6933,6964,6993,7051,7156,7157,7203,7291,7324,7444,7502,7562,7593,7668];
%% SNR CUTOFF
c = get_cell_indices(datarun000_31, on_allsnr31);
snronall = [];
for i = 1:length(c)
    r1 = sort(datarun000_31.stas.rfs{c(1,i),1}(:)', 'descend');
    snronall(1,i) = mean(r1(1:4))./std(r1);
end
on_31 = on_allsnr31(snronall > (mean(snronall) - 2.5*std(snronall)));
% 
% hax=axes; 
% hold on;
% hist(snronall)
% SP= mean(snronall) - 2.5*std(snronall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('On cutoff')

c = get_cell_indices(datarun000_31, off_allsnr31);
snroffall = [];
for i = 1:length(c)
    r1 = sort(datarun000_31.stas.rfs{c(1,i),1}(:)', 'descend');
    snroffall(1,i) = mean(r1(1:4))./std(r1);
end
off_31 = off_allsnr31(snroffall > (mean(snroffall) - 2.5*std(snroffall)));

% hax=axes; 
% hold on;
% hist(snroffall)
% SP= mean(snroffall) - 2.5*std(snroffall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('Off cutoff')

snroff_31 = off_allsnr31(snroffall < (mean(snroffall) - 2.5*std(snroffall)));
snron_31 = on_allsnr31(snronall < (mean(snronall) - 2.5*std(snronall)));
% 
% snron_31 = [1081 4038];
% snroff_31 = [1172 3121 4126];
on_31 = [4,78,108,286,362,543,618,647,781,785,903,917,962,991,1067,1203,1248,1306,1381,1471,1595,1606,1712,1756,1878,1952,2086,2087,2146,2371,2419,2494,2536,2582,2687,2719,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4113,4156,4172,4186,4248,4277,4384,4685,4697,4713,4714,4876,4892,4983,5028,5042,5132,5240,5267,5297,5328,5401,5552,5642,5746,5765,5792,5794,5866,5897,5956,5961,6002,6031,6211,6289,6497,6590,6646,6662,6692,6752,6753,6784,6812,6962,7067,7098,7114,7336,7442,7472,7640,7669];
% off_31 = [49,259,272,378,391,393,421,469,556,572,591,681,724,783,976,1006,1083,1126,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2101,2176,2177,2178,2356,2506,2539,2581,2597,2656,2659,2825,2867,2868,2884,2899,2971,2973,3046,3049,3274,3317,3482,3571,3616,3661,3871,3905,3946,4036,4066,4112,4142,4171,4204,4351,4366,4383,4413,4518,4578,4681,4712,4727,4730,4771,4786,4846,4999,5026,5058,5072,5089,5117,5281,5386,5431,5462,5495,5506,5596,5629,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6122,6212,6213,6271,6376,6422,6483,6499,6541,6542,6587,6617,6631,6796,6827,6888,6901,6933,6964,6993,7051,7156,7157,7203,7291,7324,7444,7502,7562,7593,7668];
%% ON - OFF after SNR check

nondssnr_31 = [on_31 off_31];

temp_tcs = get_time_courses_matrix(datarun000_31, nondssnr_31);
tc_fit = [];
final_params  =[];
for i = 1:length(nondssnr_31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(nondssnr_31) %fit time course
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(nondssnr_31) %or nonds
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


[tc nontc] = get_time_courses_matrix(datarun000_31, nondssnr_31); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(nondssnr_31) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

on_initsnr31 = on_31;
[C ia ib] = intersect(on_initsnr31, nondssnr_31);
vc = ones(length(nondssnr_31),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_31, nondssnr_31, tcnormnorm,0, vc);
onon_31 = nondssnr_31(idx ==2); %idx might change so be careful - on might be idx 2 and off idx 1
offoff_snr31 = nondssnr_31(idx ==1);


    plot(X(idx==2,1),X(idx==2,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [1 0 0], 'MarkerFaceColor' , [1 0.7 0.8]);
    hold on;
        plot(X(idx==1,1),X(idx==1,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]); 
            legend('Cluster 1 - On Cells','Cluster 2 - Off Cells', 'Location','NW')

set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON-OFF/', 'classfaftersnr', gcf)



%% ON T1 ------

%2012-10-31-0: 102 on 21 cells  all t1, 1t1 missing

on_31 = [4,78,108,286,362,543,618,647,781,785,903,917,962,991,1067,1203,1248,1306,1381,1471,1595,1606,1712,1756,1878,1952,2086,2087,2146,2371,2419,2494,2536,2582,2687,2719,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4113,4156,4172,4186,4248,4277,4384,4685,4697,4713,4714,4876,4892,4983,5028,5042,5132,5240,5267,5297,5328,5401,5552,5642,5746,5765,5792,5794,5866,5897,5956,5961,6002,6031,6211,6289,6497,6590,6646,6662,6692,6752,6753,6784,6812,6962,7067,7098,7114,7336,7442,7472,7640,7669];
ont1_31_init = [4,78,286,618,917,1248,1712,2087,2719,4172,4248,4277,4697,4713,4983,5401,5552,5897,6289,6590,7067,7472];
%unsure t1s not added: 2536  4892  5042 5132 - cuz they distort the clustering - add many other cells

temp_tcs = get_time_courses_matrix(datarun000_31, on_31);
tc_fit = [];
final_params  =[];
for i = 1:length(on_31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_31)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_31) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)


[tc nontc] = get_time_courses_matrix(datarun000_31, on_31); %or cellids
x = 1:1:30;
auc = [];
mx = [];
normval = [];
tcnormnorm = [];
tcnormauc = [];
tcnormmx = [];

for i = 1:length(on_31) %or nonds
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

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, on_31, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 64),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
close all;

vc = [];
[C ia ib] = intersect(ont1_31_init, on_31);
vc = ones(length(on_31),1);
vc(ib) = 2;

[COEFF,SCORE] = princomp(tcnormauc');


X = [];
X(:,1) = A32';
X(:,2) = A256';
%X(:,3) = SCORE(:,3);
X(:,3) = TCParams.dot';
%X(:,1) = TCParams.mintim';
%X(:,2) = TCParams.maxtim';
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, on_31, tcnormnorm,0, vc);
on_31(idx==2)
ont1_31 = on_31(idx==2);
on_other31 = on_31(idx==1);
ont1_31 = [4,78,286,618,917,1248,1712,2087,2719,4172,4248,4277,4697,4713,4983,5401,5552,5897,6289,6590,7067,7472];
% on_other31 = [108,362,543,647,781,785,903,962,991,1067,1203,1306,1381,1471,1595,1606,1756,1878,1952,2086,2146,2371,2419,2494,2536,2582,2687,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4113,4156,4186,4384,4685,4714,4876,4892,5028,5042,5132,5240,5267,5297,5328,5642,5746,5765,5792,5794,5866,5956,5961,6002,6031,6211,6497,6646,6662,6692,6752,6753,6784,6812,6962,7098,7114,7336,7442,7640,7669];
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T1/', 'ont1classf', gcf)

%%
plot_rf_summaries(datarun000_31, ont1_31, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T1/', 'ont1rf', gcf)

%%
plot_time_courses(datarun000_31,ont1_31, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T1/', 'ont1tc', gcf)
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T1/', 'ont1tccomp', gcf)

%%
datarun000_31 = get_interspikeinterval(datarun000_31, ont1_31);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont1_31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, ont1_31(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T1/', 'ont1isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T1/', 'ont1isifull', gcf)
%%
wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); 
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,ont1_31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont1_31)
    psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind);
    psthnorm(i,:) = psthindnorm;
end
b = 2; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

spikesbytrials{1,1} = get_raster(datarun001_31.spikes{get_cell_indices(datarun001_31, ont1_31(1,5)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T1/', 'ont1psthcell917', gcf);
close;


figure();
shadedErrorBar(binSize,mean(psthnorm(:, :)),std(psthnorm(:, :)),'k');
hold on;
plot(binSize,psthnorm(5,:), 'b');
stairs([0 3 5 8 10],[0.4 0.35 0.3 0.35 0.35], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T1/', 'ont1psth2', gcf)
%% ON T2
on_other31 = [108,362,543,647,781,785,903,962,991,1067,1203,1306,1381,1471,1595,1606,1756,1878,1952,2086,2146,2371,2419,2494,2536,2582,2687,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4113,4156,4186,4384,4685,4714,4876,4892,5028,5042,5132,5240,5267,5297,5328,5642,5746,5765,5792,5794,5866,5956,5961,6002,6031,6211,6497,6646,6662,6692,6752,6753,6784,6812,6962,7098,7114,7336,7442,7640,7669];
ont2_31_init = [108,362,647,903,962,991,1203,1756,2086,2371,2419,2582,2856,2897,3076,3289,3843,3887,4384,4685,4876,5028,5297,5642,5746,5956,6002,6031,6497,6752,6812,7098,7442,7669];

  temp_tcs = get_time_courses_matrix(datarun000_31, on_other31);
tc_fit = [];
final_params  =[];
for i = 1:length(on_other31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_other31)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_other31) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)
 
 
 
 [tc nontc] = get_time_courses_matrix(datarun000_31,  on_other31); %or cellids
x = 1:1:30;
auc = [];
mx = [];
normval = [];
tcnormnorm = [];
tcnormauc = [];
tcnormmx = [];
for i = 1:length(on_other31) %or nonds
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

% [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, on_other31, 0);
% [mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
% NS2 = [];
% A32 = [];
% A256 = [];
% NS2 = NumSpikesCell';
% A32 = sum(NS2(find(StimComb(:,2) == 64),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
% A256 = sum(NS2(find(StimComb(:,2) == 256),:));
% close all;


datarun000_31 = get_interspikeinterval(datarun000_31, on_other31);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
for i = 1:length(on_other31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, on_other31(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

% close all;
% rstd = [];
% meanpix = [];
% cellind = get_cell_indices(datarun000_31, on_other31);
% stamat = cell(length(cellind),1);
%  
% for i = 1:length(cellind)
%     B = [];
%     B = datarun000_31.stas.rfs{cellind(i), 1}(:)';
%     meanpix(i) = mean(B);
%     rstd(i) = robust_std(B, [1]);
%     stamat{i,1} = zeros(size(datarun000_31.stas.rfs{cellind(i),1},1),size(datarun000_31.stas.rfs{cellind(i),1},2));
%     for j = 1:size(datarun000_31.stas.rfs{cellind(i),1},1)
%             for k = 1:size(datarun000_31.stas.rfs{cellind(i),1},2)
%                 if (datarun000_31.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
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
%                 %hold on;
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
% radius = get_rf_fit_radius(datarun000_31, on_other31);


vc = [];
[C ia ib] = intersect(ont2_31_init, on_other31);
vc = ones(length(on_other31),1);
vc(ib) = 2;

 
[COEFF,SCORE] = princomp(tcnormauc');
[COEFF1,SCORE1] = princomp(isinormnorm');
 %maxval minval zc tc auc12 3 tc norm 1-2 pm dot
 X = [];
X(:,1) = TCParams.minval;
X(:,2) = TCParams.zerocrossing;
X(:,3) =  SCORE1(:,1);
%X(:,1) = SCORE(:,1);
%X(:,2) = SCORE(:,2);
%X(:,3) = SCORE(:,3);
%X(:,3) = pm;
% X(:,1) = 1-(v(2,:)./v(1,:));
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, on_other31, tcnormnorm,0, vc);
on_other31(idx==2)
ont2_31 = on_other31(idx==2);
on_otherother31 = on_other31(idx==1);
ont2_31 =[108,362,647,903,962,991,1203,1756,2086,2371,2419,2582,2856,2897,3076,3289,3843,3887,3931,4384,4685,4876,5028,5297,5642,5746,5956,6002,6031,6497,6752,6812,7098,7442,7669];
%on_otherother31 =[543,781,785,1067,1306,1381,1471,1595,1606,1878,1952,2146,2494,2536,2687,2747,2851,2896,3215,3244,4113,4156,4186,4714,4892,5042,5132,5240,5267,5328,5765,5792,5794,5866,5961,6211,6646,6662,6692,6753,6784,6962,7114,7336,7640];
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2classf', gcf)

%%
plot_rf_summaries(datarun000_31, ont2_31, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2rf', gcf)

%%
plot_time_courses(datarun000_31,ont2_31, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2tc', gcf)

%%
datarun000_31 = get_interspikeinterval(datarun000_31, ont2_31);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont2_31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, ont2_31(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2isifull', gcf)


%%
wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); 
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,ont2_31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont2_31)
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
plot(binSize,psthnorm(2,:), 'b');
stairs([0 3 5 8 10],[0.8 0.75 0.7 0.75 0.75], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2psth', gcf);
close;

spikesbytrials{1,1} = get_raster(datarun001_31.spikes{get_cell_indices(datarun001_31, ont2_31(1,2)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2psthcell362', gcf);
close;

%% ON T3
on_otherother31 =[543,781,785,1067,1306,1381,1471,1595,1606,1878,1952,2146,2494,2536,2687,2747,2851,2896,3215,3244,4113,4156,4186,4714,4892,5042,5132,5240,5267,5328,5765,5792,5794,5866,5961,6211,6646,6662,6692,6753,6784,6962,7114,7336,7640];
ont3_31_init = [1067 1306 1595 1878 2851 5328 6692 6753 ];

datarun000_31 = get_interspikeinterval(datarun000_31, on_otherother31);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(on_otherother31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, on_otherother31(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


  temp_tcs = get_time_courses_matrix(datarun000_31, on_otherother31);
tc_fit = [];
final_params  =[];
for i = 1:length(on_otherother31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(on_otherother31)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(on_otherother31) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)


[tc nontc] = get_time_courses_matrix(datarun000_31, on_otherother31); %or cellids
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
for i = 1:length(on_otherother31) %or nonds
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

% [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, on_otherother31, 0);
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
% cellind = get_cell_indices(datarun000_31, on_otherother31);
% stamat = cell(length(cellind),1);
%  
% for i = 1:length(cellind)
%     B = [];
%     B = datarun000_31.stas.rfs{cellind(i), 1}(:)';
%     meanpix(i) = mean(B);
%     rstd(i) = robust_std(B, [1]);
%     stamat{i,1} = zeros(size(datarun000_31.stas.rfs{cellind(i),1},1),size(datarun000_31.stas.rfs{cellind(i),1},2));
%     for j = 1:size(datarun000_31.stas.rfs{cellind(i),1},1)
%             for k = 1:size(datarun000_31.stas.rfs{cellind(i),1},2)
%                 if (datarun000_31.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
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



[C ia ib] = intersect(ont3_31_init, on_otherother31);
vc = ones(length(on_otherother31),1);
vc(ib) = 2;


[COEFF1,SCORE1] = princomp(isinormnorm');
% [COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF,SCORE] = princomp(tcnormminn');

X = [];
X(:,1) = TCParams.minval;
X(:,2) = TCParams.zerocrossing
X(:,3) = SCORE1(:,1);
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, on_otherother31, tcnormnorm,0,vc);
 ismember(ont3_31_init, on_otherother31)
on_otherother31(idx==2)
ont3_31 = on_otherother31(idx==2);
on_otherotherother31 = on_otherother31(idx==1);
%  ismember(ont3_31_init, ont3_31);
ont3_31 = [1067,1306,1595,1878,2494,2851,5240,5328,6211,6692,6753];
%1067 and 1595 in on t3 have different EIs but exact same RF
on_otherotherother31 =[543,781,785,1381,1471,1606,1952,2146,2536,2687,2747,2896,3215,3244,4113,4156,4186,4714,4892,5042,5132,5267,5765,5792,5794,5866,5961,6646,6662,6784,6962,7114,7336,7640];
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3classf', gcf)

%%
plot_rf_summaries(datarun000_31, ont3_31, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3rf', gcf)

%%
close all;
plot_time_courses(datarun000_31,ont3_31, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3tc', gcf)

 %%
datarun000_31 = get_interspikeinterval(datarun000_31, ont3_31);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(ont3_31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, ont3_31(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3isifull', gcf)

%%


wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); 
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,ont3_31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(ont3_31)
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
plot(binSize,psthnorm(2,:), 'b');
stairs([0 3 5 8 10],[0.5 0.45 0.4 0.45 0.45], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3psth2', gcf)
close;

spikesbytrials{1,1} = get_raster(datarun001_31.spikes{get_cell_indices(datarun001_31, ont3_31(1,2)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3psthcell1306', gcf);
close;



%%
datarun000_31 = get_autocorrelations(datarun000_31, ont3_31);

x2 = 0:0.0005:0.1; 
%nonds - cells not ds, tc - their time courses
acf = [];
normvalacf = [];
acfnormnorm = [];
maxacf = [];
acfmax = [];
for i = 1:length(ont3_31) %or nonds
 acf(:,i) = datarun000_31.autocorrelation{get_cell_indices(datarun000_31, ont3_31(1,i)), 1}.probabilities;
 normvalacf(1, i) = norm( acf(:,i));
 maxacf(1,i) = max(acf(:,i));
end 
normvalacf = repmat(normvalacf, size(acf, 1), 1);
maxacf = repmat(maxacf, size(acf, 1), 1);
acfnormnorm = acf./normvalacf;
acfmax = acf./maxacf;


figure();
shadedErrorBar(x2(1:100),mean(acfnormnorm(1:100, :)'),std(acfnormnorm(1:100, :)'),'k');
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3acfhalf', gcf)
figure();
shadedErrorBar(x2,mean(acfnormnorm(:, :)'),std(acfnormnorm(:, :)'),'k');
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3acffull', gcf)



 %% OFF T1 -------
off_31 = [49,259,272,378,391,393,421,469,556,572,591,681,724,783,976,1006,1083,1126,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2101,2176,2177,2178,2356,2506,2539,2581,2597,2656,2659,2825,2867,2868,2884,2899,2971,2973,3046,3049,3274,3317,3482,3571,3616,3661,3871,3905,3946,4036,4066,4112,4142,4171,4204,4351,4366,4383,4413,4518,4578,4681,4712,4727,4730,4771,4786,4846,4999,5026,5058,5072,5089,5117,5281,5386,5431,5462,5495,5506,5596,5629,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6122,6212,6213,6271,6376,6422,6483,6499,6541,6542,6587,6617,6631,6796,6827,6888,6901,6933,6964,6993,7051,7156,7157,7203,7291,7324,7444,7502,7562,7593,7668]; 
offt1_31_init = [378,572,1083,1126,2101,2176,2506,2825,3049,3317,3616,4112,4351,4413,4846,5058,5462,5629,6122,6422,6617,6933,7051,7156,7668];
 
 [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, off_31, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 64),:)); %CHANGE ACCORDING TO CORRECT TEMPORAL PERIOD
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
close all;


 temp_tcs = get_time_courses_matrix(datarun000_31, off_31);
tc_fit = [];
final_params  =[];
for i = 1:length(off_31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_31)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_31) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0)



[tc nontc] = get_time_courses_matrix(datarun000_31, off_31); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
tcnormauc = [];
for i = 1:length(off_31) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
  auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 
normval = repmat(normval, size(tc, 1), 1);
auc = repmat(auc, size(tc, 1), 1);

tcnormnorm = tc./normval;
tcnormauc = tc./auc;



vc = [];
[C ia ib] = intersect(offt1_31_init, off_31);
vc = ones(length(off_31),1);
vc(ib) = 2;

%[COEFF,SCORE] = princomp(tcnormauc');


X = [];
X(:,1) = A32';
X(:,2) = A256';
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_31, tcnormnorm,0, vc);
off_31(idx==2)
offt1_31 = off_31(idx==2);
off_other31 = off_31(idx==1);
offt1_31 =[378,572,1083,1126,2101,2176,2506,2825,3049,3317,4112,4351,4413,4846,5058,5462,5629,6122,6422,6617,6933,7051,7156];
%off_other31 =[49,259,272,391,393,421,469,556,591,681,724,783,976,1006,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2177,2178,2356,2539,2581,2597,2656,2659,2867,2868,2884,2899,2971,2973,3046,3274,3482,3571,3616,3661,3871,3905,3946,4036,4066,4142,4171,4204,4366,4383,4518,4578,4681,4712,4727,4730,4771,4786,4999,5026,5072,5089,5117,5281,5386,5431,5495,5506,5596,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6212,6213,6271,6376,6483,6499,6541,6542,6587,6631,6796,6827,6888,6901,6964,6993,7157,7203,7291,7324,7444,7502,7562,7593,7668];
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1classf2', gcf)

%%
plot_rf_summaries(datarun000_31, offt1_31, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1rf', gcf)

%%
close all;
plot_time_courses(datarun000_31,offt1_31, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1tc', gcf)

%%
plot_rf_summaries(datarun000_31, [offt1_31 6318], 'coordinates', 'monitor');
plot_time_courses(datarun000_31,[offt1_31  6318 ], 'all', true, 'bw', true);

%%
plot_rf_summaries(datarun000_31, 6318, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1extrarf', gcf)

close all;
plot_time_courses(datarun000_31,6318, 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1extratc', gcf)


 %%
datarun000_31 = get_interspikeinterval(datarun000_31, offt1_31);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt1_31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, offt1_31(1,i)), 1}.probabilities;
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
%save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1isifull', gcf)

%%
wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); 
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,offt1_31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt1_31)
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
stairs([0 3 5 8 10],[0.5 0.45 0.4 0.45 0.45], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1psth2', gcf)
close;

spikesbytrials{1,1} = get_raster(datarun001_31.spikes{get_cell_indices(datarun001_31, offt1_31(1,1)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T1/', 'offt1psthcell378', gcf);
close;





%% OFF T2

offt2_31_init =[5791,4142,1670,421,469,556,976,1006,1277,1576,1726,1922,2028,2356,2581,2656,2868,3046,3274,3482,3571,3871,4171,4366,4712,4999,5026,5281,5386,5748,6003,6271,6587,6631,6827,6964,7291,7502];
off_other31 =[49,259,272,391,393,421,469,556,591,681,724,783,976,1006,1277,1368,1430,1501,1576,1577,1670,1726,1731,1741,1772,1773,1816,1862,1922,1954,2028,2177,2178,2356,2539,2581,2597,2656,2659,2867,2868,2884,2899,2971,2973,3046,3274,3482,3571,3616,3661,3871,3905,3946,4036,4066,4142,4171,4204,4366,4383,4518,4578,4681,4712,4727,4730,4771,4786,4999,5026,5072,5089,5117,5281,5386,5431,5495,5506,5596,5748,5791,5841,5853,5926,5943,6003,6033,6062,6091,6106,6212,6213,6271,6376,6483,6499,6541,6542,6587,6631,6796,6827,6888,6901,6964,6993,7157,7203,7291,7324,7444,7502,7562,7593,7668];

[tc nontc] = get_time_courses_matrix(datarun000_31, off_other31); %or cellids
x = 1:1:30;
normval = [];
auc = [];
tcnormnorm = [];
tcnormauc = [];
for i = 1:length(off_other31) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
  auc(1,i) = trapz(abs(x), abs(tc(:,i))); %Calculate Area Under Curve forall time courses
end 
normval = repmat(normval, size(tc, 1), 1);
auc = repmat(auc, size(tc, 1), 1);

tcnormnorm = tc./normval;
tcnormauc = tc./auc;


temp_tcs = get_time_courses_matrix(datarun000_31, off_other31);
tc_fit = [];
final_params  =[];
for i = 1:length(off_other31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_other31)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_other31) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 0);
datarun000_31 = get_interspikeinterval(datarun000_31, off_other31);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_other31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, off_other31(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

radius = [];
radius = get_rf_fit_radius(datarun000_31, off_other31);

vc = [];
[C ia ib] = intersect(offt2_31_init, off_other31);
vc = ones(length(off_other31),1);
vc(ib) = 2;

[COEFF,SCORE] = princomp(tcnormnorm');

% maxt mint minval maxval  dot 
X = [];
X(:,1) = TCParams.dot;
X(:,2) = TCParams.maxtim;
X(:,3) = TCParams.mintim;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_other31, tcnormnorm,0, vc);
off_other31(idx==2)
offt2_31 = off_other31(idx==2);
off_otherother31 = off_other31(idx==1);
offt2_31 = [421,469,556,976,1006,1277,1576,1670,1726,1922,2028,2356,2581,2656,2868,3046,3274,3482,3571,3661,3871,4142,4171,4366,4578,4712,4999,5026,5281,5386,5748,5791,5853,6003,6106,6271,6587,6631,6827,6964,7291];
%off_otherother31 = [49,259,272,391,393,591,681,724,783,1368,1430,1501,1577,1731,1741,1772,1773,1816,1862,1954,2177,2178,2539,2597,2659,2867,2884,2899,2971,2973,3616,3905,3946,4036,4066,4204,4383,4518,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5926,5943,6033,6062,6091,6212,6213,6376,6483,6499,6541,6542,6796,6888,6901,6993,7157,7203,7324,7444,7502,7562,7593,7668];
%%
plot_rf_summaries(datarun000_31, [offt2_31 ], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T2/', 'offt2rf', gcf)

close all;
plot_time_courses(datarun000_31,[offt2_31], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T2/', 'offt2tc', gcf)


%% off t2 extra cells
plot_rf_summaries(datarun000_31, [offt2_31 3826 6874 ], 'coordinates', 'monitor');
plot_time_courses(datarun000_31,[offt2_31 3826 6874 ], 'all', true, 'bw', true);

%%
plot_rf_summaries(datarun000_31, [3826 6874 ], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T2/', 'offt2extrarf', gcf)

close all;
plot_time_courses(datarun000_31,[3826  6874], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T2/', 'offt2extratc', gcf)


 %%
datarun000_31 = get_interspikeinterval(datarun000_31, offt2_31);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt2_31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, offt2_31(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T2/', 'offt2isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T2/', 'offt2isifull', gcf)

%%
wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); 
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,offt2_31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt2_31)
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
stairs([0 3 5 8 10],[0.5 0.45 0.4 0.45 0.45], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T2/', 'offt2psth2', gcf)
close;

spikesbytrials{1,1} = get_raster(datarun001_31.spikes{get_cell_indices(datarun001_31, offt2_31(1,1)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T2/', 'offt2psthcell421', gcf);
close;





%% OFF T4

off_otherother31 = [49,259,272,391,393,591,681,724,783,1368,1430,1501,1577,1731,1741,1772,1773,1816,1862,1954,2177,2178,2539,2597,2659,2867,2884,2899,2971,2973,3616,3905,3946,4036,4066,4204,4383,4518,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5926,5943,6033,6062,6091,6212,6213,6376,6483,6499,6541,6542,6796,6888,6901,6993,7157,7203,7324,7444,7502,7562,7593,7668];
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

% wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
% gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); %Might change with dataset
% [h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,off_otherother31), 0, '/0', wh, gr, 10, false);
% binSize = 0.1:0.1:10; %change depending on length of trial
% pulsePSTH = [];
% normvalpulsePSTH = [];
% pulsenormnormPSTH = [];
% maxpulse = [];
% maxpulsetime =[];
%  for a = 1:length(off_otherother31)
%  pulsePSTH(:,a) = sum(nhist{a,1})./50; %change depending on num of trials
% end
%  
% for i = 1:length(off_otherother31)
%  normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
% end
% normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
% pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;
% 
% [maxpulse maxpulsetime] = max(pulsePSTH);
% maxpulsetime = maxpulsetime*0.1;

 
[C ia ib] = intersect(offt4_31_init, off_otherother31);
vc = ones(length(off_otherother31),1);
vc(ib) = 2;


[COEFF1,SCORE1] = princomp(isinormnorm');
% [COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF2,SCORE2] = princomp(tcnormminn');

X = [];
X(:,1) = SCORE1(:,1);
X(:,2) = TCParams.minval;
X(:,3) = TCParams.mintim;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_otherother31, tcnormnorm,0,vc);
 offt4_31 = off_otherother31(idx==2);
 off_otherotherother31 = off_otherother31(idx==1);
ismember(offt4_31_init, off_otherother31(idx==2))
offt4_31 = [591,681,1501,1772,2178,2899,2973,3946,4036,6499,7444];
%off_otherotherother31 = [49,259,272,391,393,724,783,1368,1430,1577,1731,1741,1773,1816,1862,1954,2177,2539,2597,2659,2867,2884,2971,3616,3905,4066,4204,4383,4518,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5926,5943,6033,6062,6091,6212,6213,6376,6483,6541,6542,6796,6888,6901,6993,7157,7203,7324,7502,7562,7593,7668];
 
%% off t4 extra cells
plot_rf_summaries(datarun000_31, [offt4_31 7280 7201 3514 1216 96 ], 'coordinates', 'monitor');
plot_time_courses(datarun000_31,[offt4_31 7280 7201 3514 1216 96], 'all', true, 'bw', true);
%%
plot_rf_summaries(datarun000_31, [7280 7201 3514 1216 96 ], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T4/', 'offt4extrarf', gcf)

close all;
plot_time_courses(datarun000_31,[7280 7201 3514 1216 96], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T4/', 'offt4extratc', gcf)

%%
plot_rf_summaries(datarun000_31, [offt4_31], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T4/', 'offt4rf', gcf)

close all;
plot_time_courses(datarun000_31,[offt4_31], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T4/', 'offt4tc', gcf)

 %%
datarun000_31 = get_interspikeinterval(datarun000_31, offt4_31);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt4_31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, offt4_31(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T4/', 'offt4isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T4/', 'offt4isifull', gcf)

%%
wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); 
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,offt4_31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt4_31)
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
stairs([0 3 5 8 10],[0.5 0.45 0.4 0.45 0.45], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T4/', 'offt4psth2', gcf)
close;

spikesbytrials{1,1} = get_raster(datarun001_31.spikes{get_cell_indices(datarun001_31, offt4_31(1,2)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T4/', 'offt4psthcell681', gcf);
close;



 %% OFF T3
 
off_otherotherother31 = [49,259,272,391,393,724,783,1368,1430,1577,1731,1741,1773,1816,1862,1954,2177,2539,2597,2659,2867,2884,2971,3616,3905,4066,4204,4383,4518,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5926,5943,6033,6062,6091,6212,6213,6376,6483,6541,6542,6796,6888,6901,6993,7157,7203,7324,7502,7562,7593,7668];
 offt3_31_init =[1577,1954,2659,7324,6033,4383,2884];

[tc nontc] = get_time_courses_matrix(datarun000_31, off_otherotherother31); %or cellids
x = 1:1:30;
 normval = [];
 mx = [];
tcnormnorm = [];
tcnormmx = [];
for i = 1:length(off_otherotherother31) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);

tcnormnorm = tc./normval;
[mx mxt] = max(tc);
mx = repmat(mx, size(tc, 1), 1);
tcnormmx = tc./mx;

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

[TCParams] = time_course_parameters(tcfittednormnorm, 0);


vc = [];
[C ia ib] = intersect(offt3_31_init, off_otherotherother31);
vc = ones(length(off_otherotherother31),1);
vc(ib) = 2;


X = [];
X(:,1) = TCParams.maxmingrad;
X(:,2) = TCParams.maxtim;
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_otherotherother31, tcnormnorm,0, vc);
%off_otherotherother31(idx==2)
%ismember(off_otherotherother31(idx==2),offt3_31_init)
offt3_31 = off_otherotherother31(idx==2);
off_otherotherotherother31 = off_otherotherother31(idx==1);
offt3_31 =[1577,1954,2659,2884,4383,6033,7324]; 
off_otherotherotherother31 = [49,259,272,391,393,724,783,1368,1430,1731,1741,1773,1816,1862,2177,2539,2597,2867,2971,3616,3905,4066,4204,4518,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5926,5943,6062,6091,6212,6213,6376,6483,6541,6542,6796,6888,6901,6993,7157,7203,7502,7562,7593,7668];

%% off t3 extra cells
plot_rf_summaries(datarun000_31, [offt3_31  6766  5657  3797  2929], 'coordinates', 'monitor');
plot_time_courses(datarun000_31,[offt3_31 6766 5657 3797 2929], 'all', true, 'bw', true);

%%
plot_rf_summaries(datarun000_31, [6766  5657  3797  2929], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/', 'offt3rfextra', gcf)

%%
close all;
plot_time_courses(datarun000_31,[6766  5657  3797  2929], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/', 'offt3tcextra', gcf)


%% OFF T5

off_otherotherotherother31 = [49,259,272,391,393,724,783,1368,1430,1731,1741,1773,1816,1862,2177,2539,2597,2867,2971,3616,3905,4066,4204,4518,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5926,5943,6062,6091,6212,6213,6376,6483,6541,6542,6796,6888,6901,6993,7157,7203,7502,7562,7593,7668];
%offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];
%plot_time_courses(datarun000_31, [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157], 'all', true, 'bw', true);

offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376]
%offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];

[tc nontc] = get_time_courses_matrix(datarun000_31, off_otherotherotherother31); %or cellids
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
for i = 1:length(off_otherotherotherother31) %or nonds
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

temp_tcs = get_time_courses_matrix(datarun000_31, off_otherotherotherother31);
tc_fit = [];
final_params  =[];
for i = 1:length(off_otherotherotherother31)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(off_otherotherotherother31)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(off_otherotherotherother31) %or nonds
 normval(1, i) = norm( tcfitted(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tcfitted, 1), 1);
tcfittednormnorm = tcfitted./normval;   

[TCParams] = time_course_parameters(tcfittednormnorm, 1);

radius = [];
radius = get_rf_fit_radius(datarun000_31, off_otherotherotherother31);


datarun000_31 = get_interspikeinterval(datarun000_31, off_otherotherotherother31);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherotherotherother31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, off_otherotherotherother31(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;



[C ia ib] = intersect(offt5_31_init, off_otherotherotherother31);
vc = ones(length(off_otherotherotherother31),1);
vc(ib) = 2;


[COEFF,SCORE] = princomp(tcnormnorm');

X = [];
X(:,1) = TCParams.maxtim;
X(:,2) = TCParams.mintim;
X(:,3) = TCParams.dot;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_otherotherotherother31, tcnormnorm,0, vc);
offt5_31 = off_otherotherotherother31(idx==2);
off_otherotherotherotherother31 = off_otherotherotherother31(idx==1);%%
offt5_31 = [724,1731,1773,1816,2867,4204,5117,5431,6376,6483,6542,6796];
off_otherotherotherotherother31 = [49,259,272,391,393,783,1368,1430,1741,1862,2177,2539,2597,2971,3616,3905,4066,4518,4681,4727,4730,4771,4786,5072,5089,5495,5506,5596,5841,5926,5943,6062,6091,6212,6213,6541,6888,6901,6993,7157,7203,7502,7562,7593,7668];
%% off t5 extra cells
plot_rf_summaries(datarun000_31, [offt5_31   1009 3288 3361 3678 5326 5959 ], 'coordinates', 'monitor');
plot_time_courses(datarun000_31,[offt5_31 1009 3288 3361 3678 5326 5959 ], 'all', true, 'bw', true);  

%%
plot_rf_summaries(datarun000_31, [offt5_31], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/', 'offt5rf', gcf)

%%
close all;
plot_time_courses(datarun000_31,[offt5_31], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/', 'offt5tc', gcf)




%% off t6
off_otherotherotherotherother31 = [49,259,272,391,393,783,1368,1430,1741,1862,2177,2539,2597,2971,3616,3905,4066,4518,4681,4727,4730,4771,4786,5072,5089,5495,5506,5596,5841,5926,5943,6062,6091,6212,6213,6541,6888,6901,6993,7157,7203,7502,7562,7593,7668];
offt6_31_init = [ 4786 5596 6888];%1501 2177];

datarun000_31 = get_interspikeinterval(datarun000_31, off_otherotherotherotherother31 );
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(off_otherotherotherotherother31 ) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, off_otherotherotherotherother31 (1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;


[tc nontc] = get_time_courses_matrix(datarun000_31, off_otherotherotherotherother31 ); %or cellids
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
for i = 1:length(off_otherotherotherotherother31 ) %or nonds
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
cellind = get_cell_indices(datarun000_31, off_otherotherotherotherother31 );
stamat = cell(length(cellind),1);
 
for i = 1:length(cellind)
    B = [];
    B = datarun000_31.stas.rfs{cellind(i), 1}(:)';
    meanpix(i) = mean(B);
    rstd(i) = robust_std(B, [1]);
    stamat{i,1} = zeros(size(datarun000_31.stas.rfs{cellind(i),1},1),size(datarun000_31.stas.rfs{cellind(i),1},2));
    for j = 1:size(datarun000_31.stas.rfs{cellind(i),1},1)
            for k = 1:size(datarun000_31.stas.rfs{cellind(i),1},2)
                if (datarun000_31.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
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



 num_rgcs = length(off_otherotherotherotherother31 );

% initialize some variables for the look
rf_areas = zeros(num_rgcs,1);
abs_mean_pixel_val = zeros(num_rgcs,1);
snrs = zeros(num_rgcs, 1);
contrast_index = zeros(num_rgcs, 1);
[Y, X] = meshgrid(1:1:40, 1:1:80);


for rgc = 1:num_rgcs
    
    temp_index= get_cell_indices(datarun000_31, off_otherotherotherotherother31 (rgc));
    
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
        temp_rf = get_rf(datarun000_31, off_otherotherotherotherother31 (rgc));
 
      
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


[C ia ib] = intersect(offt6_31_init, off_otherotherotherotherother31 );
vc = ones(length(off_otherotherotherotherother31 ),1);
vc(ib) = 2;


X = [];
X(:,1) = SCORE(:,1);
X(:,2) = 1 - v(2,:)./v(1,:);
X(:,3) = contrast_index;
[idx obj] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_otherotherotherotherother31 , tcnormnorm,0, vc);
offt6_31  = off_otherotherotherotherother31 (idx==2);
off_otherotherotherotherotherother31  = off_otherotherotherotherother31 (idx==1);%%
offt6_31 = offt6_31_init;
%%off_otherotherotherotherotherother31 = off_otherotherotherotherother31(~ismember( off_otherotherotherotherother31, offt6_31));

%% off t6 extra cells
plot_rf_summaries(datarun000_31, [offt6_31  1262 1592 ], 'coordinates', 'monitor');
plot_time_courses(datarun000_31,[offt6_31  1262 1592 ], 'all', true, 'bw', true);  
%%
plot_rf_summaries(datarun000_31, [on_otherotherother31 ], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/', 'onunclassfrf', gcf)

close all;
plot_time_courses(datarun000_31,[on_otherotherother31], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/', 'onuncassftc', gcf)

%%
plot_rf_summaries(datarun000_31, [offt6_31], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T6/', 'offt6rf', gcf)

close all;
plot_time_courses(datarun000_31,[offt6_31], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T6/', 'offt6tc', gcf)



 %%
datarun000_31 = get_interspikeinterval(datarun000_31, offt6_31);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
maxisi = [];
isimax = [];
for i = 1:length(offt6_31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, offt6_31(1,i)), 1}.probabilities;
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
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T6/', 'offt6isihalf', gcf)
figure();
shadedErrorBar(x2,mean(isinormnorm(:, :)'),std(isinormnorm(:, :)'),'k');
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T6/', 'offt6isifull', gcf)

%%
wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); 
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,offt6_31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; 
psthnorm = [];
psthind = [];
psthindnorm = [];
for i = 1:length(offt6_31)
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
stairs([0 3 5 8 10],[0.5 0.45 0.4 0.45 0.45], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T6/', 'offt6psth2', gcf)
close;

spikesbytrials{1,1} = get_raster(datarun001_31.spikes{get_cell_indices(datarun001_31, offt6_31(1,1)), 1}, wh(1:2), 'tic_color', [0 0 0], 'axis_range', [0 10 0 4]);
hold on;
stairs([0 3 5 8 10],[w g b g g], 'Color', 'k', 'LineWidth',1);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/T6/', 'offt6psthcell4786', gcf);
close;


%%
plot_rf_summaries(datarun000_31, [offt5_31], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/', 'offt5rf', gcf)

%%
close all;
plot_time_courses(datarun000_31,[offt5_31], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/OFF/', 'offt5tc', gcf)



%% nov 26th finding extra cells for off t4, t3, t5, t6
temp_tcs = get_time_courses_matrix(datarun000_31, datarun000_31.cell_ids);
tc_fit = [];
final_params  =[];
for i = 1:length(datarun000_31.cell_ids)
[tc_fit(i,:), final_params(i,:)] = fit_time_course(temp_tcs(:,i), 'verbose', false);
end
tcfitted = [];
for i = 1:length(datarun000_31.cell_ids)
    params = final_params(i,:);
    t_points = (1:0.1:params(6))-1;
    t_filter_one = params(1) .* (t_points ./ params(3)).^params(5) .* exp(-params(5)*((t_points ./ params(3)) - 1));
    t_filter_two = params(2) .* (t_points ./ params(4)).^params(5) .* exp(-params(5)*((t_points ./ params(4)) - 1));
    tcbef = t_filter_one + t_filter_two;
    tcfitted(:,i) = fliplr(tcbef);
end
normval = [];
tcfittednormnorm = [];
for i = 1:length(datarun000_31.cell_ids) %or nonds
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


[tc nontc] = get_time_courses_matrix(datarun000_31, datarun000_31.cell_ids); %or cellids
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(datarun000_31.cell_ids) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

on_init31 = [4,78,108,286,362,543,618,647,781,785,903,917,962,991,1067,1081,1203,1248,1306,1381,1471,1595,1606,1712,1756,1878,1952,2086,2087,2146,2371,2419,2494,2536,2582,2687,2719,2747,2851,2856,2896,2897,3076,3215,3244,3289,3843,3887,3931,4038,4113,4156,4172,4186,4248,4277,4384,4685,4697,4713,4714,4876,4892,4983,5028,5042,5132,5240,5267,5297,5328,5401,5552,5642,5746,5765,5792,5794,5866,5897,5956,5961,6002,6031,6211,6289,6497,6590,6646,6662,6692,6752,6753,6784,6812,6962,7067,7098,7114,7336,7442,7472,7640,7669];
[C ia ib] = intersect(on_init31, datarun000_31.cell_ids);
vc = ones(length(datarun000_31.cell_ids),1);
vc(ib) = 2; %initializing on cells to cluster 2, everything else cluster 1

X = [];
X(:,1) = t_points(minnt);
X(:,2) = extrval;
[idx] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun000_31, datarun000_31.cell_ids, tcnormnorm,0, vc);
on_allsnr31_all = datarun000_31.cell_ids(idx ==2);
off_allsnr31_all = datarun000_31.cell_ids(idx ==1);



c = get_cell_indices(datarun000_31, on_allsnr31_all);
snronall = [];
for i = 1:length(c)
    r1 = sort(datarun000_31.stas.rfs{c(1,i),1}(:)', 'descend');
    snronall(1,i) = mean(r1(1:4))./std(r1);
end
on_31_all = on_allsnr31_all(snronall > (mean(snronall) - 2.5*std(snronall)));

% hax=axes; 
% hold on;
% hist(snronall)
% SP= mean(snronall) - 2.5*std(snronall); %your point goes here 
% line([SP SP],get(hax,'YLim'),'Color',[1 0 0])
% title('On cutoff')

c = get_cell_indices(datarun000_31, off_allsnr31_all);
snroffall = [];
for i = 1:length(c)
    r1 = sort(datarun000_31.stas.rfs{c(1,i),1}(:)', 'descend');
    snroffall(1,i) = mean(r1(1:4))./std(r1);
end
off_31_all = off_allsnr31_all(snroffall > (mean(snroffall) - 2.5*std(snroffall)));

%% end of nov 26th finding extra cells for off t4, t3, t5, t6

offt1_31_all = [offt1_31 6318];
offt2_31_all = [offt2_31 3826 6874 ];
offt4_31_all = [offt4_31 7280 7201 3514 1216 96 ];
offt3_31_all = [offt3_31 6766 5657 3797 2929];
offt5_31_all = [offt5_31 1009 3288 3361 3678 5326 5959 ];
offt6_31_all = [offt6_31];
off_otherotherotherotherotherother31_all = off_31_all(~ismember(off_31_all, [offt1_31_all offt2_31_all offt4_31_all offt3_31_all offt5_31_all offt6_31_all]));




%% OTHER STUFF



addpath('C:\Users\Sneha\Dropbox\Lab\Development\matlab-standard\code\lab')
 addpath('C:\Users\Sneha\Dropbox\Lab\Development\matlab-standard\code\lab\utilities')
 addpath('C:\Users\Sneha\Dropbox\Lab\Development\matlab-standard\code\lab\snlepaths');
 
labPath = 'C:\Users\Sneha\Dropbox\Lab';
matlabPath = [labPath '\Development\matlab-standard\code'];
addpath([matlabPath '\lab\utilities'])
addpath(genpath2(matlabPath, {'.svn'}));
visionPath = [labPath '\Applications\Vision.app\Contents\Resources\Java\Vision.jar'];
javaaddpath(visionPath);
set(0, 'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Helvetica')


 datarun001_31.names.rrs_neurons_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-31-0\data000-1800-7200\data001-map\data001-map.neurons';
 datarun001_31 = load_neurons(datarun001_31);
 datarun001_31.names.rrs_params_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-31-0\data000-1800-7200\data001-map\data001-map.params';
datarun001_31 = load_params(datarun001_31);

 
 datarun002_31.names.rrs_neurons_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-31-0\data000-1800-7200\data002-map\data002-map.neurons';
  datarun002_31 = load_neurons(datarun002_31);
  datarun002_31.names.rrs_params_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-31-0\data000-1800-7200\data002-map\data002-map.params';
 datarun002_31 = load_params( datarun002_31);
  datarun002_31.names.stimulus_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-31-0\data000-1800-7200\stimuli\s02';
    datarun002_31 = load_stim(datarun002_31);
    
    
 datarun000_31.names.rrs_neurons_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-31-0\data000-1800-7200\data000-map\data000-map.neurons';
 datarun000_31 = load_neurons(datarun000_31);
 datarun000_31.names.rrs_sta_path = 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-31-0\data000-1800-7200\data000-map\data000-map.sta';
 datarun000_31 = load_sta(datarun000_31, 'load_sta', 'all');
    marks_params.thresh = 4.0;
    datarun000_31 = get_sta_summaries(datarun000_31, 'all', 'marks_params', marks_params);
    %params and get sta fits and ei dont work
 
  
 
   datarun000_31 = get_interspikeinterval(datarun000_31, on_otherother31);

offt5_31 = [1731 1816 2867 4204 5117 5431 6542 6796 724 1773  6376 6993 1368 6483 7203 7157];
allplotsinone(datarun000_31, datarun001_31, datarun002_31, [1606 3931 4113 4156 6211], 0,0,0, 0,0,0, 0, 0, 'C:\Users\Sneha\Dropbox\Retina lab work\Analysis\2012-10-31-0\', wh, gr)

 ont3_31 = [1067 1306 1595 1878 2851 5328 6692 6753 2494 781 785 2747 5240];
plot_time_courses(datarun000_31,[1067 1306 1595 1878 2851 5328 6692 6753  5240], 'all' , true, 'bw', true, 'normalize', true);


datarun000_31 = get_interspikeinterval(datarun000_31, on_otherother31);
x2 = 0:0.001:0.1; 
isi = [];
normvalisi = [];
for i = 1:length(on_otherother31) %or nonds
 isi(:,i) = datarun000_31.interspikeinterval{get_cell_indices(datarun000_31, on_otherother31(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

plot(x2, isinormnorm(:,~ismember(on_otherother31,ont3_31)), 'b')
hold on;
plot(x2, isinormnorm(:,ismember(on_otherother31,ont3_31)), 'r')
title('ISI-norm - T3 in red')







[tc nontc] = get_time_courses_matrix(datarun000_31, on_otherother31); %or cellids
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

for i = 1:length(on_otherother31) %or nonds
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

plot(x, tcnormnorm(:,~ismember(on_otherother31,ont3_31)), 'b')
hold on;
plot(x, tcnormnorm(:,ismember(on_otherother31,ont3_31)), 'r')
title('TC-norm - T5 in red')


wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); %For 2012-10-10-1 dataset, that is how the triggers are arranges
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,on_otherother31), 0, '/0', wh, gr, 10, false);
binSize = 0.1:0.1:10; %change depending on length of trial
pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(on_otherother31)
 pulsePSTH(:,a) = sum(nhist{a,1})./50; %change depending on num of trials
end
 
plot(binSize, pulsePSTH(:,~ismember(on_otherother31,ont3_31)), 'b')
hold on;
plot(binSize, pulsePSTH(:,ismember(on_otherother31,ont3_31)), 'r')
title('Pulse PSTH - T4 in red')

for i = 1:length(on_otherother31)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
end
normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;

[maxpulse maxpulsetime] = max(pulsePSTH);
maxpulsetime = maxpulsetime*0.1;

plot(binSize, pulsenormnormPSTH(:,~ismember(on_otherother31,ont3_31)), 'b')
hold on;
plot(binSize, pulsenormnormPSTH(:,ismember(on_otherother31,ont3_31)), 'r')
 

 
 
 [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31,on_otherother31,0);
 
 NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 64),:));
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
plot(A32(~ismember(on_otherother31,ont3_31)), A256(~ismember(on_otherother31,ont3_31)), 'ob')
hold on;
plot(A32(ismember(on_otherother31,ont3_31)), A256(ismember(on_otherother31,ont3_31)), 'or')


plot(binSize, pulsePSTH(:,~ismember(on_otherother31,ont3_31)), 'b')
hold on;
plot(binSize, pulsePSTH(:,ismember(on_otherother31,ont3_31)), 'r')
title('Pulse PSTH - T4 in red')




 %ont3_31 = [1067 1306 1595 1878 2851 5328 6692 6753 2494 781 785 2747 5240];

 ont3_31 = [1067 1306 1595 1878 2851 5328 6692 6753 2494 2747 5240];
ont3_31 = [1067 1306 1595 1878 2851 5328 6692 6753 ];

vc = [];
[C ia ib] = intersect(ont3_31, on_otherother31);
vc = ones(length(on_otherother31),1);
vc(ib) = 2;
[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsenormnormPSTH');
[COEFF2,SCORE2] = princomp(tcnormmx');

X = [];
X(:,1) = SCORE2(:,1);
X(:,2) = SCORE2(:,2);
X(:,3) = SCORE2(:,3);
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, on_otherother31, tcnormnorm,0,vc);
 on_otherother31(idx==2)
 ismember(ont3_31, on_otherother31(idx==2))



%%
%% 2012-10-10-1
0 right  706 1023 1666 3574 4323 4366
90 up   1879 3470 3828 4549 5465 5627 7381
180 left  1369 1595 1681 1712 3306 3769 4577 4833 5375 6186 6633 7114
270 down   395 695 995 1265 1416 2555 2631 2975 3931 5419 6110 6138 6544 6754


 
%% 2012-10-15-0
0 right 257  1097 1683 1895 2898 3512 3736 4069  4157 4846 5632 6321 6751 6797
90 up  708 3815 6229
180 left 301 2449 3934 4173 4789 5719 6332 7308
270 down 333 438 467 1595 1685  2042 3636 3842 4353 4731 4985 5150 5702
 
%% 2012-01-31-0
0 right 92 484 619 1487 1683 1996 2074 2613 3019 3077 6242 7503
90 up 316 1037 2118 3125 4548 5148 6602 6781
180 left 559 1442 2433 3260 4159 5000 6290 6708 7518 7159 7083
270 down 2942 4625 5210
 
 %% Plot EIs
close all;
zero_degs = [92 484 619 1487 1683 1996 2074 2613 3019 3077 6242 7503];
color_palette = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 0.5 0.5; 0 0 0; 1 0.5 0; 0.5 1 0; 0 0.5 1; 0.25 0.5 1; 1 0.25 0.5; 0.5 1 0.25; 0.5 0.25 1; 1 0.5 0.25];
 eiplots(datarun000_31, zero_degs, 'color_palette', color_palette, 'cutoff', 0.03, 'alpha', 1, 'scale',1,'exclude_axon', false, 'plotei', false, 'pause', false,'plot_sourceneighb', true, 'neighbor_size', 12, 'save', false, 'pathname', '/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/DS/', 'filename', 'test');

 figure();
  eiplots(datarun000_31, zero_degs, 'color_palette', color_palette, 'plotei', true,'cutoff', 0.05,'exclude_axon',true, 'threshold', 4, 'alpha', 0.5, 'scale',1, 'pause', false, 'fit_Gaussian', true, 'contour_type','ezcontour', 'contourgrid', 250, 'save', false, 'pathname', '/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/DS/', 'filename', 'test');

   eiplots(datarun000_31, zero_degs, 'color_palette', color_palette, 'plotei', true,'cutoff', 0.05,'exclude_axon',true, 'threshold', 4, 'alpha', 0.5, 'scale',1, 'pause', false, 'fit_Gaussian', true, 'contour_type','simplecontour', 'contourgrid', 250, 'save', false);
 
%%
%% 2012-10-10-1 EI DG + WN run
0 right  706 1023 1666 3574 4323 4366
90 up   1879 3470 3828 4549 5465 5627 7381
180 left  1369 1595 1681 1712 3306 3769 4577 4833 5375 6186 6633 7114
270 down   395 695 995 1265 1416 2555 2631 2975 3931 5419 6110 6138 6544 6754


close all;
zero_degs = [395 695 995 1265 1416 2555 2631 2975 3931 5419 6110 6138 6544 6754];
color_palette = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 0.5 0.5; 0 0 0; 1 0.5 0; 0.5 1 0; 0 0.5 1; 0.25 0.5 1; 1 0.25 0.5; 0.5 1 0.25; 0.5 0.25 1; 1 0.5 0.25];
 eiplots(datarun002_10, zero_degs, 'color_palette', color_palette, 'cutoff', 0.05, 'alpha', 0.5, 'scale',1,'exclude_axon', true, 'threshold', 4, 'plotei', true, 'save', true, 'pathname', '/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/DS/No axons/', 'filename', 'EI-DG-270-scale1co005a05at4');

 
 
 
  eiplots(datarun000_10, 6186, 'color_palette', color_palette, 'cutoff', 0.05, 'alpha', 0.5, 'scale',1,'exclude_axon', true, 'threshold', 7, 'plotei', true, 'save', false, 'pathname', '/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/DS/No axons/', 'filename', 'EI-270-scale1co005a05at4');

  
  [datarun002_31_nc] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-31-0/data000-1800-7200/', 'data002-map/data002-map', 1,'stimuli/s02',0);

  


  

43 ds cells

ds_31_wn = [92 484 619 1487 1683 1996 2074 2613 3019 3077 6242 7503 316 1037 2118 3125 4548 5148 6602 6781 559 1442 2433 3260 4159 5000 6290 6708 7518 7159 7083 2942 4625 5210];
  [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, datarun002_31.cell_ids, 0);
[mag dsindex magmax magave angle rho theta num U V] = dscellanalysis(NumSpikesCell, StimComb);





X = [];
N = [];
p = [];
X(:,1) = log(mag{1,1})';
X(:,2) = log(mag{2,1})';
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun002_31, datarun002_31.cell_ids, 0,0, 0);

plot(mag{1,1}, mag{2,1}, 'o');
ds_31_dg = datarun002_31.cell_ids(idx==1);
ismember(ds_31_wn, datarun002_31.cell_ids(idx==1));


[cell_list_map, failed_cells] = map_ei(datarun002_31,datarun002_31);



ismember(,datarun002_31_nc.cell_ids(idx==1));
ismember([92 484 619 1487 1683 1996 2074 2613 3019 3077 6242 7503 316 1037 2118 3125 4548 5148 6602 6781 559 1442 2433 3260 4159 5000 6290 6708 7518 7159 7083 2942 4625 5210], datarun002_31_nc.cell_ids);

color_palette = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 0.5 0.5; 0 0 0; 1 0.5 0; 0.5 1 0; 0 0.5 1; 0.25 0.5 1; 1 0.25 0.5; 0.5 1 0.25; 0.5 0.25 1; 1 0.5 0.25];
color_palette = rand(length(ds_31_dg),3);
eiplots(datarun002_31, ds_31_dg, 'pause', true, 'plotei', true, 'color_palette', color_palette);
temp_indices = get_cell_indices(datarun002_31, ds_31_dg);


for j = 1:length(temp_indices)
for i = 1:length(temp_indices)
    eico = max(max(corr(datarun002_31.ei.eis{temp_indices(j),1},datarun002_31.ei.eis{temp_indices(i),1})));
    if(eico >= 0.95)
        if(i ~= j)
            i
            j
        eico
        pause;
        end
    end
end
end


max(max(corr(datarun002_31.ei.eis{temp_indices(1),1},datarun002_31.ei.eis{49,1})))




i is cell index

 mat = zeros(size(master_datarun.ei.eis{master_indices(1)},1),length(master_indices)+length(slave_indices)) %512 by 940+940 (1880)
 for i = 1: all indices
st = max(datarun.ei.eis{i},[],2);
temp = find(st < 5)
st(temp)=0;
if length(find(st > 5)) < 10
            exclude_list = [exclude_list i];
end
        mat(:,i) = st;
 end
 
 
  for i = 1:length(slave_indices)
        st = max(slave_datarun.ei.eis{slave_indices(i)},[],2);
        st(find(st < electrode_threshold)) = 0;
        mat(:,length(master_indices)+i) = st;
  end
    
  
  temp = corrcoef(mat);
    corr = temp(1:length(master_indices),length(master_indices)+1:end);
    corr(isnan(corr)) = 0;

 %%
 [datarun002_31_nc] = load_dsdata('/Analysis/sravi/Rat/WildType/2012-10-31-0/data000-1800-7200/', 'data002-map/data002-map', 1,'stimuli/s02',0);
sc = datarun002_31_nc.cell_ids(~ismember(1:1:length(datarun002_31_nc.cell_ids),(get_cell_indices(datarun002_31_nc,ds_31_wn(ismember(ds_31_wn, datarun002_31_nc.cell_ids))))));
f = ds_31_wn(~ismember(ds_31_wn, datarun002_31_nc.cell_ids)); %cells not in neuron cleaned DG run - only in WN
c = map_ei(datarun000_31, datarun002_31_nc, 'master_cell_type', f, 'slave_cell_type', sc, 'corr_threshold', 0.5); %map EIs
c(get_cell_indices(datarun000_31, f)) %cells in DG run
s = ds_31_wn(ismember(ds_31_wn, datarun002_31_nc.cell_ids));
ds_31_dg_in = [s cell2mat(c(get_cell_indices(datarun000_31, f)))];

  
map_ei(datarun002_31, datarun002_31, 'master_cell_type', 7518, 'slave_cell_type',ds_31_dg(69:end));

[c f] = map_ei(datarun000_31, datarun002_31_nc, 'master_cell_type', [ds_31_wn(~ismember(ds_31_wn, datarun002_31_nc.cell_ids))], 'slave_cell_type', '', 'electrode_threshold',4);

%%

  [NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31_nc, datarun002_31_nc.cell_ids, 0);
[mag dsindex magmax magave angle rho theta num U V] = dscellanalysis(NumSpikesCell, StimComb);

[C ia ib] = intersect(ds_31_dg_in, datarun002_31_nc.cell_ids);
vc = ones(length(datarun002_31_nc.cell_ids),1);
vc(ib) = 2;

X = [];
N = [];
p = [];
X(:,1) = log(mag{1,1})';
X(:,2) = log(mag{2,1})';
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, datarun002_31, datarun002_31.cell_ids, 0,0, vc);
ds31dg = datarun002_31_nc.cell_ids(idx==2);
moreds = ds31dg(~ismember(datarun002_31_nc.cell_ids(idx==2),ds_31_dg_in));

for i = 1:length(moreds)
 polar_plots_one(rho,theta,U,V,num, get_cell_indices(datarun002_31_nc, moreds(i)))
 pause;
end

%







43 ds cells

ds_31_wn = [92 484 619 1487 1683 1996 2074 2613 3019 3077 6242 7503 316 1037 2118 3125 4548 5148 6602 6781 559 1442 2433 3260 4159 5000 6290 6708 7518 7159 7083 2942 4625 5210];

0 right 92 484 619 1487 1683 1996 2074 2613 3019 3077 6242 7503
90 up 316 1037 2118 3125 4548 5148 6602 6781
180 left 559 1442 2433 3260 4159 5000 6290 6708 7518 7159 7083
270 down 2942 4625 5210

close all;
zero_degs = [632 1442 2433 4159 5000 6290 6708 7518 7159 7083];
color_palette = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 0.5 0.5; 0 0 0; 1 0.5 0; 0.5 1 0; 0 0.5 1; 0.25 0.5 1; 1 0.25 0.5; 0.5 1 0.25; 0.5 0.25 1; 1 0.5 0.25; 0 0.75 0.3; 0.3 0.75 0; 0 0.3 0.75; 0.3 0 0.75; 0.75 0 0.3; 0.75 0.3 0];
 eiplots(datarun002_31_nc, zero_degs, 'color_palette', color_palette, 'cutoff', 0.05, 'alpha', 0.5, 'scale',1,'exclude_axon', true, 'threshold', 4, 'plotei', true, 'pause', false, 'save', false, 'pathname', '/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/DS/', 'filename', 'test');
hold on;
%pause;
 eiplots(datarun002_31_nc, [1021 2855 3335 3363 6469], 'color_palette', color_palette, 'cutoff', 0.05, 'alpha', 0.5, 'scale',1,'exclude_axon', true, 'threshold', 4, 'plotei', true, 'pause', false, 'save', true, 'pathname', '/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/DS/No axons/', 'filename', 'EI-DG-180-NCv1-scale1co005a05at4');







plot(mag{1,1}, mag{2,1}, 'o');
ds_31_dg = datarun002_31.cell_ids(idx==1);
ismember(ds_31_wn, datarun002_31.cell_ids(idx==1));


[cell_list_map, failed_cells] = map_ei(datarun002_31,datarun002_31);



ismember(,datarun002_31_nc.cell_ids(idx==1));
ismember([92 484 619 1487 1683 1996 2074 2613 3019 3077 6242 7503 316 1037 2118 3125 4548 5148 6602 6781 559 1442 2433 3260 4159 5000 6290 6708 7518 7159 7083 2942 4625 5210], datarun002_31_nc.cell_ids);

color_palette = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 0.5 0.5; 0 0 0; 1 0.5 0; 0.5 1 0; 0 0.5 1; 0.25 0.5 1; 1 0.25 0.5; 0.5 1 0.25; 0.5 0.25 1; 1 0.5 0.25];
color_palette = rand(length(ds_31_dg),3);
eiplots(datarun002_31, ds_31_dg, 'pause', true, 'plotei', true, 'color_palette', color_palette);
temp_indices = get_cell_indices(datarun002_31, ds_31_dg);


for j = 1:length(temp_indices)
for i = 1:length(temp_indices)
    eico = max(max(corr(datarun002_31.ei.eis{temp_indices(j),1},datarun002_31.ei.eis{temp_indices(i),1})));
    if(eico >= 0.95)
        if(i ~= j)
            i
            j
        eico
        pause;
        end
    end
end
end


max(max(corr(datarun002_31.ei.eis{temp_indices(1),1},datarun002_31.ei.eis{49,1})))




i is cell index

 mat = zeros(size(master_datarun.ei.eis{master_indices(1)},1),length(master_indices)+length(slave_indices)) %512 by 940+940 (1880)
 for i = 1: all indices
st = max(datarun.ei.eis{i},[],2);
temp = find(st < 5)
st(temp)=0;
if length(find(st > 5)) < 10
            exclude_list = [exclude_list i];
end
        mat(:,i) = st;
 end
 
 
  for i = 1:length(slave_indices)
        st = max(slave_datarun.ei.eis{slave_indices(i)},[],2);
        st(find(st < electrode_threshold)) = 0;
        mat(:,length(master_indices)+i) = st;
  end
    
  
  temp = corrcoef(mat);
    corr = temp(1:length(master_indices),length(master_indices)+1:end);
    corr(isnan(corr)) = 0;
[t_corr, t_corr_index] = sort(corr(i,:));
    largest_corr=t_corr(end);
    max_corr_index=t_corr_index(end);
 
 %%
 %%

cellids = [ 1595 3842 2042 4985 333 438 467 1685 3636 4353 4731 5150 5702];

tc = [];
[tc nontc] = get_time_courses_matrix(datarun000_15, cellids); 
x = 1:1:30;
normval = [];
tcnormnorm = [];
for i = 1:length(cellids) %or nonds
 normval(1, i) = norm( tc(:,i)); %Calculate norm (magnitude) for all time courses
end 
normval = repmat(normval, size(tc, 1), 1);
tcnormnorm = tc./normval;

plot(1:1:30, tcnormnorm(:,1:4), 'b');
hold on;
plot(1:1:30, tcnormnorm(:,5:13), 'r');

plot(1:1:30, tcnormnorm(:,3:11), 'b');
plot(1:1:30, tcnormnorm(:,8), 'b');
plot(1:1:30, tcnormnorm(:,9:12), 'b');
plot(1:1:30, tcnormnorm(:,6:7), 'r');
plot(1:1:30, tcnormnorm(:,8), 'r');

title('2012-10-15 - 270 degrees - look diff - 333 438 467 1685 3636 4353 4731 5150 5702');


[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, cellids, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell, StimComb);
NS2 = [];
A32 = [];
A256 = [];
NS2 = NumSpikesCell';
A32 = sum(NS2(find(StimComb(:,2) == 64),:)); % CHANGE ACCORDING TO WHAT YOUR 2 TEMPORAL PERIODS ARE!
A256 = sum(NS2(find(StimComb(:,2) == 256),:));
close all;
cellids = [1879 3470 4549 5465 7381 3828 5627];

plot(A32(1:4), A256(1:4), 'ob');
hold on;
plot(A32(5:13), A256(5:13), 'or');

plot(A32(6:11), A256(6:11), 'ob');
plot(A32(8), A256(8), 'ob');
plot(A32(9:12), A256(9:12), 'ob');
plot(A32(6:7), A256(6:7), 'or');
plot(A32(8), A256(8), 'or');

title('2012-10-15 - 270 degrees - look diff - 333 438 467 1685 3636 4353 4731 5150 5702 - Drifting Grating Spike Rate');


cellids = [1879 3470 4549   5465 7381 3828 5627];

datarun000_15 = get_interspikeinterval(datarun000_15, cellids);
x2 = 0:0.001:0.1; 
%nonds - cells not ds, tc - their time courses
isi = [];
normvalisi = [];
isinormnorm = [];
for i = 1:length(cellids) %or nonds
 isi(:,i) = datarun000_15.interspikeinterval{get_cell_indices(datarun000_15, cellids(1,i)), 1}.probabilities;
 normvalisi(1, i) = norm( isi(:,i));
end 
normvalisi = repmat(normvalisi, size(isi, 1), 1);
isinormnorm = isi./normvalisi;

%cellids = [1879 3470 4549  3828 5627 5465 7381];

plot(0:0.001:0.1, isinormnorm(:,1:4), 'b');
hold on;
plot(0:0.001:0.1, isinormnorm(:,5:13), 'r');

plot(0:0.001:0.1, isinormnorm(:,8), 'b');
%plot(0:0.001:0.1, isinormnorm(:,9:12), 'b');
plot(0:0.001:0.1, isinormnorm(:,5), 'r');
plot(0:0.001:0.1, isinormnorm(:,6:7), 'r');
plot(0:0.001:0.1, isinormnorm(:,8), 'r');

title('2012-10-15 - 270 degrees - look diff - 333 438 467 1685 3636 4353 4731 5150 5702  - ISI');


cellids = [1879 3470 4549   3828 5627 5465 7381 ];

wh = datarun001_15.triggers(1:4:length(datarun001_15.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
gr = datarun001_15.triggers(2:4:length(datarun001_15.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_15, get_cell_indices(datarun001_15,cellids), 0, '/0', wh, gr, 10, false);
binSize = 0.1:0.1:10; %change depending on length of trial
pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(cellids)
 pulsePSTH(:,a) = sum(nhist{a,1})./25; %change depending on num of trials
end
 
for i = 1:length(cellids)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
end
normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;


plot(0.1:0.1:10, pulsenormnormPSTH(:,1:4), 'b');
hold on;
plot(0.1:0.1:10, pulsenormnormPSTH(:,5:13), 'r');
plot(0.1:0.1:10, pulsenormnormPSTH(:,8), 'b');
%plot(0.1:0.1:10, pulsenormnormPSTH(:,9:12), 'b');
plot(0.1:0.1:10, pulsenormnormPSTH(:,5), 'r');
plot(0.1:0.1:10, pulsenormnormPSTH(:,6:7), 'r');
plot(0.1:0.1:10, pulsenormnormPSTH(:,8), 'r');

title('2012-10-15 - 270 deg - look diff - 333 438 467 1685 3636 4353 4731 5150 5702 - ');
 

%%
[COEFF,SCORE] = princomp(tcnormauc');
subplot(3,1,1)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
legend('PC1', 'PC2', 'PC3');
title('TC-AUC');
[COEFF,SCORE] = princomp(tcnormnorm');
subplot(3,1,2)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('TC-NORM');
[COEFF,SCORE] = princomp(tcnormmx');
subplot(3,1,3)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('TC-MAX');





[COEFF,SCORE] = princomp(tcnormauc');
subplot(4,1,1)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
legend('PC1', 'PC2', 'PC3');
title('TC-AUC');
[COEFF,SCORE] = princomp(tcnormnorm');
subplot(4,1,2)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('TC-NORM');
[COEFF,SCORE] = princomp(tcnormmx');
subplot(4,1,3)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('TC-MAX');
[COEFF,SCORE] = princomp(isinormnorm');
subplot(4,1,4)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('ISI-NORM');










[COEFF,SCORE] = princomp(tcnormauc');
subplot(5,1,1)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
legend('PC1', 'PC2', 'PC3');
title('TC-AUC');
[COEFF,SCORE] = princomp(tcnormnorm');
subplot(5,1,2)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('TC-NORM');
[COEFF,SCORE] = princomp(tcnormmx');
subplot(5,1,3)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('TC-MAX');
[COEFF,SCORE] = princomp(tcnormminn');
subplot(5,1,4)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('TC-MIN');
[COEFF,SCORE] = princomp(isinormnorm');
subplot(5,1,5)
plot(COEFF(:,1),'r');
hold on;
plot(COEFF(:,2),'b');
plot(COEFF(:,3),'g');
title('ISI-NORM');

%% off t3 after t4

off_otherotherotherother31 = [49,259,272,391,393,724,783,1368,1430,1577,1731,1741,1773,1816,1862,1954,2177,2539,2597,2659,2867,2884,2971,3616,3905,4066,4204,4383,4518,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5926,5943,6033,6062,6091,6212,6213,6376,6483,6541,6542,6796,6888,6901,6993,7157,7203,7324,7502,7562,7593,7668];
offt3_31_init =[1577,1954,2659,7324,6033,4383,2884];


%% off t5 after t4
% off_otherotherotherother31 = [49,259,272,391,393,724,783,1368,1430,1577,1731,1741,1773,1816,1862,1954,2177,2539,2597,2659,2867,2884,2971,3616,3905,4066,4204,4383,4518,4681,4727,4730,4771,4786,5072,5089,5117,5431,5495,5506,5596,5841,5926,5943,6033,6062,6091,6212,6213,6376,6483,6541,6542,6796,6888,6901,6993,7157,7203,7324,7502,7562,7593,7668];
% offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];

%% oct 7th 2-14
offt3_10_init =[512,1368,2750,4188,4487,5732,6170,6511,6588];
plot_time_courses(datarun000_10, offt5_10_init, 'all', true, 'bw', true);

offt5_10_init = [800,2029,2597,3423,3617,4321,4519,5536,5614,6949,7353,3874];

[tc nontc] = get_time_courses_matrix(datarun000_15, [offt5_15_init]); %or cellids
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
    plot(1:30, tcnormnorm(:,i), 'r');
    pause;
    hold on;
end

%offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376]
offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];
offt3_31_init =[1577,1954,2659,7324,6033,4383,2884];

offt5_15_init = [1246,2253,3695,5116,6260,4771,6380,226];
offt3_15_init = [62,991,1156,4234,4278,4487,5733,6286,6931];


offt5_15_init = [1246,2253,3695,5116,6260]%4771,6380,226];

offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376]%6993,6483,7203,7157];


[tc nontc] = get_time_courses_matrix(datarun000_15, [offt5_15_init]); %or cellids
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

%% more
A = [];
A(:,1) = TCParams.maxmingrad;
A(:,2) = TCParams.maxval;
A(:,3) = TCParams.minval;
A(:,4) = TCParams.maxtim;
A(:,5) = TCParams.mintim;
A(:,6) = TCParams.dot;
A(:,7) = TCParams.zerocrossing;

B = {'maxmingrad'; 'maxval';'minval';'maxtim';'mintim';'dot';'zerocrossing'};

D = [];
D = [1 2 3; 1 2 4; 1 2 5; 1 2 6; 1 2 7; 1 3 4; 1 3 5; 1 3 6; 1 3 7; 1 4 5; 1 4 6; 1 4 7; 1 5 6; 1 5 7; 1 6 7; 2 3 4; 2 3 5; 2 3 6; 2 3 6; 2 4 5; 2 4 6; 2 4 7; 2 5 6; 2 5 7; 2 6 7; 3 4 5; 3 4 6; 3 4 7; 3 5 6; 3 5 7; 3 6 7; 4 5 6; 4 5 7; 4 6 7; 5 6 7;];


for i = 1:length(D)
    close all;
X = [];
X(:,1) = A(:,D(i,1));
X(:,2) = A(:,D(i,2));
X(:,3) = A(:,D(i,3));
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_otherotherotherother31, tcnormnorm,0,vc);
 %offt5_31 = off_otherotherotherother31(idx==2);
 off_otherotherotherother31(idx==2)
 ismember(off_otherotherotherother31(idx==2),offt5_31_init)
 sum(ismember(off_otherotherotherother31(idx==2),offt5_31_init))
 strcat(B(D(i,1)) , B(D(i,2))  ,B(D(i,3)))
 pause;
end

%%

addpath(path);
spatPer = find(unique(StimComb(:,1))==sp);
tempPer = find(unique(StimComb(:,2))==temp);


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
subpin = [12 5 4 3 10 17 18 19]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(3,7,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 'axis_range', [0 8 0 10], 'stop', 12, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,7,11)
%set(gca, 'FontName', 'Helvetica', 'FontSize', 11)
set(gca, 'FontSize', 9)
p = polar(T{tempPer,1}, R{tempPer,1});
h = findall(gca, 'type', 'line'); % find all lines in the current figure
h(h==p) = []; % remove the line you ploted from the list.
set(h, 'LineStyle', '--');
hold on;
h1 = compass(Unew(tempPer,1),Vnew(tempPer,1), 'r'); %Vector average plot
set(h1,'linewidth',2) 
hold off;

 
 
%%

to be continued oct 15th

offt5_10_init = [2029,2597,3423,3617,4321,4519,5536,6949,7353];%800 5614 3874 are not t5 - could possibly use pulse information
offt5_15_init = [1246,2253,3695,5116,6260];%4771 (pulses are good, tc not so much),6380,226 (pulses not good) 6886 (pulses good, tc not)
offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,6376]; %1773 and 1731 are duplicates; 6483 and 6993 are prob t5 but tc slightly off
%offt5_31_init = [1731,1816,2867,4204,5117,5431,6542,6796,724,1773,6376,6993,1368,6483,7203,7157];


offt3_10_init =[512,1368,2750,4188,4487,5732,6170,6511,6588]; %1895 is prob t3
 offt3_15_init = [62,991,1156,4234,4278,4487,5733,6286,6931]; %[6886 same tc but not same pulse and dg]
 offt3_31_init =[1577,1954,2659,7324,6033,4383,2884]; %5853 is prob t3 but mosaic viol - just need to correct fit [6290 also? 7157? 7203?]

 %%
 plot_time_courses(datarun000_31, [offt5_31_init],'all',true,'bw',true);
 
 %%
 
 plot_rf_summaries(datarun000_31, [ offt5_31_init 6483 6993],'label', true)
 
 %%
 
 offt5_10_init = [2029,2597,3423,3617,4321,4519,5536,6949,7353];%800 5614 3874 are not t5 - could possibly use pulse information

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

% radius = [];
% radius = get_rf_fit_radius(datarun000_10, off_otherotherother10);
% 
% datarun000_10 = get_interspikeinterval(datarun000_10, off_otherotherother10);
% x2 = 0:0.001:0.1; 
% isi = [];
% normvalisi = [];
% for i = 1:length(off_otherotherother10) %or nonds
%  isi(:,i) = datarun000_10.interspikeinterval{get_cell_indices(datarun000_10, off_otherotherother10(1,i)), 1}.probabilities;
%  normvalisi(1, i) = norm( isi(:,i));
% end 
% normvalisi = repmat(normvalisi, size(isi, 1), 1);
% isinormnorm = isi./normvalisi;

wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,off_otherotherother10), 0, '/0', wh, gr, 24, false);
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

vc = [];
[C ia ib] = intersect(offt5_10_init, off_otherotherother10);
vc = ones(length(off_otherotherother10),1);
vc(ib) = 2;


[COEFF,SCORE] = princomp(pulsePSTH');

X = [];
X(:,1) = TCParams.maxtim;
X(:,2) = TCParams.mintim;
X(:,3) = SCORE(:,3);
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherotherother10, tcnormnorm,0, vc);
off_otherotherother10(idx==2)
ismember(off_otherotherother10(idx==2), offt5_10_init)
%offt3_10 = off_otherotherother10(idx==2);
%off_otherotherotherother10 = off_otherotherother10(idx==1);

 %% Calculate OS index

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, off_otherotherother31, 0);
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
for i = 1:length(off_otherotherother31)
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
figure()
hist(maxAxFir./minAxFir)
xlabel('maximum firing / minimum firing');
title('S 64 T 256')

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
for i = 1:length(off_otherotherother31)
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
hist(maxAxFir./minAxFir)
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
hist(xxx.*yyy);
xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');

%%
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
off_otherotherother31(1,chos)

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


wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,off_otherotherother10), 0, '/0', wh, gr, 6, false,0.1);


pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(off_otherotherother10)
 pulsePSTH(:,a) = sum(nhist{a,1})./length(wh); %change depending on num of trials
end
 
for i = 1:length(off_otherotherother10)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
end
normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;

% [maxpulse maxpulsetime] = max(pulsePSTH);
% maxpulsetime = maxpulsetime*0.5;

vc = [];
[C ia ib] = intersect(offt5_10_init, off_otherotherother10);
vc = ones(length(off_otherotherother10),1);
vc(ib) = 2;


[COEFF,SCORE] = princomp(pulsenormnormPSTH');

X = [];
X(:,1) = TCParams.maxtim;
X(:,2) = TCParams.mintim;
X(:,3) = SCORE(:,1);
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherotherother10, tcnormnorm,0, vc);
off_otherotherother10(idx==2)
ismember(off_otherotherother10(idx==2), offt5_10_init)
sum(ans)
%offt3_10 = off_otherotherother10(idx==2);
%off_otherotherotherother10 = off_otherotherother10(idx==1);


%%
oct 16th onwards


[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, off_otherotherother15, 0);
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
for i = 1:length(off_otherotherother15)
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
figure()
hist(maxAxFir./minAxFir, 20)
xlabel('maximum firing / minimum firing');
title('S 64 T 256')

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
for i = 1:length(off_otherotherother15)
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
hist(maxAxFir./minAxFir,20)
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
hist(xxx.*yyy,100);
xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');

%%
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
off_otherotherother15(1,chos)

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


close all;
wh = datarun001_10.triggers(1:4:length(datarun001_10.triggers), 1); %For 2012-31-31-1 dataset, that is how the triggers are arranges - need to change with dataset
gr = datarun001_10.triggers(2:4:length(datarun001_10.triggers),1); %Might change with dataset
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_10, get_cell_indices(datarun001_10,off_otherotherother10), 0, '/0', wh, gr, 8, false,0.5);


pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(off_otherotherother10)
 pulsePSTH(:,a) = sum(nhist{a,1})./length(wh); %change depending on num of trials
end
 
for i = 1:length(off_otherotherother10)
 normvalpulsePSTH(1, i) = norm( pulsePSTH(:,i)); %Calculate norm (magnitude) for all time courses
end
normvalpulsePSTH = repmat(normvalpulsePSTH, size(pulsePSTH, 1), 1);
pulsenormnormPSTH = pulsePSTH./normvalpulsePSTH;

[maxpulse maxpulsetime] = max(pulsePSTH);
maxpulsetime = maxpulsetime*0.1;

vc = [];
[C ia ib] = intersect(offt5_10_init, off_otherotherother10);
vc = ones(length(off_otherotherother10),1);
vc(ib) = 2;
% 
% 
[COEFF,SCORE] = princomp(pulsenormnormPSTH');
% 


%maxtim mintim maxval minval dot
X = [];
X(:,1) = TCParams.maxval;
X(:,2) = TCParams.maxtim;
X(:,3) = SCORE(:,2);
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_10, off_otherotherother10, tcnormnorm,0, vc);
off_otherotherother10(idx==2)
ismember(off_otherotherother10(idx==2), offt5_10_init)
sum(ans)



%offt3_10 = off_otherotherother10(idx==2);
%off_otherotherotherother10 = off_otherotherother10(idx==1);


%% NEXT

 off_otherotherother15 = [226,407,649,843,872,1130,1411,1427,1564,1817,1921,2206,2236,2326,2343,2373,2521,2555,2732,2794,3091,3286,3857,4097,4235,4697,4732,4771,5148,5672,6139,6257,6376,6380,6886,6992,7517];
 off_otherotherother31 = [49,259,272,391,393,783,1368,1430,1741,1862,2177,2539,2597,2971,3616,3905,4066,4518,4681,4727,4730,4771,4786,5072,5089,5495,5506,5596,5841,5926,5943,6062,6091,6212,6213,6483,6541,6888,6901,6993,7157,7203,7502,7562,7593,7668];
off_otherotherother10 = [17,454,496,1231,1459,1486,1624,1701,1801,1803,1895,2057,2073,2462,2794,3016,3093,3363,4367,4384,4580,4996,5449,5746,5836,5959,6378,6483,6856,7336,7561];

offt6_10_init = [5836 6856 1486 496]; %17];
offt6_15_init = [1817 2343 3286 6992];
offt6_31_init = [ 4786 5596]; %1501 2177


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
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,off_otherotherother31), 0, '/0', wh, gr, 24, false,0.1);
binSize = 0.1:0.1:31; %change depending on length of trial
pulsePSTH = [];
normvalpulsePSTH = [];
pulsenormnormPSTH = [];
maxpulse = [];
maxpulsetime =[];
 for a = 1:length(off_otherotherother31)
 pulsePSTH(:,a) = sum(nhist{a,1})./31; %change depending on num of trials
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

%%
sumsumsp = [];
for i = 1: length(sumSpTrTrig)
    sumsumsp(i,:) = sum(sumSpTrTrig{i, 1});
end
  
%%

datarun000_31 = get_significant_stixels(datarun000_31, off_otherotherother31);
stavar = [];
for i = 1:length(off_otherotherother31)
        ii = get_cell_indices(datarun000_31, off_otherotherother31(1,i));
        stavar(i,1) = var(datarun000_31.stas.rfs{ii,1}(datarun000_31.stas.significant_stixels{ii, 1}));
        if(isnan(stavar(i,1)))
            stavar(i,1) = 0;
        end
end


%%
[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsenormnormPSTH');
%[COEFF2,SCORE2] = princomp(tcnormminn');

X = [];
X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
X(:,3) = stavar;
[idx N p] = clustering_analysis_plots(X, 0,1, 2, 1, 0, datarun000_31, off_otherotherother31, tcnormnorm,0,vc);
 %offt4_31 = off_otherotherother31(idx==2);
 %off_otherotherother31 = off_otherotherother31(idx==1);
 off_otherotherother31(idx==2)
ismember(offt6_31_init, off_otherotherother31(idx==2))

%% oct 27th
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
  


datarun000_31 = get_significant_stixels(datarun000_31, off_otherotherother31);
stavar = [];
for i = 1:length(off_otherotherother31)
        ii = get_cell_indices(datarun000_31, off_otherotherother31(1,i));
        stavar(i,1) = var(datarun000_31.stas.rfs{ii,1}(datarun000_31.stas.significant_stixels{ii, 1}));
        if(isnan(stavar(i,1)))
            stavar(i,1) = 0;
        end
end

%%
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, off_otherotherother31, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V spave] = dscellanalysis(NumSpikesCell, StimComb);

minangle = [];
maxangle = [];
minAxis = [];
maxAxis = [];
minAxFir = [];
maxAxFir = [];
minFir = cell(2,1);
maxFir = cell(2,1);
[CC,minAxis] = min(spave{2,1}');
for i = 1:length(off_otherotherother31)
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

minFir{2,1} = minAxFir;
maxFir{2,1} = maxAxFir;


minangle = [];
maxangle = [];
minAxis = [];
maxAxis = [];
minAxFir = [];
maxAxFir = [];
[CC,minAxis] = min(spave{1,1}')
for i = 1:length(off_otherotherother31)
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

% scatter(minAxFir,maxAxFir);
% hold on;
% xlabel('minimum firing');
% ylabel('maximum firing')
% title('S 64 T 32')
% %plot(0:1:250, 0:1:250);
% hold off;
% figure();
% hist(maxAxFir./minAxFir,20)
% xlabel('maximum firing / minimum firing');
% title('S 64 T 32')

minFir{1,1} = minAxFir;
maxFir{1,1} = maxAxFir;

%%
xxx = maxFir{1,1}./minFir{1,1};
yyy = maxFir{2,1}./minFir{2,1};
scatter(xxx,yyy);
% xlabel('maximum firing / minimum firing - T 32');
% ylabel('maximum firing / minimum firing - T 256')
% title('S 64 T 32-256')
% figure()
% hist(xxx.*yyy,100);
% xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');

%% ISI, TC PARAMS, TC PC, PULSE PC, PULSE PARAMS, STA
%isi 1 2 3  pulsesum pulsenorm 1 pulsenorm 2 pulsenorm 3 pulse 1 pulse 2 pulse 3  stavar
%X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));

[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsePSTH');
[COEFF2,SCORE2] = princomp(tcnormminn');

X = [];

X(:,1) = yyy;
X(:,2) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,3) = SCORE(:,1);
%scatter3 (X(:,1), X(:,2), X(:,3))



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
  


datarun000_10 = get_significant_stixels(datarun000_10, off_otherotherother10);
stavar = [];
for i = 1:length(off_otherotherother10)
        ii = get_cell_indices(datarun000_10, off_otherotherother10(1,i));
        stavar(i,1) = var(datarun000_10.stas.rfs{ii,1}(datarun000_10.stas.significant_stixels{ii, 1}));
        if(isnan(stavar(i,1)))
            stavar(i,1) = 0;
        end
end



[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsenormnormPSTH');
%[COEFF2,SCORE2] = princomp(tcnormminn');
%%
X = [];
X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
X(:,3) = stavar;
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
  


datarun000_15 = get_significant_stixels(datarun000_15, off_otherotherother15);
stavar = [];
for i = 1:length(off_otherotherother15)
        ii = get_cell_indices(datarun000_15, off_otherotherother15(1,i));
        stavar(i,1) = var(datarun000_15.stas.rfs{ii,1}(datarun000_15.stas.significant_stixels{ii, 1}));
        if(isnan(stavar(i,1)))
            stavar(i,1) = 0;
        end
end


[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_15, off_otherotherother15, 0);
[mag  dsindex  magmax  magave  angle  rho  theta  num  U  V spave] = dscellanalysis(NumSpikesCell, StimComb);

minangle = [];
maxangle = [];
minAxis = [];
maxAxis = [];
minAxFir = [];
maxAxFir = [];
minFir = cell(2,1);
maxFir = cell(2,1);
[CC,minAxis] = min(spave{2,1}');
for i = 1:length(off_otherotherother15)
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

minFir{2,1} = minAxFir;
maxFir{2,1} = maxAxFir;


minangle = [];
maxangle = [];
minAxis = [];
maxAxis = [];
minAxFir = [];
maxAxFir = [];
[CC,minAxis] = min(spave{1,1}')
for i = 1:length(off_otherotherother15)
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

% scatter(minAxFir,maxAxFir);
% hold on;
% xlabel('minimum firing');
% ylabel('maximum firing')
% title('S 64 T 32')
% %plot(0:1:250, 0:1:250);
% hold off;
% figure();
% hist(maxAxFir./minAxFir,20)
% xlabel('maximum firing / minimum firing');
% title('S 64 T 32')

minFir{1,1} = minAxFir;
maxFir{1,1} = maxAxFir;

%%
xxx = maxFir{1,1}./minFir{1,1};
yyy = maxFir{2,1}./minFir{2,1};
scatter(xxx,yyy);
% xlabel('maximum firing / minimum firing - T 32');
% ylabel('maximum firing / minimum firing - T 256')
% title('S 64 T 32-256')
% figure()
% hist(xxx.*yyy,100);
% xlabel('maximum firing / minimum firing - T 32 * maximum firing / minimum firing - T 256');

%%

[COEFF1,SCORE1] = princomp(isinormnorm');
[COEFF,SCORE] = princomp(pulsenormnormPSTH');
%[COEFF2,SCORE2] = princomp(tcnormminn');

X = [];
X(:,1) = (sumsumsp(:,2)+sumsumsp(:,3))./(sumsumsp(:,1)+sumsumsp(:,4));
X(:,2) = SCORE1(:,1);
X(:,3) = yyy;

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


 plot(tcnormnorm(:,ismember(on_10,ont1_10)), 'r')
 hold on;
 plot(tcnormnorm(:,~ismember(on_10,ont1_10)), 'b')
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('T1 Comparison', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1comp-tc', gcf)

%%
y = [];
y(1,1) = mean(A32(ismember(on_10,ont1_10)));
y(1,2) = mean(A32(~ismember(on_10,ont1_10)));

y(2,1) = mean(A256(ismember(on_10,ont1_10)));
y(2,2) = mean(A256(~ismember(on_10,ont1_10)));


b = bar(y, 0.5);
set(b(1), 'LineWidth', 2, 'EdgeColor' , [1 0 0], 'FaceColor' , [1 0.7 0.8])
set(b(2), 'LineWidth', 2, 'EdgeColor' , [0 0 1], 'FaceColor' , [.7 0.8 1])


h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('T1 Comparison  - Drifting Grating Spikes', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('Drifting Grating Speed', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('Average total number of spikes', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-10-1/ON/T1/', 'ont1comp-dg', gcf)

%%
%DS
% ON T1
%on t2
%on t3
%off t1
y(5,3) = mean(log(mvnpdf(X(idx==2,:),obj.mu(2,:), obj.Sigma(:,:,2))./(mvnpdf(X(idx==2,:),obj.mu(1,:), obj.Sigma(:,:,1)))));

%%
b = bar(y, 0.5);

set(b(1), 'LineWidth', 2, 'EdgeColor' , [1 0 0], 'FaceColor' , [1 0.7 0.8])
set(b(2), 'LineWidth', 2, 'EdgeColor' , [0 0 1], 'FaceColor' , [.7 0.8 1])
set(b(3), 'LineWidth', 2, 'EdgeColor' , [0 1 0], 'FaceColor' , [.7 1 0.8])
set(gca,'XTickLabel',{'DS', 'ON T1', 'ON T2', 'ON T3', 'OFF T1'})
legend('10-10', '10-15', '10-31')
%%
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Likelihood ratios for all datasets, all types', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);

%%
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/', 'LR', gcf)

%% likelihood ratios
%DS
% ON T1
%on t2
%on t3
%off t1
[7.39750624371868,7.60350298589891,7.48906818841408;26.6735133892520,12.9800814106315,12.0814571640969;7.91109935433406,15.2857514694647,11.3081717855923;5.77024142295366,7.08370473256285,6.56816113091963;20.0159889011665,24.7095978509001,8.42683262291685]

%%
%%
for i = 1:length(ont1_15)
    subplot(6,6,i)
    plot_rf(datarun000_15, ont1_15(i), 'title', false)
end

pause;
close;
figure();
plot_rf_summaries(datarun000_15, ont1_15)

%%

for i = 1:5
    subplot(5,5,i)
    plot_rf(datarun000_31, ont1_31(i), 'title', false)
end
for i = 6:10
    subplot(5,5,i)
    plot_rf(datarun000_31, ont1_31(i), 'title', false)
end
for i = 11:15
    subplot(5,5,i)
    plot_rf(datarun000_31, ont1_31(i), 'title', false)
end
for i = 16:20
    subplot(5,5,i)
    plot_rf(datarun000_31, ont1_31(i), 'title', false)
end
for i = 21:22
    subplot(5,5,i)
    plot_rf(datarun000_31, ont1_31(i), 'title', false)
end
pause;
close;
figure();
plot_rf_summaries(datarun000_31, ont1_31)
%%
plot_rf_summaries(datarun000_10, ont1_10)
pause;
close;
for i = 1:length(ont1_10)
    subplot(5,4,i)
    plot_rf(datarun000_10, ont1_10(i), 'title', false)
end


%%

[1067,1306,1595,1878,2494,2851,5240,5328,6211,6692,6753]

%% on t3
plot_rf_summaries(datarun000_31, [ont3_31, 781,  1201, 1382, 2524, 2747, 2851, 3092, 3182, 3259, 3512, 3858, 4232, 4847, 4877, 6872, 7218 ], 'coordinates', 'monitor');

plot_rf_summaries(datarun000_31, [ont3_31 ], 'coordinates', 'monitor');
%% on t3
[1385,1549,2136,3452,4892,4941,5179,5405,6093,6980,7157,7520]
plot_rf_summaries(datarun000_15, [ont3_15,  800,1368,1996,2237,2941,3017,3185,3410,3755, 6242,6618,4966,3811,3979,4143,4337,7592], 'coordinates', 'monitor');

plot_rf_summaries(datarun000_15, [ont3_15 ], 'coordinates', 'monitor');

%% ont3
[1759,2134,3213,4443,4534,6275,6752,7564]
plot_rf_summaries(datarun000_10, [ont3_10,6425,7384, 468, 842, 1129,1278,1384,1414,1925, 1967,2374,2687,2943,3125,3721,3887,4819,5148,5493,6137,], 'coordinates', 'monitor');
plot_rf_summaries(datarun000_10, [ont3_10], 'coordinates', 'monitor');

%%
color_palette = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 0.5 0.5; 0 0 0; 1 0.5 0; 0.5 1 0; 0 0.5 1; 0.25 0.5 1; 1 0.25 0.5; 0.5 1 0.25; 0.5 0.25 1; 1 0.5 0.25; 0 0.75 0.3; 0.3 0.75 0; 0 0.3 0.75; 0.3 0 0.75; 0.75 0 0.3; 0.75 0.3 0; 0.2 0 0.7;1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 0.5 0.5; 0 0 0; 1 0.5 0; 0.5 1 0; 0 0.5 1; 0.25 0.5 1; 1 0.25 0.5; 0.5 1 0.25; 0.5 0.25 1; 1 0.5 0.25; 0 0.75 0.3; 0.3 0.75 0; 0 0.3 0.75; 0.3 0 0.75; 0.75 0 0.3; 0.75 0.3 0; 0.2 0 0.7;];


%ont1_10 =[77,211,244,1144,1517,2402,2716,3287,3331,4712,5087,5176,5642,5671,5821,6197,6767,6977,7668];


ont1_15 = ont1_15(8:10);

rf1 = zeros(40,80,3);
for i = 1:length(ont1_15)
    rf = get_rf(datarun000_15, ont1_15(i), 'color_transform', color_palette(i,:));
    rf = matrix_scaled_up(rf,1);
    tform = coordinate_transform(datarun000_15,'sta','input','sta scaled','scale',1);
    [xx,yy] = tformfwd(tform,[1 size(rf,2)],[1 size(rf,1)]);
    cell_index  = get_cell_indices(datarun000_15, ont1_15(i));
    if isfield(datarun000_15.stas, 'polarities') && ~isempty(datarun000_15.stas.polarities{cell_index})
        polarity = datarun000_15.stas.polarities{cell_index};
        if polarity == 0
            polarity = 1;
        end
        rf = rf * polarity;
    else
        warning('polarity information unavailable, OFF cells might look like ON cells')
    end
    rf1 = rf1+rf;
        imagesc(xx,yy,norm_image(rf1))%'parent',plot_axes);
        %pause;
end
pause;

hold on;
 plot_rf_summaries(datarun000_15, ont1_15)
 
 %%
 %%
eiplots(datarun000_31, ont1_31, 'color_palette', color_palette, 'cutoff', 0.05, 'pause', true, 'plotei', true)

%%
datarun = datarun000_31;
datarun.cell_types{1,15}.name = 'ont1_31';
datarun.cell_types{1,15}.cell_ids = int32(ont1_31);

%%
cell_spec = {15};

% compute and plot contours for a mosaic of RFs

% filter summary frames (saved to datarun.stas.summaries_filt{})
datarun = get_rfs_filtered(datarun,cell_spec,'verbose',1,'filt_params',struct('filt_type','gauss','radius',.4));

% Do Delaunay Triangulation; for determining RoI for calculating UI
datarun = do_delaunay_tri(datarun, cell_spec, true);
    

% Cut triangles with excessively long edges; indicates where cells are
% likely missing in mosaic.
datarun = cull_delaunay_tri(datarun, cell_spec, 1.9);

% Build Region of Interest from culled Delaunay Triangulation
datarun = build_rf_roi(datarun, cell_spec);

% Get contour threshold that maximizes uniformity index
datarun = maximize_uniformity_index(datarun, cell_spec, 0.3, 0.05, 0.01);

% compute contours (summary frames normalized, then saved to datarun.stas.contours)
datarun = get_rf_contours(datarun, cell_spec, datarun.stas.mosaics{cell_spec{:}}.best_thresh, 'norm_params', struct('method', 'peak'), 'verbose', 1);

% plot contours
plot_rf_summaries(datarun, cell_spec, 'plot_contours', 1, 'foa', 1, 'label', 0);
lock;

% plot simplified contours, with alpha fills
datarun = simplify_rf_contours(datarun, cell_spec);
plot_rf_summaries(datarun, cell_spec, 'plot_contours', 1, 'contours_field', 'rf_contours_simple', 'foa', 1, 'label', 0, 'contour_fill', 'r', 'contour_alpha', 0.5);


% ALTERNATIVELY...
figure
plot_rf_coloring(datarun, cell_spec, 'rfs', 'summaries_filt');

%%

close all;
rstd = [];
meanpix = [];
cellind = get_cell_indices(datarun000_10, ont1_10);
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
[x y] = meshgrid(1:1:size(datarun000_10.stas.rfs{1,1},2), 1:1:size(datarun000_10.stas.rfs{1,1},1));
stddev = 1;
for i = 1:length(cellind)
    rm = [];
    rn = [];
    [rm,rn] = find(stamat{i,1});
    X = [];
    X(:,1) = rn;
    X(:,2) = rm;
    rf = get_rf(datarun000_10, ont1_10(i));
    rf = matrix_scaled_up(rf,1);
    cell_index  = get_cell_indices(datarun000_10, ont1_10(i));
    if isfield(datarun000_10.stas, 'polarities') && ~isempty(datarun000_10.stas.polarities{cell_index})
    polarity = datarun000_10.stas.polarities{cell_index};
        if polarity == 0
            polarity = 1;
        end
        rf = rf * polarity;
    else
        warning('polarity information unavailable, OFF cells might look like ON cells')
    end
    rf = norm_image(rf);
    rf1 = rf(:,:,1);
    mass = rf1(stamat{i,1}~= 0);
    pts = centroid(X, mass);
    S.mu = pts;
    S.Sigma = cov(X);
    obj = gmdistribution.fit(X,1, 'Replicates', 1, 'CovType', 'full', 'Start', S);
    plot(pts(1), pts(2), '+b');
    hold on;
            for a = 1:size(y,2)
                YY(:,a) = mvnpdf([x(:,a) y(:,a)], obj.mu, obj.Sigma); %List of probabilities at each position in the 2-d gaussian
            end
            if (round(obj.mu(1,1)+ stddev*sqrt(obj.Sigma(1,1))) < xx(1,2)) %Check if standard deviation limit is within bounds of the array for plotting
                [rx cx] = find(x==round(obj.mu(1,1)+  stddev*sqrt(obj.Sigma(1,1))));
            else
                [rx cx] = find(x==round(obj.mu(1,1)-  stddev*sqrt(obj.Sigma(1,1))));
            end
            if (round(obj.mu(1,2)+  stddev*sqrt(obj.Sigma(2,2))) < yy(1,2))
                [ry cy] = find(y==round(obj.mu(1,2)+  stddev*sqrt(obj.Sigma(2,2))));
            else
                [ry cy] = find(y==round(obj.mu(1,2)-  stddev*sqrt(obj.Sigma(2,2))));
            end
            v = [YY(ry(1,1),cx(1,1)), YY(ry(1,1),cx(1,1))];
            [C h] = contour(x,y,YY,  v, 'b', 'LineWidth', 2);   
end
                set(gca,'Ydir','reverse');%'Ydir','reverse')



%         ei(abs(ei)<params.cutoff) = 0; %Set all electrodes below cutoff to zero
%         ei_frame = get_ei_max_frame(ei, 1);
%         elecofei = find(ei_frame);
%         ei(abs(ei)<1) = 0;
%         [r2 c2] = find(ei);
%         elec_cent = r2;
%      
%         allx = datarun.ei.position(elecofei,1); %get all x and y positions of electrodes above cutoff
%         ally = datarun.ei.position(elecofei,2);
%         X = [];
%         X(:,1) = allx;
%         X(:,2) = ally;
        %mu = [datarun.ei.position(elec_cent,1), datarun.ei.position(elec_cent,2)]; %Source electrode is the mean
%         if(strcmp(params.fit_type, 'gmdist'))
%             obj = gmdistribution(mu,cov(X));
%         elseif(strcmp(params.fit_type, 'fitgm'))
%              obj = gmdistribution.fit(X,1, 'Replicates',1,'CovType', 'full'); % If you want to fit the EI - no initial conditions
%         elseif(strcmp(params.fit_type, 'fitgmic')) %Fit EI with initial conditions
            %S.mu = mu;
            %S.mu = mean(X);
            
            
%         end  
%         plot(datarun.ei.position(elec_cent,1), datarun.ei.position(elec_cent,2), '+r');
%         hold on;
%         %plot(pts(1), pts(2), '+b');
%         if (strcmp(params.contour_type, 'ezcontour'))
%             ezcontour(@(x,y)pdf(obj,[x y]),[datarun.ei.array_bounds_x(1,1) datarun.ei.array_bounds_x(1,2)],[datarun.ei.array_bounds_y(1,1) datarun.ei.array_bounds_y(1,2)], params.contourgrid);
%             hold on;
%         elseif (strcmp(params.contour_type, 'simplecontour'))
%             for a = 1:size(y,2)
%             YY(:,a) = mvnpdf([x(:,a) y(:,a)], obj.mu, obj.Sigma); %List of probabilities at each position in the 2-d gaussian
%             end
%             contour(x,y,YY, params.numcontours);
%             hold on;
%         elseif (strcmp(params.contour_type, 'sdcontour'))
%             for a = 1:size(y,2)
%                 YY(:,a) = mvnpdf([x(:,a) y(:,a)], obj.mu, obj.Sigma); %List of probabilities at each position in the 2-d gaussian
%             end
%             if (round(obj.mu(1,1)+ params.stddev*std(allx)) < datarun.ei.array_bounds_x(1,2)) %Check if standard deviation limit is within bounds of the array for plotting
%                 [rx cx] = find(x==round(obj.mu(1,1)+ params.stddev*std(allx)));
%             else
%                 [rx cx] = find(x==round(obj.mu(1,1)- params.stddev*std(allx)));
%             end
%             if (round(obj.mu(1,2)+ params.stddev*std(ally)) < datarun.ei.array_bounds_y(1,2))
%                 [ry cy] = find(y==round(obj.mu(1,2)+ params.stddev*std(ally)));
%             else
%                 [ry cy] = find(y==round(obj.mu(1,2)- params.stddev*std(ally)));
%             end
%             v = [YY(ry(1,1),cx(1,1)), YY(ry(1,1),cx(1,1))];
%             [C h] = contour(x,y,YY,  v, 'b', 'LineWidth', 2);
%         end
    end


 %%
 
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
                hold on;
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

%%

 [x y] = meshgrid(datarun.ei.array_bounds_x(1,1):1:datarun.ei.array_bounds_x(1,2), datarun.ei.array_bounds_y(1,1):1:datarun.ei.array_bounds_y(1,2));
        ei = temp_ei/max(abs(temp_ei(:)));
        ei(abs(ei)<params.cutoff) = 0; %Set all electrodes below cutoff to zero
        ei_frame = get_ei_max_frame(ei, 1);
        elecofei = find(ei_frame);
        ei(abs(ei)<1) = 0;
        [r2 c2] = find(ei);
        elec_cent = r2;
        allx = datarun.ei.position(elecofei,1); %get all x and y positions of electrodes above cutoff
        ally = datarun.ei.position(elecofei,2);
        X = [];
        X(:,1) = allx;
        X(:,2) = ally;
        %mu = [datarun.ei.position(elec_cent,1), datarun.ei.position(elec_cent,2)]; %Source electrode is the mean
%         if(strcmp(params.fit_type, 'gmdist'))
%             obj = gmdistribution(mu,cov(X));
%         elseif(strcmp(params.fit_type, 'fitgm'))
%              obj = gmdistribution.fit(X,1, 'Replicates',1,'CovType', 'full'); % If you want to fit the EI - no initial conditions
%         elseif(strcmp(params.fit_type, 'fitgmic')) %Fit EI with initial conditions
            %S.mu = mu;
            %S.mu = mean(X);
            mass = abs(ei_frame(ei_frame~= 0));
            pts = centroid(X, mass);
            S.mu = pts;
            S.Sigma = cov(X);
            obj = gmdistribution.fit(X,1, 'Replicates', 1, 'CovType', 'full', 'Start', S);
        end  
        plot(datarun.ei.position(elec_cent,1), datarun.ei.position(elec_cent,2), '+r');
        hold on;
        %plot(pts(1), pts(2), '+b');
        if (strcmp(params.contour_type, 'ezcontour'))
            ezcontour(@(x,y)pdf(obj,[x y]),[datarun.ei.array_bounds_x(1,1) datarun.ei.array_bounds_x(1,2)],[datarun.ei.array_bounds_y(1,1) datarun.ei.array_bounds_y(1,2)], params.contourgrid);
            hold on;
        elseif (strcmp(params.contour_type, 'simplecontour'))
            for a = 1:size(y,2)
            YY(:,a) = mvnpdf([x(:,a) y(:,a)], obj.mu, obj.Sigma); %List of probabilities at each position in the 2-d gaussian
            end
            contour(x,y,YY, params.numcontours);
            hold on;
        elseif (strcmp(params.contour_type, 'sdcontour'))
            for a = 1:size(y,2)
                YY(:,a) = mvnpdf([x(:,a) y(:,a)], obj.mu, obj.Sigma); %List of probabilities at each position in the 2-d gaussian
            end
            if (round(obj.mu(1,1)+ params.stddev*std(allx)) < datarun.ei.array_bounds_x(1,2)) %Check if standard deviation limit is within bounds of the array for plotting
                [rx cx] = find(x==round(obj.mu(1,1)+ params.stddev*std(allx)));
            else
                [rx cx] = find(x==round(obj.mu(1,1)- params.stddev*std(allx)));
            end
            if (round(obj.mu(1,2)+ params.stddev*std(ally)) < datarun.ei.array_bounds_y(1,2))
                [ry cy] = find(y==round(obj.mu(1,2)+ params.stddev*std(ally)));
            else
                [ry cy] = find(y==round(obj.mu(1,2)- params.stddev*std(ally)));
            end
            v = [YY(ry(1,1),cx(1,1)), YY(ry(1,1),cx(1,1))];
            [C h] = contour(x,y,YY,  v, 'b', 'LineWidth', 2);
        end
        end

    %% 2012-10-31 extra cells
%%on t3
plot_time_courses(datarun000_31, [[ ont3_31, 1201, 1382, 2524,3182, 3259, 3512, 3858,  4847, 4877, 6872, 7218]],'all', true, 'bw' ,true);
plot_rf_summaries(datarun000_31, [[ ont3_31, 1201, 1382, 2524,3182, 3259, 3512, 3858,  4847, 4877, 6872, 7218]],'coordinates', 'monitor');
%1067 and 1595 in on t3 have different EIs but exact same RF
%3092 4232 removed because tcs don't fit


plot_rf_summaries(datarun000_31, [1201, 1382, 2524,3182, 3259, 3512, 3858,  4847, 4877, 6872, 7218], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3extrarf', gcf)

plot_time_courses(datarun000_31,[1201, 1382, 2524,3182, 3259, 3512, 3858,  4847, 4877, 6872, 7218], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T3/', 'ont3extratc', gcf)

%% on t2 - 10-31
close all;
plot_rf_summaries(datarun000_31, [[ ont2_31, 541, 1561, 2015, 3362, 3903, 4205, 4516, 4531, 4576,5150, 5162,5611, 6439, 7381 ]],'coordinates', 'monitor');
%5716 and 6002 in on t2 are the same cell - double overlapping part of rf
%6588 looks a little different - tc
figure();
plot_time_courses(datarun000_31,[ ont2_31, 541, 1561, 2015, 3362, 3903, 4205, 4516, 4531, 4576,5150, 5162,5611, 6439,  7381 ],'all', true, 'bw' ,true);
ismember([ 541, 1561, 2015, 3362, 3903, 4205, 4516, 4531, 4576, 5162,5611, 6439, 7381, 5150], cellids_31)

%%
plot_rf_summaries(datarun000_31, [541, 1561, 2015, 3362, 3903, 4205, 4516, 4531, 4576,5150, 5162,5611, 6439, 7381 ], 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2extrarf', gcf)

plot_time_courses(datarun000_31,[541, 1561, 2015, 3362, 3903, 4205, 4516, 4531, 4576,5150, 5162,5611, 6439, 7381 ], 'all', true, 'bw', true);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sneharavi/Documents/from sravi/My Documents/My Publications/Rat Classification 2014/2012-10-31-0/ON/T2/', 'ont2extratc', gcf)






%% ont1 - 10-31
close all;
plot_rf_summaries(datarun000_31, [[ ont1_31,3106, 3979,2552,3391 ]],'coordinates', 'monitor');
figure();
plot_time_courses(datarun000_31,[ ont1_31, 3106, 3979,2552,3391],'all', true, 'bw' ,true);
plot_rf_summaries(datarun000_31, [[ ont1_31]],'coordinates', 'monitor');

%%

close all;
a = [ont1_10(8:10)];
rstd = [];
meanpix = [];
cellind = get_cell_indices(datarun000_10, [a]);
stamat = cell(length(cellind),1);

color_palete = [1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;];
 
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
[x y] = meshgrid(1:1:size(datarun000_10.stas.rfs{1,1},2), 1:1:size(datarun000_10.stas.rfs{1,1},1));
rf2 = zeros(40,80,3);
for i = 1:length(cellind)
    rm = [];
    rn = [];
    [rm,rn] = find(stamat{i,1});
    X = [];
    X(:,1) = rn;
    X(:,2) = rm;
    rf = get_rf(datarun000_10, a(i));
        rf = get_rf(datarun000_10, a(i)); %'color_transform', color_palete(i,:));
        tform = coordinate_transform(datarun000_10,'sta','input','sta scaled','scale',1);
    [xx,yy] = tformfwd(tform,[1 size(rf,2)],[1 size(rf,1)]);
    rf = matrix_scaled_up(rf,1);
    cell_index  = cellind(i);
    if isfield(datarun000_10.stas, 'polarities') && ~isempty(datarun000_10.stas.polarities{cell_index})
    polarity = datarun000_10.stas.polarities{cell_index};
        if polarity == 0
            polarity = 1;
        end
        rf = rf * polarity;
    else
        warning('polarity information unavailable, OFF cells might look like ON cells')
    end
        rf1 = norm_image(rf);
        rf1 = rf1(:,:,1);
    rf1(find(~stamat{i,1})) = 0.5;
    rf1 = rf1*2;
    rf1 = rf1-1;
    
    [num_rows, num_cols, num_pages] = size(rf1);
    reshaped_rf = reshape(rf1, [], num_pages);
    transformed_rf = reshaped_rf *color_palete(i,:);
    rf1 = reshape(transformed_rf, num_rows, num_cols, []);
  rf2 = rf2 + rf1;
        pause;


end

       imagesc(xx,yy,rf2);

       
       %%

close all;
a = [ont1_31 ];
rstd = [];
meanpix = [];
cellind = get_cell_indices(datarun000_31, [ont1_31 3106, 3979,2552,3391]);
stamat = cell(length(cellind),1);
 
for i = 1:length(cellind)
    B = [];
    B = datarun000_31.stas.rfs{cellind(i), 1}(:)';
    meanpix(i) = mean(B);
    rstd(i) = robust_std(B, [1]);
    stamat{i,1} = zeros(size(datarun000_31.stas.rfs{cellind(i),1},1),size(datarun000_31.stas.rfs{cellind(i),1},2));
    for j = 1:size(datarun000_31.stas.rfs{cellind(i),1},1)
            for k = 1:size(datarun000_31.stas.rfs{cellind(i),1},2)
                if (datarun000_31.stas.rfs{cellind(i),1}(j,k) >= meanpix(i) + 5*rstd(i))
                    stamat{i,1}(j,k) = 1;
                end
            end
    end
 
    
end
    
 

close all;
[x y] = meshgrid(1:1:size(datarun000_31.stas.rfs{1,1},2), 1:1:size(datarun000_31.stas.rfs{1,1},1));
stddev = 1;
for i = 1:length(cellind)
    rm = [];
    rn = [];
    [rm,rn] = find(stamat{i,1});
    X = [];
    X(:,1) = rn;
    X(:,2) = rm;
    rf = get_rf(datarun000_31, a(i));
    rf = matrix_scaled_up(rf,1);
    cell_index  = cellind(i);
    if isfield(datarun000_31.stas, 'polarities') && ~isempty(datarun000_31.stas.polarities{cell_index})
    polarity = datarun000_31.stas.polarities{cell_index};
        if polarity == 0
            polarity = 1;
        end
        rf = rf * polarity;
    else
        warning('polarity information unavailable, OFF cells might look like ON cells')
    end
    rf = norm_image(rf);
    rf1 = rf(:,:,1);
    mass = rf1(stamat{i,1}~= 0);
    pts = centroid(X, mass);
    S.mu = pts;
    S.Sigma = cov(X);
    obj = gmdistribution.fit(X,1, 'Replicates', 1, 'CovType', 'full', 'Start', S);
    plot(pts(1), pts(2), 'r+');
    hold on;
    plot(pts(1), pts(2), 'Marker', 'o', 'MarkerSize', 20);
    hold on;
%             for a = 1:size(y,2)
%                 YY(:,a) = mvnpdf([x(:,a) y(:,a)], obj.mu, obj.Sigma); %List of probabilities at each position in the 2-d gaussian
%             end
%             if (round(obj.mu(1,1)+ stddev*sqrt(obj.Sigma(1,1))) < xx(1,2)) %Check if standard deviation limit is within bounds of the array for plotting
%                 [rx cx] = find(x==round(obj.mu(1,1)+  stddev*sqrt(obj.Sigma(1,1))));
%             else
%                 [rx cx] = find(x==round(obj.mu(1,1)-  stddev*sqrt(obj.Sigma(1,1))));
%             end
%             if (round(obj.mu(1,2)+  stddev*sqrt(obj.Sigma(2,2))) < yy(1,2))
%                 [ry cy] = find(y==round(obj.mu(1,2)+  stddev*sqrt(obj.Sigma(2,2))));
%             else
%                 [ry cy] = find(y==round(obj.mu(1,2)-  stddev*sqrt(obj.Sigma(2,2))));
%             end
%             v = [YY(ry(1,1),cx(1,1)), YY(ry(1,1),cx(1,1))];
%             [C h] = contour(x,y,YY,  v, 'b', 'LineWidth', 2);   
end
                set(gca,'Ydir','reverse');%'Ydir','reverse')
                
                %%
                
                close all;
a = [ont1_10(8:10)];
rstd = [];
meanpix = [];
cellind = get_cell_indices(datarun000_10, [a]);
stamat = cell(length(cellind),1);

color_palete = [1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1;];
 
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
[x y] = meshgrid(1:1:size(datarun000_10.stas.rfs{1,1},2), 1:1:size(datarun000_10.stas.rfs{1,1},1));
xq = reshape(x, 3200,1);
yq = reshape(y, 3200,1);
rf2 = zeros(40,80,3);
for i = 1:length(cellind)
    rm = [];
    rn = [];
    rmm = [];
    rnn = [];
    DT = [];
    kr = [];
    points = [];
    rf1 = [];
    [rm,rn] = find(stamat{i,1}); %finds all the significant pixels
        [rmm,rnn] = pix_border(rm,rn); %adding everything by .5 and subtracting by 0.5 - i think it is calculating the 4 coordinate values for each pixel - edge correction
        DT = delaunayTriangulation(rmm,rnn); %prediction is still correct with inputs rm and rn, but edges need correcting
        if ~isempty(DT.ConnectivityList)
            [kr v(1,i)] = convexHull(DT);
        end
        [in on] = inpolygon(xq, yq, DT.Points(kr,2), DT.Points(kr,1));
        X = [];
        X(:,1) = xq(in);
        X(:,2) = yq(in);
        %figure();
        %plot(DT.Points(kr,2), DT.Points(kr,1), 'Color',[rand(1) rand(1) rand(1)]);
        %hold on
        %plot(xq(in),yq(in),'r+') % points inside
        %plot(xq(~in),yq(~in),'wo') % points outside
        %hold off
        %pause;
        rf = get_rf(datarun000_10, a(i)); %'color_transform', color_palete(i,:));
        tform = coordinate_transform(datarun000_10,'sta','input','sta scaled','scale',1);
        [xx,yy] = tformfwd(tform,[1 size(rf,2)],[1 size(rf,1)]);
        rf = matrix_scaled_up(rf,1);
        cell_index  = cellind(i);
        if isfield(datarun000_10.stas, 'polarities') && ~isempty(datarun000_10.stas.polarities{cell_index})
            polarity = datarun000_10.stas.polarities{cell_index};
            if polarity == 0
                polarity = 1;
            end
            rf = rf * polarity;
        else
            warning('polarity information unavailable, OFF cells might look like ON cells')
        end
        rf1 = norm_image(rf);
        rf1 = rf1(:,:,1);
    rf1(find(~stamat{i,1})) = 0;
        %rf1(yq(~in),xq(~in)) = 0.5;
        %rf1 = rf1*2;
        %rf1 = rf1-1;
        [num_rows, num_cols, num_pages] = size(rf1);
        reshaped_rf = reshape(rf1, [], num_pages);
%         reshaped_rf(~in) = 0.5;
%         reshaped_rf = reshaped_rf*2;
%         reshaped_rf  = reshaped_rf -1;
        %reshaped_rf(~in) = 0;
        transformed_rf = reshaped_rf *color_palete(i,:);
        rf1 = reshape(transformed_rf, num_rows, num_cols, []);
        rf2 = rf2 + rf1;
        imagesc(xx,yy,rf2);
        pause;
end



%%

%%
for i = 1:length(ont1_15)
[T, psth01, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 12, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 12, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 12, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 12, 'bin_size', 0.5);
subplot(2,2,1);plot(0:0.1:12, psth01);title('Bin Size 0.1');
subplot(2,2,2);plot(0:0.2:12, psth02);title('Bin Size 0.2');
subplot(2,2,3);plot(0:0.25:12, psth025);title('Bin Size 0.25');
subplot(2,2,4);plot(0:0.5:12, psth05);title('Bin Size 0.5');
pause;
end
close all;
%%
for i = 1:length(ont1_31)
[T, psth01, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont1_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont1_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont1_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont1_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.5);
subplot(2,2,1);plot(0:0.1:10, psth01);title('Bin Size 0.1');
subplot(2,2,2);plot(0:0.2:10, psth02);title('Bin Size 0.2');
subplot(2,2,3);plot(0:0.25:10, psth025);title('Bin Size 0.25');
subplot(2,2,4);plot(0:0.5:10, psth05);title('Bin Size 0.5');
pause;
end
close all;
%%
for i = 1:length(ont1_10)
[T, psth01, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont1_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont1_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont1_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont1_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.5);
subplot(2,2,1);plot(0:0.1:10, psth01);title('Bin Size 0.1');
subplot(2,2,2);plot(0:0.2:10, psth02);title('Bin Size 0.2');
subplot(2,2,3);plot(0:0.25:10, psth025);title('Bin Size 0.25');
subplot(2,2,4);plot(0:0.5:10, psth05);title('Bin Size 0.5');
pause;
end
close all;
%%
for i = 1:length(ont2_15)
[T, psth01, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont2_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont2_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont2_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont2_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.5);
subplot(2,2,1);plot(0:0.1:8, psth01);title('Bin Size 0.1');
subplot(2,2,2);plot(0:0.2:8, psth02);title('Bin Size 0.2');
subplot(2,2,3);plot(0:0.25:8, psth025);title('Bin Size 0.25');
subplot(2,2,4);plot(0:0.5:8, psth05);title('Bin Size 0.5');
pause;
end
close all;
%%
for i = 1:length(ont2_31)
[T, psth01, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont2_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont2_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont2_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont2_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.5);
subplot(2,2,1);plot(1:0.1:7, psth01(11:71));title('Bin Size 0.1');
subplot(2,2,2);plot(1:0.2:7, psth02(6:36));title('Bin Size 0.2');
subplot(2,2,3);plot(1:0.25:7, psth025(5:29));title('Bin Size 0.25');
subplot(2,2,4);plot(1:0.5:7, psth05(3:15));title('Bin Size 0.5');
pause;
end
close all;
%%
for i = 1:length(ont2_10)
[T, psth01, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont2_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont2_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont2_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont2_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.5);
subplot(2,2,1);plot(0:0.1:10, psth01);title('Bin Size 0.1');
subplot(2,2,2);plot(0:0.2:10, psth02);title('Bin Size 0.2');
subplot(2,2,3);plot(0:0.25:10, psth025);title('Bin Size 0.25');
subplot(2,2,4);plot(0:0.5:10, psth05);title('Bin Size 0.5');
pause;
end
close all;
%%
for i = 1:length(ont3_15)
[T, psth01, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont3_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont3_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont3_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont3_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.5);
subplot(2,2,1);plot(0:0.1:8, psth01);title('Bin Size 0.1');
subplot(2,2,2);plot(0:0.2:8, psth02);title('Bin Size 0.2');
subplot(2,2,3);plot(0:0.25:8, psth025);title('Bin Size 0.25');
subplot(2,2,4);plot(0:0.5:8, psth05);title('Bin Size 0.5');
pause;
end
close all;
%%
for i = 1:length(ont3_31)
[T, psth01, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont3_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont3_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont3_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont3_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,5)), 'stop', 10, 'bin_size', 0.5);
subplot(2,2,1);plot(0:0.1:10, psth01);title('Bin Size 0.1');
subplot(2,2,2);plot(0:0.2:10, psth02);title('Bin Size 0.2');
subplot(2,2,3);plot(0:0.25:10, psth025);title('Bin Size 0.25');
subplot(2,2,4);plot(0:0.5:10, psth05);title('Bin Size 0.5');
pause;
end
close all;
%%
for i = 1:length(ont3_10)
[T, psth01, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont3_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.1);
[T, psth02, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont3_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.2);
[T, psth025, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont3_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.25);
[T, psth05, bins] = get_psth_sr(datarun002_10.spikes{get_cell_indices(datarun002_10, ont3_10(1,i)),1},datarun002_10.stimulus.triggers(ismember(datarun002_10.stimulus.trial_list,2)), 'stop', 10, 'bin_size', 0.5);
subplot(2,2,1);plot(0:0.1:10, psth01);title('Bin Size 0.1');
subplot(2,2,2);plot(0:0.2:10, psth02);title('Bin Size 0.2');
subplot(2,2,3);plot(0:0.25:10, psth025);title('Bin Size 0.25');
subplot(2,2,4);plot(0:0.5:10, psth05);title('Bin Size 0.5');
pause;
end
close all;


       %%
       for i = 1:length(ont1_15)
[T, psth01, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.1);
[T, psth0123, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,7)), 'stop', 8, 'bin_size', 0.1);

% [T, psth02, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 12, 'bin_size', 0.2);
% [T, psth025, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 12, 'bin_size', 0.25);
% [T, psth05, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 12, 'bin_size', 0.5);
% subplot(2,2,1);plot(2:0.1:6, psth01(21:61));title('Bin Size 0.1');
% hold on;
% plot(2:0.1:6, psth0123(21:61), 'r');%title('Bin Size 0.1');
% hold off;
y = psth01(21:61);
[ymx,imx]=max(y);
yy(1:length(imx:length(y))) = y(imx:end);
yy(length(imx:length(y))+1: length(y)) = y(1:imx-1);
plot(2:0.1:6,yy);

% subplot(2,2,2);plot(0:0.2:12, psth02);title('Bin Size 0.2');
% subplot(2,2,3);plot(0:0.25:12, psth025);title('Bin Size 0.25');
% subplot(2,2,4);plot(0:0.5:12, psth05);title('Bin Size 0.5');
pause;
end
close all;


%%
for i = 1:length(ont1_15)
[T, psth01, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.1);
[T, psth0123, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,7)), 'stop', 8, 'bin_size', 0.1);

[T, psth0125, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,19)), 'stop', 8, 'bin_size', 0.25);
[T, psth012325, bins] = get_psth_sr(datarun002_15.spikes{get_cell_indices(datarun002_15, ont1_15(1,i)),1},datarun002_15.stimulus.triggers(ismember(datarun002_15.stimulus.trial_list,7)), 'stop', 8, 'bin_size', 0.25);

y = [];ymx=[];imx=[];yy = [];
y = psth01(21:61);
[ymx,imx]=max(y);
yy = circshift(y,[0,length(y)-imx+1]);

x = [];xmx=[];imx = []; xx = [];
x = psth0123(21:61);
[xmx,imx]=max(x);
xx = circshift(x,[0,length(x)-imx+1]);

subplot(2,3,1);plot(2:0.1:6,yy, 'b');
hold on;
plot(2:0.1:6,xx, 'r');
hold off; legend('direction 1', 'direction 2');
title('Max-BinSize 0.1');

subplot(2,3,3);plot(2:0.1:6,((xx+yy)/2), 'b');hold on;

y = [];ymx=[];imx=[];yy = [];
y = psth01(21:61);
[ymx,imx]=min(y);
yy = circshift(y,[0,length(y)-imx+1]);

x = [];xmx=[];imx = []; xx = [];
x = psth0123(21:61);
[xmx,imx]=min(x);
xx = circshift(x,[0,length(x)-imx+1]);

subplot(2,3,2);plot(2:0.1:6,yy, 'b');
hold on;
plot(2:0.1:6,xx, 'r');legend('direction 1', 'direction 2');
hold off;
title('Min-BinSize 0.1');

subplot(2,3,3);plot(2:0.1:6,((xx+yy)/2), 'r'); legend('max-ave', 'min-ave');hold off;


y = [];ymx=[];imx=[];yy = [];
y = psth0125(9:25);
[ymx,imx]=max(y);
yy = circshift(y,[0,length(y)-imx+1]);

x = [];xmx=[];imx = []; xx = [];
x = psth012325(9:25);
[xmx,imx]=max(x);
xx = circshift(x,[0,length(x)-imx+1]);

subplot(2,3,4);plot(2:0.25:6,yy, 'b');
hold on;
plot(2:0.25:6,xx, 'r');legend('direction 1', 'direction 2');
hold off;
title('Max-BinSize 0.25');

subplot(2,3,6);plot(2:0.25:6,((xx+yy)/2), 'b');hold on;


y = [];ymx=[];imx=[];yy = [];
y = psth0125(9:25);
[ymx,imx]=min(y);
yy = circshift(y,[0,length(y)-imx+1]);

x = [];xmx=[];imx = []; xx = [];
x = psth012325(9:25);
[xmx,imx]=min(x);
xx = circshift(x,[0,length(x)-imx+1]);

subplot(2,3,5);plot(2:0.25:6,yy, 'b');
hold on;
plot(2:0.25:6,xx, 'r');legend('direction 1', 'direction 2');
hold off;
title('Min-BinSize 0.25')

subplot(2,3,6);plot(2:0.25:6,((xx+yy)/2), 'r'); legend('max-ave', 'min-ave');hold off;

pause;

end
close all;

%% Drifitng grating PSTH code
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002_31, ont3_31, 0);

ind = find(StimComb(:,2)==256);
zrind = find(StimComb(:,3)==0);
zeroind = intersect(ind, zrind);
psthall = [];
for i = 1:length(ont3_31)
    psthzero = [];
    psth = [];
    [T, psthzero, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont3_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,zeroind)), 'stop', 8, 'bin_size', 0.1);
    psthzerocut = psthzero(21:61);
    for k = 1:length(ind)
        psth01 = [];
        [T, psth01, bins] = get_psth_sr(datarun002_31.spikes{get_cell_indices(datarun002_31, ont3_31(1,i)),1},datarun002_31.stimulus.triggers(ismember(datarun002_31.stimulus.trial_list,ind(k))), 'stop', 8, 'bin_size', 0.1);
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
dgresp_31(1,:) = mean(psthall);
dgresp_31(2,:) = mean(psthall);
dgresp_31(3,:) = mean(psthall);

plot(dgresp_31','DisplayName','dgresp')
legend('t1', 't2', 't3');
title('2012-10-31');

%% Calculate Response Dominance Index (RDI) for cells
% RDI = Ron - Roff / Ron + Roff 
%Check if cells are on or off or on-off

close all;

wh = datarun001_31.triggers(1:4:length(datarun001_31.triggers), 1); 
gr = datarun001_31.triggers(2:4:length(datarun001_31.triggers),1); 
[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun001_31, get_cell_indices(datarun001_31,ds_31), 0, '/0', wh, gr, 10, false,0.1);
binSize = 0.1:0.1:10; %Bin Size 100 ms
psthnorm = [];
psthind = [];
psthindnorm = [];
b = 0.5; 
g = b+0.5;
w = g+0.5;
rdi = [];
ron = [];
roff = [];

for i = 1:length(ds_31)
     psthind = sum(nhist{i,1})/length(wh);
    psthindnorm = psthind./norm(psthind); %Normalized PSTH
    ron = max(psthindnorm(1:10)) + max(psthindnorm(51:60)); %Take maximum spike rate during 1st 1 second of light pulse
    roff = max(psthindnorm(31:40))+max(psthindnorm(81:90));
    rdi(i) = (ron-roff)/(ron+roff);
    psthnorm(i,:) = psthindnorm;
    %plot(binSize,psthindnorm);
    %pause;
end
histogram(rdi, 20)
title('RDI distribution - 2012-10-31')
length(find(rdi > 0.6)) + length(find(rdi < -0.6))
    

