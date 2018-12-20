[datarun002_11] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-11-13-0/data001-3600-7200s/', 'data002-map/data002-map', 1, '/stimuli/s02', 0);
[datarun000_11] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-11-13-0/data001-3600-7200s/', 'data001-map/data001-map', 0, 0, 1);

[datarun000_03] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-03-01-0/data000-3600-7200s/', 'data000-map/data000-map', 0, 0, 1);
[datarun002_03] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-03-01-0/data000-3600-7200s/', 'data003-map/data003-map', 1, '/stimuli/s03', 0);

[datarun000_10] = load_dsdata('/Analysis/sravi/2012-10-15-0/data000-3600-7200s/datamaps-sr-model/', 'data000-map/data000-map', 0, 0, 1);

datarun000_12 = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2011-12-23-0/data001-da/', 'data001-map/data001-map', 0, 0, 1);

datarun000_03 = get_autocorrelations(datarun000_03, 'all');
datarun000_11 = get_autocorrelations(datarun000_11, 'all');
datarun000_10 = get_autocorrelations(datarun000_10, 'all');
datarun000_12 = get_autocorrelations(datarun000_12, 'all');


ONT1_03 = [127 140 349 2990 3318 4428 4684 4773 5207 5584 7640];
ONT1_11 = [200 456 798 2508 2569 2732 2973 3065 3976 4248 4760 5285 5912 5944 6362 6723 6767 7144 7385];
ONT1_10 = [4 154 692 860 1111 1339 1786 1877 1892 2011 2161 2208 2461 2851 3002 3319 3422 3691 3994 4774 5073 5088 5359 5567 5853 6542 6826 7261 7442 7487 7532];
ONT1_12 = [49 124 197 407 484 558 602 858 993 1370 1623 1847 1954 2011 2419 2927 3196 3346 3541 3827 4412 4682 4816 5236 5313 5852 6226 6363 6468 6693 6886 6963 7144 7201 7276 7474 7518];



ONT2_11 = [49 678 721 1189 1336 1351 1758 1771 2449 2869 2941 3303 3574 3664 4187 4219 4729 4741 4954 5345 5446 5659 5807 5853 6363 6527 6679 6946 7037 7236 7277 7622];
ONT2_03 = [259 331 616 1022 1126 1381 1531 2716 3091 3286 4306 4727 4997 6499 6917 6992 7066 7276];
ONT2_10 = [272 454 783 888 979 1067 1398 2026 2417 2929 3198 3244 3303 3559 3634 3917 3946 3991 4246 4562 4998 5657 5898 5927 6106 6170 6455 6512 6722 6903 7067 7203 7354 7562];
ONT2_12 = [16  302  647  856  887 1129 1143 1231 1456 1667 1712 2283 2476 2521 2568 2838 2987 3181 3527 3677 3872 4128 4220 4232 4621 5312 5447 5538 6048 6138 6393 6664 6991 7381];

ONT3_11 = [333 484 586 918 2165 2179 2626 3273 3556 3635 3766 4876 4982 5134 5375 5447 5461 5749 6034 7397];
ONT3_03 = [286 320 603 902 1788 2567 2866 3797 5329 5551 5613  6289 6647 6814];
ONT3_10 = [1385 1549 2136 3452 4892 4941 5179 5405 6093 6980 7157 7520];

ONT4_11 = [512 875 935 1098 1414 1876 3004 4007 4052 4054 4068 4070 4247 4516 4592 4622 5671 6334 6786 6918 6979 7114];
ONT4_03 = [126 137 243 364 1655 1818 2089 2117 2146 2180 2869 3169 3711 4084 4473 5462 5763 6124 7040];
ONT4_10 = [182 528 857 1532 1591 3530 3586 4231 6152 6155 6304 6439 6737];


ONT5_03 = [1966 2042 2883 3140 3768];
ONT5_11 = [2074 2763 3244 6002];


OFFT1_03 = [469 543 904 1426 1713 2013 2463 2974 3541 3721 4128 4594 5266 5597 5686 6211 6378 6661 7531 7562];
OFFT1_11 = [648 662 1488 1699 1894 2686 2779 4491 4906 7322];

OFFT1_10 = [46 751 782 1051 1381 1637 2401 3933 4006 4324 4503 5223 5446 5641 5658 5836 6451 6811 7306 768 1081 1652 2177 2686 3287 3589 4591 6391 6976 7021 7471];

OFFT1_12 = [3 214 363 497 694 782 1157 1518 1981 2432 2747 2883 2942 3211 3408 3603 3708 4052 4067 4398 4413 4546 4639 4848 5057 5418 5628 6301 6406 6575 6632 6663 7022 7203 7353 7413 7591];



OFFT2_03 = [275 708 2147 2601 2674 2854 3185 3560 3918 4566 4759 5194 5596 5886 6153 6425 6786 6874 7041 7068 7174 7204 7383];
OFFT2_11 = [109 241 528 887 968 997 1163 1417 1549 1893 2451 3034 3347 3544 3723 3842 4021 4352 4505 4744 5104 5194 5450 5599 6106 6256 6466 6586 6845 6874 7250];

OFFT2_10 = [94 347 514 586 872 1098 1310 1384 1581 1893 1908 1999 2192 2357 2462 2522 2746 2809 2868 2881 3061 3139 3226 3258 3601 3692 3813 4098 4442 4459 4486 4864 4999 5071 5433 5464 5569 5851 6034 6140 6196 6422 6589 6721 6828 7234 7503 7667];

 OFFT2_12 = [46 64 80 603 679 1008 1430 1516 2313 2689 3091 3362 3514 3858 4023 4249 4292 4652 4982 5027 5221 5253 5523 5778 6077 6122 6288 6662 6738 7068 7368];


OFFT3_03 = [138 661 1652 1757 3046 3047 5432 5478 5868 7006 7325 7367 7368 7503 7623];
OFFT3_11 = [467 1051 1263 2357 2911 2926 4488 5281 5701 6722 7666];

OFFT4_03 = [365 1156 1400 1473 3769 3961 6214 6421 6571 6753];
OFFT4_11 = [93 485 661 1910 2297 2418 2419 2642 3151 4295 4879 5147 5297 6077 6079 6333 7339];

OFFT4_10 = [991 1156 4234 4278 4487 5733 6286 6931];
OFFT4_12 = [2238 2446 2869 3647 4352 4475 5012 5314 7174 7354];

OFFT5_03 = [122 529 1906 3362 3497 4190 4351 5116 5296 5599 6258 6694 7232 7578];
OFFT5_11 = [962 1831 2476 2851 2855 3316 4457 6782 6933];

timcs_03 = [];
timcs_10 = [];
timcs_11 = [];
timcs_12 = [];
  timcs_03 = get_time_courses_matrix(datarun000_03, OFFT2_03);
 timcs_10 = get_time_courses_matrix(datarun000_10, OFFT4_10);
 timcs_11 = get_time_courses_matrix(datarun000_11, OFFT2_11);
timcs_12 = get_time_courses_matrix(datarun000_12, OFFT4_12);
 

 
  for b = 1:length(OFFT2_03)
   timcs_03(:,b) = timcs_03(:,b)./norm(timcs_03(:,b));    %normalize time course
  end 
   for b = 1:length(OFFT4_10)
   timcs_10(:,b) = timcs_10(:,b)./norm(timcs_10(:,b));    %normalize time course
   end 
   for b = 1:length(OFFT2_11)
   timcs_11(:,b) = timcs_11(:,b)./norm(timcs_11(:,b));    %normalize time course
   end 
    for b = 1:length(OFFT4_12)
   timcs_12(:,b) = timcs_12(:,b)./norm(timcs_12(:,b));    %normalize time course
   end 
  for b = 1:length(OFFT4_10)
   plot(1:1:30, timcs_10(:,b), 'Color', [1 0.3 0.3], 'LineWidth', 0.5)
   hold on;
  end
  hold on;
  for b = 1:length(OFFT2_11)
   plot(1:1:30, timcs_11(:,b), 'Color', [0.3 0.3 1], 'LineWidth', 0.5)
   hold on;
  end
   for b = 1:length(OFFT2_03)
   plot(1:1:30, timcs_03(:,b), 'Color', [0.3 1 0.3], 'LineWidth', 0.5)
   hold on;
   end
  for b = 1:length(OFFT4_12)
   plot(1:1:30, timcs_12(:,b), 'Color', [1 0 1], 'LineWidth', 0.5)
   hold on;
  end

    xlabel('frame number');
    ylabel('contrast');
    title('Time courses - OFFT3 - 2012-11-13, 2012-03-01 (Glaucoma) & 2012-10-15 (Normal)');
    legend('OFFT3- 2012-10-15-0 - Normal');

    
   
   
  




% h = figure;
% plot_time_courses(datarun000_03, ONT2_03, 1, 1, 'figure',-1);
% hold all;
% plot_time_courses(datarun000_11, ONT2_11, 1, 1, 'figure', h);
% 
% plot_time_courses(datarun000_03, ONT2_03, 1, 1)
% hold all;
% plot_time_courses(datarun000_11, ONT2_11, 1, 1)



[259 331 616 1022 1126 1381 1531 2716 3091 3286 4306 4727 4997 6499 6917 6992 7066 7276]; 

plot_rf_summaries(datarun000, ONT1_03,'label', true)
plot_rf_portraits(datarun000_03,OFFT5_03)
plot_time_courses(datarun000, ONT1_03, 1, 1);
plot_autocorrelograms(datarun000_03, OFFT5_03, 'foa', 0);

a = subplot(2,1,1);
plot_time_courses(datarun000_03, ONT2_03, 1, 1, 'figure', a);
title('Time courses - ON T2 -  2012-03-01-0');
b = subplot(2,1,2);
plot_time_courses(datarun000_11, ONT2_11, 1, 1, 'figure', b);
title('Time courses - ON T2 - 2012-11-13-0');

plot_time_courses(datarun000_03, ONT1_03, 1, 1);
hold all;
plot_time_courses(datarun000_11, ONT1_11, 1, 1, 'figure', gcf);

plot_time_courses(datarun000_03, ONT1_03, 1, 1)
plot_time_courses(datarun000_11, ONT1_11, 1, 1)



cellids = intersect(datarun002.cell_ids, datarun000.cell_ids);
[NumSpikesCell, StimComb] = get_spikescellstim(datarun002, datarun002.cell_ids, 0);
NumSpikesCellIn = NumSpikesCell(get_cell_indices(datarun002, cellids), :);
[mag dsindex magmax magave angle rho theta num U V] = dscellanalysis(NumSpikesCell, StimComb);
[magin dsindexin magmaxin magavein anglein rhoin thetain numin Uin Vin] = dscellanalysis(NumSpikesCellIn, StimComb);

plot(magin{1,1}, magin{2,1}, 'o');

NumSpikesCellDS = NumSpikesCell(get_cell_indices(datarun002, A), :);
[magds dsindexds magmaxds magaveds angleds rhods thetads numds] = dscellanalysis(NumSpikesCellDS, StimComb);

plot(angleds{1,1}, angleds{2,1}, '+')


path = '/Analysis/sravi/Rat/Glaucoma/2012-11-13-0/data001-3600-7200s/data002-map/All Plots/';
addpath(path);
sp = 64;
temp = 256;
spatPer = find(unique(StimComb(:,1))==sp);
tempPer = find(unique(StimComb(:,2))==temp);


for aa = 1:length(cellids)
cell_id = cellids(1, aa); 
 
h3 = figure;
%Receptive Field from data000
ax1 = subplot(3,5,[1 2], 'Parent', h3);
plot_rf(datarun001, cell_id, 'foa', ax1, 'scale', 5);
title ('Spatial Receptive Field', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0])
set(ax1,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)

%Time course from data000
ax4 = subplot(3,5,[6 7], 'Parent', h3);
plot_time_course(datarun001,cell_id, 'figure', ax4, 'clear', false); 
title ('Temporal Receptive Field', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax4,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax4, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


%EI from data000
ax5 = subplot(3,5,[11 12], 'Parent', h3);
plot_ei(datarun001, cell_id, 'foa', ax5, 'coordinates' , 'array');
title ('Electrophysiological Image', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax5,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax5, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


%Polar Plot & Spike Raster from data002
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
subpin = [10 5 4 3 8 13 14 15]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
 trigpre = ismember(datarun002.stimulus.trial_list,A(j));
 destaxes=subplot(3,5,subpin(1,j),'Parent', h3);
 spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 0, 8, 0, 10, 'stop', 12, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,5,9)
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
title ('Drifting Grating Rasters', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0], 'Position', [0.16, 5.3, 18.23])


save_figure_pdf(path, num2str(cellids(1,aa)), h3) %uncommnet

% pathname = ['/Users/sravi/matlab/DS cell analysis/2012-10-15-0/AllPlotsInOne-512/'  num2str(cellids(1,aa))];
% set(h3,'PaperOrientation','landscape');
% set(h3,'PaperUnits','normalized');
% set(h3,'PaperPosition', [0 0 1 1]);
% print(h3,'-dpdf', pathname);

close(h3); %uncomment
%set(gcf, 'Color', [1 1 1]);

end


plot_rf_summaries(datarun000_03, ONT1_03,'label', true)
plot_rf_portraits(datarun000_03,ONT1_03, 'scale_factor',5)

plot_rfs(datarun000_03,ONT1_03)


plot_time_courses(datarun000, A, 1, 1);

clear rf
for i = 1:length(ONT1_11)
 rf{i,1} = datarun000_11.stas.rfs{get_cell_indices(datarun000_11, ONT1_11(1,i)),1};
 rf{i,1} =  rf{i,1}./std(std(rf{i,1}));
end

mapp11 = zeros(40,80);
for i = 1:length(ONT1_11)
 mapp11 = mapp11 + rf{i,1};
 imagesc(mapp11); colormap(gray)
 pause;
end
 
clear rf
for i = 1:length(ONT1_03)
 rf{i,1} = datarun000_03.stas.rfs{get_cell_indices(datarun000_03, ONT1_03(1,i)),1};
 rf{i,1} =  rf{i,1}./std(std(rf{i,1}));
end

figure(2)
mapp03 = zeros(40,80);
for i = 1:length(ONT1_03)
 mapp03 = mapp03 + rf{i,1};
 imagesc(mapp03); colormap(gray)
 pause;
end

clear rf
for i = 1:length(ONT1_10)
 rf{i,1} = datarun000_10.stas.rfs{get_cell_indices(datarun000_10, ONT1_10(1,i)),1};
 rf{i,1} =  rf{i,1}./std(std(rf{i,1}));
end

figure(3)
mapp10 = zeros(40,80);
for i = 1:length(ONT1_10)
 mapp10 = mapp10 + rf{i,1};
 imagesc(mapp10); colormap(gray)
 pause;
end


clear rf
for i = 1:length(ONT1_12)
 rf{i,1} = datarun000_12.stas.rfs{get_cell_indices(datarun000_12, ONT1_12(1,i)),1};
 rf{i,1} =  rf{i,1}./std(std(rf{i,1}));
end

figure(4)
mapp12 = zeros(40,80);
for i = 1:length(ONT1_12)
 mapp12 = mapp12 + rf{i,1};
 imagesc(mapp12); colormap(gray)
 pause;
end



figure(5)
subplot(2,2,1)
imagesc(mapp11);colormap(gray);
title('GL 11')
subplot(2,2,2)
imagesc(mapp03);colormap(gray);
title('GL 03')
subplot(2,2,3)
imagesc(mapp12);colormap(gray);
title('GL 12')
subplot(2,2,4)
imagesc(mapp10);colormap(gray);
title('WT 10')

por = plot_rf_portraits(datarun000_11,875);


imagesc(mapp);colormap(gray)

mapp = zeros(40,80);
for i = 1:length(ONT1_11)
 imagesc(rf{i,1}); colormap(gray)
 hold on;
end


[datarun000] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2011-12-23-0/data001-da/', 'data001-map/data001-map', 0, 0, 1);
[datarun002] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2011-12-23-0/data001-da/', 'data004-map/data004-map', 1, '/stimuli/s04', 0);
[datarun001] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2011-12-23-0/data001-da/', 'data005-map/data005-map', 0, 0, 0);

cellids = intersect((intersect(datarun000.cell_ids, datarun001.cell_ids)), datarun002.cell_ids);

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002, datarun002.cell_ids, 0);
NumSpikesCellIn = NumSpikesCell(get_cell_indices(datarun002, cellids), :);
[mag dsindex magmax magave angle rho theta num U V] = dscellanalysis(NumSpikesCell, StimComb);
[magin dsindexin magmaxin magavein anglein rhoin thetain numin Uin Vin] = dscellanalysis(NumSpikesCellIn, StimComb);
save_figure_pdf('/Analysis/sravi/Rat/Glaucoma/2012-03-01-0/DS analysis/', 'All Polar Plots', gcf)

plot(magin{1,1}, magin{2,1}, 'o');

allplotsinone(datarun000, datarun001, datarun002, cellids, rho,theta,U,V,num,StimComb, 64, 256, '/Analysis/sravi/Rat/Glaucoma/2011-12-23-0/data001-da/AllPlots/'); 




[datarun002] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-11-21-1/data001-3600-7200s/', 'data002-map/data002-map', 1, '/stimuli/s02', 0);
[datarun000] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-11-21-1/data001-3600-7200s/', 'data001-map/data001-map', 0, 0, 1);

[datarun001] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2011-12-23-0/data001-da/', 'data005-map/data005-map', 0, 0, 0);

   
 listing = dir('/Analysis/sravi/Rat/Glaucoma/2011-12-23-0/data001-da/DS/')
  

cellids = intersect(datarun000.cell_ids, datarun002.cell_ids);

DS_03 = [51 409 451 452 453 995 1366 1520 1727 1741 1865 2060 2824 3021 3334 4338 4368 4458 4576 4592 4789 5180 5507 5612 5690 5842 5899 5945 5975 6094 6275 6723 6785 6887 7039 7505 7592 7594];
DS_10 = [1595 2042 301 3736 4069 438 4846 5632 6321 708 1683 2449 333 3815 4157 467 4985 5702 6332 7308 1685 257 3512 3842 4173 4731 5148 5719 6751 7475 1097 1895 2898 3636 3934 4353 4789 5150 6229 6797 766]; 
DS_11 = [1115 1896 2044 4384 4532 5252 5342 6619 2747 3437 468 6228 1277 1699 2659 3242 4052 5463 1595 2553 3468 4037 4548 5059 5387 6274 904 1294 1788 2642 3603 4280 4730 5177 5643 6377 920 1368 1956 3167 380 4369 4757 5225 5881 6979 107 1399 2387 3425 3812 4411 4922 5374 5915 736]; 
DS_11_2 = [167 1970 2882 3631 4548 5597 6050 6514 7202 7518 1714 2013 2914 4142 5390 5630 6094 6602 7323 827 1908 2795 2975 4336 5552 5976 6288 6813 7353 1640 1925 2870 3512 4533 557 6046 6468 7054 7412];
DS_12 = [1051 1863 2148 2597 2881 3123 3649 3785 4981 5941 1007 1159 2119 2311 2687 2974 319  3782 4112 5764 7610];

DS_04 = [1321 2432 2476 3348 4517 4531 4533 4712 5161 5209 5238 5899 6125 6393];

plot(magin{5,1}, magin{6,1}, 'o', 'Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('64')
ylabel('256')
save_figure_pdf('/Analysis/sravi/Rat/WildType/2013-04-01-0/DSPlots/Analysis/', 'DS Classf Plot 2', gcf)



ds03 = get_cell_indices(datarun002, DS_04);
plot(magin{3,1}(1,ds03), magin{5,1}(1,ds03), 'o', 'Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
xlabel('64')
ylabel('256')
hold on;
plot(1:0.1:2, 1:0.1:2);

save_figure_pdf('/Analysis/sravi/Rat/WildType/2013-04-01-0/DSPlots/Analysis/', 'DS only 2', gcf)


subplot(1,3,1);
hist(mag{2,1}(1,ds03))
title('32')
subplot(1,3,2);
hist(magin{3,1}(1,ds03))
title('64')
subplot(1,3,3);
hist(magin{5,1}(1,ds03))
title('256')
save_figure_pdf('/Analysis/sravi/Rat/WildType/2013-04-01-0/DSPlots/Analysis/', 'Hist', gcf)



hist((mag{2,1}(1,ds03))./(mag{1,1}(1,ds03)))



 magr_112101 = [1.8412 1.0938 1.1847 1.1869 1.2034 1.0463 1.2055 1.0933 1.7237 1.1317 1.2142 0.9243 1.0549 1.2360 1.3553 1.9415 1.3091 1.1267 0.9455 0.9709 0.9299 1.4389 0.7172 1.2057 1.1045 1.0635 1.0190 1.2834 0.9539 0.7988 0.9281 1.2170 1.1089 1.1651 1.1669 1.6257 1.1727 0.9151 0.7978 1.0443 1.1704 1.1367 1.2086 1.5032 1.0071 1.0705 1.2403 1.3283 1.2616 0.9729 0.7236 1.2057 1.3492 0.9708 1.5293 0.8373];
 magr_040100 = [1.1587 1.1329 1.1256 1.0200 1.1061 1.1665 1.3363 1.3165 0.9898 1.5721 1.1558];
magr_101500 = [1.0377 0.8858 1.0828 1.0364 1.1610 0.9076 1.1168 1.3256 0.8894 1.1076 1.0939 0.8904 0.9687 1.0981 1.0701 0.9505 0.8349 1.0687 0.9788 1.3735 0.9111 1.1638 1.0286 1.0625 0.9349 1.0185 1.9208 0.9669 1.1904 1.6484 1.0522 0.9847 1.1457 1.0572 0.8447 0.9349 1.0511 1.0065 1.0289 1.1130 2.2255];
magr_122300 =  [1.3891 1.1971 1.2223 1.2455 0.8651 0.9776 0.8193 1.4773 1.5001 2.1235 1.5170 1.5338 1.2653 1.4007 1.3688 1.2257 1.3797 0.5489 1.5251 1.3287 0.9885];
magr_111300 = [1.7586 1.1435 0.9509 2.0626 1.0479 1.2988 1.3056 1.2248 1.4137 0.8987 1.2783 0.8685 1.1611 1.7233 1.3545 2.1426 1.7418 1.1743 1.4072 1.1323 1.1434 1.1771 1.0399 0.9610 0.9194 1.2104 0.9827 2.6368 1.1926 1.4479 1.2858 1.3464 1.0461 1.0287 1.0287 1.0814 1.0670 0.9811 1.4063];
magr_030100 = [1.1176 1.4996 0.9237 1.6729 1.2283 1.4449 1.0909 2.7177 1.3558 1.0636 1.0061 1.0605 0.9616 0.9614 1.0088 1.1786 0.9805 1.1747 1.3899 1.0902 1.0543 1.0760 1.3623 1.2130 1.1384 0.8583 1.8047 2.9807 1.8269 3.0769 1.3612 1.8041 1.1429 0.9406 0.9959 1.0688 0.8445 1.4935];



ksdensity(magr_040100(1,1:11));
hold on;
 ksdensity(magr_101500(1,1:41));
hold off;
hold all
ksdensity(magr_112101(1,1:56));
hold on;
 ksdensity(magr_122300(1,1:21));
 ksdensity(magr_111300(1,1:39));
 ksdensity(magr_030100(1,1:38));


% max(magr_122300(1,1:21))/std(magr_122300(1,1:21))
% hold on;
% max(magr_111300(1,1:21))/std(magr_111300(1,1:39))
%  max(magr_030100(1,1:21))/std(magr_030100(1,1:38))
%  max(magr_112101(1,1:21))/std(magr_112101(1,1:56))
% 
% hold all;
%  max(magr_101500(1,1:21))/std(magr_101500(1,1:41))
% max(magr_040100(1,1:21))/std(magr_040100(1,1:11))

 
h(1) = subplot(3,2,1);
ksdensity(magr_112101(1,1:56));
hold on;
 ksdensity(magr_122300(1,1:21));
 ksdensity(magr_111300(1,1:39));
 ksdensity(magr_030100(1,1:38));
 hold off;
 h(2) = subplot(3,2,2);
ksdensity(magr_040100(1,1:11));
hold on;
 ksdensity(magr_101500(1,1:41));
 h(3) = subplot(3,2,3);
ksdensity(magr_112101(1,1:56));
hold on;
 ksdensity(magr_122300(1,1:21));
 ksdensity(magr_111300(1,1:39));
 ksdensity(magr_030100(1,1:38));
 hold all;
ksdensity(magr_040100(1,1:11));
h(4) =subplot(3,2,4);
ksdensity(magr_112101(1,1:56));
hold on;
 ksdensity(magr_122300(1,1:21));
 ksdensity(magr_111300(1,1:39));
 ksdensity(magr_030100(1,1:38));
 hold all;
 ksdensity(magr_101500(1,1:41));
h(5) = subplot(3,2,5);
ksdensity(magr_112101(1,1:56));
hold on;
 ksdensity(magr_122300(1,1:21));
 ksdensity(magr_111300(1,1:39));
 ksdensity(magr_030100(1,1:38));
 hold all;
 ksdensity(magr_040100(1,1:11));
 ksdensity(magr_101500(1,1:41));

linkaxes(h,'y');

ylim(h(1),[0 4.2])

save_figure_pdf('//Analysis/sravi/Rat/', 'KS density', gcf)


subplot(1,3,1);
t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
compass(Uin{2,1}(1,ds03),Vin{2,1}(1,ds03))
title('32')
subplot(1,3,2);
t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
compass(Uin{3,1}(1,ds03),Vin{3,1}(1,ds03))
title('64')
subplot(1,3,3);
t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
compass(Uin{5,1}(1,ds03),Vin{5,1}(1,ds03))
title('256')


save_figure_pdf('/Analysis/sravi/Rat/WildType/2013-04-01-0/DSPlots/Analysis/', 'Polar DS', gcf)

subplot(1,3,1);
t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
compass(Uin{2,1},Vin{2,1})
title('32')
subplot(1,3,2);
t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
compass(Uin{3,1},Vin{3,1})
title('64')
subplot(1,3,3);
t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
compass(Uin{5,1},Vin{5,1})
title('256')



save_figure_pdf('/Analysis/sravi/Rat/WildType/2013-04-01-0/DSPlots/Analysis/', 'Polar All', gcf)


%Direction Tuning
ds = get_cell_indices(datarun002, DS_10);
xi = 0:0.1:6;
for i = 1:length(ds)
    x1 = theta{1,1}(ds(1,i),:);
    y1 = rho{1,1}(ds(1,i),:);
    x2 = theta{2,1}(ds(1,i),:);
    y2 = rho{2,1}(ds(1,i),:);    
    yi1 = interp1(x1,y1,xi, 'pchip'); 
    plot(xi,yi1);
    hold all;
    yi2 = interp1(x2,y2,xi, 'pchip');
    plot(xi,yi2, 'r');
    legend('32', '256');
    pause;
    hold off;
end

for i = 1:length(ds)
    x1 = theta{1,1}(ds(1,i),:);
    y1 = rho{1,1}(ds(1,i),:);
    yi1 = interp1(x1,y1,xi, 'pchip'); 
    plot(xi,yi1);
    hold on;
    pause;
end


for i = 1:length(ds)
    x1 = theta{1,1}(ds(1,i),:);
    y1 = rho{1,1}(ds(1,i),:);
    x2 = theta{2,1}(ds(1,i),:);
    y2 = rho{2,1}(ds(1,i),:);    
    yi1 = interp1(x1,y1,xi, 'pchip'); 
    plot(xi,yi1);
    hold all;
    yi2 = interp1(x2,y2,xi, 'pchip');
    plot(xi,yi2, 'r');
    legend('32', '256');
    pause;
    hold off;
end



misespara = cell(length(ds),1);
response = cell(length(ds),1);
for i = 1:length(ds)
A0 = 1;
u0 = angle{1,1}(1,ds(1,i));
if(u0 < 0)
    u0 = u0 + (2*pi);
end
kappa0 = 1;
func1 = @ComputeMises;
misespara{i,1} = fminsearch(@(params)sseerror(func1, theta{1,1}(ds(1,i),:),rho{1,1}(ds(1,i),:), params), [A0,u0,kappa0],optimset('TolX',1e-10,'MaxFunEvals',50000,'MaxIter',12000)); 
response{i,1} = ComputeMises(misespara{i,1}, 0:0.1:6);
plot(theta{1,1}(ds(1,i),:), rho{1,1}(ds(1,i),:), 'o')
hold all;
plot(0:0.1:6, response{i,1});
hold off;
pause;
params_64_1015(i,1) = max(response{i,1}) - min(response{i,1}); %A
params_64_1015(i,2) = misespara{i,1}(1,3); %std, k
end


plot(params_256_1113(:,1), params_256_1113(:,2), 'o')
hold all;
plot(params_32_1113(:,1), params_32_1113(:,2), 'o')
xlabel('A');
ylabel('k');

ksdensity(params_256_1113(:,2))
hold all;
ksdensity(params_32_1113(:,2))

plot(params_256_0301(:,1), params_256_0301(:,2), 'o')
hold all;
plot(params_32_0301(:,1), params_32_0301(:,2), 'o')
xlabel('A');
ylabel('k');

ksdensity(params_256_0301(:,2))
hold all;
ksdensity(params_32_0301(:,2))

plot(params_256_1223(:,1), params_256_1223(:,2), 'o')
hold all;
plot(params_32_1223(:,1), params_32_1223(:,2), 'o')
xlabel('A');
ylabel('k');

ksdensity(params_256_1223(:,2))
hold all;
ksdensity(params_32_1223(:,2))



plot(params_256_1121(:,1), params_256_1121(:,2), 'o')
hold all;
plot(params_32_1121(:,1), params_32_1121(:,2), 'o')
xlabel('A');
ylabel('k');

ksdensity(params_256_1121(:,2))
hold all;
ksdensity(params_32_1121(:,2))


plot(params_256_1015(:,1), params_256_1015(:,2), 'o',  'Color', [0 0 0], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 0], 'MarkerFaceColor' , [0.9 0.8 0.1])
hold all;
plot(params_64_1015(:,1), params_64_1015(:,2), 'o')
xlabel('A');
ylabel('k');

ksdensity(params_256_1015(:,2))
hold all;
ksdensity(params_64_1015(:,2))

plot(params_256_0401(:,1), params_256_0401(:,2), 'o',  'Color', [0 0 0], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 0], 'MarkerFaceColor' , [0.9 0.3 0.1])
hold all;
plot(params_32_0401(:,1), params_32_0401(:,2), 'o')
plot(params_64_0401(:,1), params_64_0401(:,2), 'o')

xlabel('A');
ylabel('k');

ksdensity(params_256_0401(:,2))
hold all;
ksdensity(params_64_0401(:,2))
ksdensity(params_32_0401(:,2))



params_256_0401(:,1) = [0.8551 0.5589 0.4666 1.0297 0.9920 0.8657 0.9472 1.0367 0.7987 0.8815 0.8681 0.8867 0.3474 0.8105];
params_256_0401(:,2) = [ 0.8003 0.5107 0.3737 1.8015 1.5415 0.7971 1.2213 1.4148 0.7056 0.8541 1.1849 0.9789 0.2867 0.9278];
   
params_32_0401(:,1) = [0.8408 0.3281 0.5684 0.9808 0.7996 0.9282 0.8750 0.9891 0.6240 0.6815 0.8063 0.5768 0.3215 0.6776];
params_32_0401(:,2) = [0.9763 0.3139 0.4819 1.7873 0.7602 1.0395 1.0644 1.4170 0.4792 0.6463 0.7568 0.4525 0.2240 0.5467];

params_64_0401(:,1) = [ 1.0377 0.5821 0.6346 1.0011 0.9870 1.0312 0.9479 0.9584 0.6909 0.8795 0.8724 0.8325 0.1712 0.7424];
params_64_0401(:,2) = [1.5381 0.8004 0.5465 1.6835 1.3277 1.5981 1.1552 1.3668 0.5632 1.1984 1.1264 0.7727 0.1308 0.6775];


params_256_1113(:,1), 


y(1,1) = mean(params_256_1113(:,2))
x(1,1) = mean(params_256_1113(:,1))
e(1,1) = std(params_256_1113(:,2))
e1(1,1) = std(params_256_1113(:,1))

y(1,2) = mean(params_256_0301(:,2))
x(1,2) = mean(params_256_0301(:,1))
e(1,2) = std(params_256_0301(:,2))
e1(1,2) = std(params_256_0301(:,1))


y(1,3) = mean(params_256_1121(:,2))
x(1,3) = mean(params_256_1121(:,1))
e(1,3) = std(params_256_1121(:,2))
e1(1,3) = std(params_256_1121(:,1))


y(1,4) = mean(params_256_1223(:,2))
x(1,4) = mean(params_256_1223(:,1))
e(1,4) = std(params_256_1223(:,2))
e1(1,4) = std(params_256_1223(:,1))

y(1,5) = mean(params_256_0401(:,2))
x(1,5) = mean(params_256_0401(:,1))
e(1,5) = std(params_256_0401(:,2))
e1(1,5) = std(params_256_0401(:,1))


y(1,6) = mean(params_256_1015(:,2))
x(1,6) = mean(params_256_1015(:,1))
e(1,6) = std(params_256_1015(:,2))
e1(1,6) = std(params_256_1015(:,1))

plot(x,y,'o')
errorbar(x(1,1:4),y(1,1:4),e(1,1:4), 'o')
hold on;
errorbar(x(1,5:6),y(1,5:6),e(1,5:6), 'or')
xlabel('A')
ylabel('k')
title('Plot of mean of amplitude and tuning parameters with error bar for tuning width')

plot(y,x,'o')
errorbar(y(1,1:4),x(1,1:4),e1(1,1:4), 'o')
hold on;
errorbar(y(1,5:6),x(1,5:6),e1(1,5:6), 'or')
xlabel('k')
ylabel('A')
title('Plot of mean of amplitude and tuning parameters with error bar for amplitude')

errorbar(mean(x(1,1:4)),mean(y(1,1:4)),std(y(1,1:4)), 'o')
hold on;
errorbar(mean(x(1,5:6)),mean(y(1,5:6)),std(y(1,5:6)), 'or')
xlabel('A')
ylabel('k')
title('Plot of mean of amplitude and tuning parameters with error bar for tuning width for glaucoma vs WT')

for i = 1:length(ds)
    a = []; b = [];
    a(1,1) = U{1,1}(1,ds(1,i));
    a(1,2) = V{1,1}(1,ds(1,i));
    b(1,1) = U{2,1}(1,ds(1,i));
    b(1,2) = V{2,1}(1,ds(1,i));
    dp(1,i) = dot(a,b);
end

 dp_0401 = [2.6918 0.4754 0.8267 3.2140 2.5351 3.0013 2.9681 3.4797 1.9451 2.4903 2.4646 1.8430 0.5448 2.1330];
%11 out of 141 cells
%dp_1015 - 40 out of 354 cells are ds

 dp_1113 = [1.2660 2.4013 1.1023 0.8896 2.4939 1.2552 1.1370 2.4121 2.6476 2.3555 2.3897 0.8461 1.5570 0.9462 1.3087 1.2883 1.0620 2.1080 1.0860 1.6663 2.5422 1.4331 1.3060 1.8578 1.2392 1.1368 1.9870 1.1029 1.7301 1.1282 1.4943 1.4493 1.1165 1.4108 1.0526 1.0742 2.2047 2.5839 1.1695];
dp_0301 = [1.7242 0.8151 0.5329 0.7213 0.6324 0.2281 2.0269 0.9578 0.8234 2.5735 1.6039 0.8256 1.3677 1.1013 1.1295 2.1386 1.8849 0.7380 0.7458 1.5725 0.9451 1.1997 2.1277 0.8311 2.1534 0.6019 1.1072 0.4531 0.5184 0.3998 0.7015 1.0469 0.9885 1.4371 1.3352 2.0614 0.7421 0.4658];
dp_1223 = [2.5232 2.4898 1.5566 1.7262 0.8828 0.9732 1.9532 0.9253 1.0603 0.5075 1.1429 1.2043 2.2043 0.8129 1.4110 2.5496 1.0679 1.1856 2.5803 0.9536 -0.1536];
dp_1121 = [1.7909 2.3849 1.9055 2.9639 3.0998 4.0066 2.7364 3.8154 2.0146 3.3306 1.5310 2.8770 2.7461 1.6919 2.4724 1.9893 2.2737 1.7343 2.9803 3.5416 2.0338 1.6416 1.1837 2.4371 1.2859 1.6674 2.0063 1.6904 2.9518 1.4467 3.6701 2.0235 1.7620 1.2872 2.4580 2.0350 2.5993 2.0726 1.0752 3.1891 2.7376 1.9375 1.9660 2.1580 3.4763 1.9137 3.6980 1.1087 1.4051 1.2096 1.8116 1.6963 1.8201 3.0271 1.7984 1.9905];
%56 / 354 DS cells
 
 
 
 
 