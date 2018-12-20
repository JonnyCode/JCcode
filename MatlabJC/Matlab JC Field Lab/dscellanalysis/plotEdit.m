%DS plots


plot(magin{1,1}, magin{2,1},'o','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
title ('Normalized vector average in preferred direction of cell', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
hXLabel = xlabel('Drifitng grating with fast temporal period - 32', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
hYLabel = ylabel('Drifting grating with slow temporal period - 256', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);   
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/My Publications/Rat Classification 2014/2013-05-31-1/DS/', '1', gcf)



plot(magin{1,1}, magin{2,1},'o');
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

UU{1,1} = Uin{2,1}(1,in);
VV{1,1} = Vin{2,1}(1,in);
magmag{1,1} = magin{2,1}(1,in);
polar_plots_all(UU, VV, 256, magmag)
set(gcf, 'Color', [1 1 1]);
hTitle  = title ('Preferred Direction of Direction Selective Cells');
set(hTitle, 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
%PUT MARKERS BEFORE SAVING
save_figure_pdf('/Users/sravi/Documents/My Documents/My Publications/Rat Classification 2014/2013-05-31-1/DS/', '2', gcf)

%Angle plots
plot(anglein{1,1}(1, chos(1,:)), anglein{2,1}(1, chos(1,:)), '+')
ginput(1);
Y = get(get(gca,'Children'),'YData');
X = get(get(gca,'Children'),'XData');
h = impoly;
accepted_pos = wait(h);
hold on;
drawPolygon(accepted_pos(:,1), accepted_pos(:,2));
in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));
plot(X(in),Y(in),'r+');
cidcid = 1:length(X);
choschos = ds(in);


plot_rf_summaries(datarun000, chos90, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/My Publications/Rat Classification 2014/2013-05-31-1/DS/', '3-90degrees', gcf)

temp_indices = get_cell_indices(datarun000, chos270);
color_palette = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 0; 1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5; 1 0.5 0; 1 0 0.5; 0.5 1 0; 0.5 0 1; 0 0.5 1; 0 1 0.5; 0.25 0.5 1; 0.25 1 0.5; 1 0.25 0.5; 1 0.5 0.25; 0.5 1 0.25; 0.5 0.25 1; 0 0.8 0.4; 0 0.4 0.8; 0.8 0 0.4; 0.8 0.4 0; 0.4 0 0.8];
 
for cc = 1:length(temp_indices)
%    plot_ei(datarun, zero_degs(1))
 
    temp_ei = datarun000.ei.eis{temp_indices(cc)};
    
    plot_ei_(temp_ei, datarun000.ei.position, 0, 'cutoff', 0.05, 'elec_colors', repmat(color_palette(cc,:), 512, 1), 'alpha', 0.5)
 
end
axis off
print(1, '/Users/sravi/Documents/My Documents/My Publications/Rat Classification 2014/2013-05-31-1/DS/4-EI270deg.pdf', '-dpdf')






[mag dsindex magmax magave angle rho theta num U V] = dscellanalysis(NumSpikesCell, StimComb);

sp = 64;
temp = 256;
spatPer = find(unique(StimComb(:,1))==sp);
tempPer = find(unique(StimComb(:,2))==temp);

cell_id = 438;   
    
%Receptive Field from data000
plot_rf(datarun000, cell_id, 'scale', 5);
title ('Spatial Receptive Field', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/DS/', 'rfoo', gcf)

%Time course from data000
plot_time_course(datarun000,cell_id); %'clear', false needed at varargin param?
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Temporal Receptive Field', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/DS/', 'tcoo', gcf)


%EI from data000
plot_ei(datarun000, cell_id, 'coordinates' , 'array');
title ('Electrophysiological Image', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('microns', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('microns', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/DS/', 'eioo', gcf)


%ACF from data000
datarun000 = get_autocorrelations(datarun000, 'all');
plot_autocorrelograms(datarun000, cell_id, 'foa', 0)
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Autocorrelation Function', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('Time Difference (s)', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('# of pairs', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/DS/', 'acfoo', gcf)


h3 = figure;
p = plot(datarun000.autocorrelation{get_cell_indices(datarun000, 7533),1}.bins, datarun000.autocorrelation{get_cell_indices(datarun000, 7533),1}.probabilities);
set(p,'Color','black','LineWidth',2)



figure;
[T R Unew Vnew] = polar_plots_one(rho(:,spatPer),theta(:,spatPer),U(:,spatPer),V(:,spatPer),num,get_cell_indices(datarun002, cell_id));
close;
h3 = figure;
all = 1:length(datarun002.stimulus.combinations);
[an indic] = sortrows(StimComb, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [6 3 2 1 4 7 8 9]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(3,3,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 0, 8, 0, 8, 'stop', 12, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,3,5)
p = polar(T{tempPer,1}, R{tempPer,1});
th = findall(gca,'Type','text');
for i = 1:length(th),
set(th(i),'FontName', 'AvantGarde', 'FontSize', 14)
end
h = findall(gca, 'type', 'line'); % find all lines in the current figure
h(h==p) = []; % remove the line you ploted from the list.
set(h, 'LineStyle', '--');
hold on;
h1 = compass(Unew(tempPer,1),Vnew(tempPer,1), 'r'); %Vector average plot
set(h1,'linewidth',2) 
hold off;
title ('Drifting Grating Rasters and Polar Plot', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0], 'Position' , [0, 4.8, 18.23]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/DS/', 'dgoo', gcf)

%Pulse Raster
h3 = figure;
wh = datarun001.triggers(21:4:length(datarun001.triggers), 1); %For 2012-10-10-1 dataset, that is how the triggers are arranges
gr = datarun001.triggers(22:4:length(datarun001.triggers),1); %Might change with dataset
[hp ap] = pulse_analysis(datarun001, get_cell_indices(datarun001, 7533), 0, path, wh, gr); %change raster plot axes

pulse_analysis(datarun001, get_cell_indices(datarun001, cell_id), 0, path); %change raster plot axes
title ('Full Field Flashing Pulse Raster', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]); 
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
xlabel('seconds', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('trials', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/DS/', 'pulsesoo', gcf)

%Static Nonlinearity
datarun000 = load_java_movie(datarun000,'/Volumes/lab/acquisition/movie-xml/BW-10-1-0.48-11111-60x60-60.35.xml');
datarun000 = get_snls(datarun000, 7533);
%  x = -1:0.001:1;
%  y = exp(-3.2795+(4.2972*x)); %e(-b+ag)
% plot(x,y)

[X, Y] = curve_from_binning(datarun000.stas.snls{418,1}.gen_signal, datarun000.stas.snls{418,1}.spikes, 'average_y', 'mean','average_x', 'mean', 'num_bins', 20);
Y(:, 1) = Y(:, 1)*60.35/1;
plot(X,Y, 'b') 
title ('Static Non Linearity', 'FontName','Tahoma', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]); 
set(gca, 'FontName', 'Helvetica', 'FontSize', 16);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [.005 .005],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
xlabel('Generator Signal', 'FontName','AvantGarde', 'FontSize', 16, 'Color', [0 0 0]);
ylabel('Spikes', 'FontName','AvantGarde', 'FontSize', 16, 'Color', [0 0 0]);
save_figure_pdf('/Users/sravi/Documents/My Documents/Grants/Sep 2013/Single cell figures/', 'SNL-7533-1', gcf)

%plot_snl_(datarun000.stas.snls{9,1}.gen_signal,datarun000.stas.snls{9,1}.spikes,'foa',0,'fit',datarun000.stas.snls{9,1}.fit_params)





%off
%type 1

[COEFF,SCORE] = princomp(tcnormnorm(:,b)');
plot(SCORE(:,1), SCORE(:,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
title ('Principal Components of Normalized Time Courses for OFF Cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
xlabel('Principal Component 1', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('Principal Component 2', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'classfplot', gcf)

%mosaics
plot_rf_summaries(datarun000, c, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'mosaict1', gcf)


%tcs
plot_time_courses(datarun000,c, 1, 1);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'tcallt1', gcf)



cell_id = 3933;
%Receptive Field from data000
plot_rf(datarun000, cell_id, 'scale' ,5);
title ('Spatial Receptive Field', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'rft1', gcf)

%Time course from data000
plot_time_course(datarun000,cell_id); %'clear', false needed at varargin param?
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Temporal Receptive Field', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'tct1', gcf)


%EI from data000
plot_ei(datarun000, cell_id, 'coordinates' , 'array');
title ('Electrophysiological Image', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('microns', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('microns', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'eit1', gcf)


%ACF from data000
plot_autocorrelograms(datarun000, cell_id, 'foa', 0)
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Autocorrelation Function', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('Time Difference (s)', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('# of pairs', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'acft1', gcf)

figure;
[T R Unew Vnew] = polar_plots_one(rho(:,spatPer),theta(:,spatPer),U(:,spatPer),V(:,spatPer),num,get_cell_indices(datarun002, cell_id));
close;
h3 = figure;
all = 1:length(datarun002.stimulus.combinations);
[an indic] = sortrows(StimComb, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [6 3 2 1 4 7 8 9]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(3,3,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 0, 8, 0, 8, 'stop', 12, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,3,5)
p = polar(T{tempPer,1}, R{tempPer,1});
th = findall(gca,'Type','text');
for i = 1:length(th),
set(th(i),'FontName', 'AvantGarde', 'FontSize', 14)
end
h = findall(gca, 'type', 'line'); % find all lines in the current figure
h(h==p) = []; % remove the line you ploted from the list.
set(h, 'LineStyle', '--');
hold on;
h1 = compass(Unew(tempPer,1),Vnew(tempPer,1), 'r'); %Vector average plot
set(h1,'linewidth',2) 
hold off;
title ('Drifting Grating Rasters and Polar Plot', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0], 'Position' , [0, 4.8, 18.23]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'dgt1', gcf)



pulse_analysis(datarun001, get_cell_indices(datarun001, cell_id), 0, path); %change raster plot axes
title ('Full Field Raster', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]); 
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
xlabel('seconds', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('trials', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'pulsest1', gcf)



%mosaics
plot_rf_summaries(datarun000, c, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'mosaict2', gcf)


%tcs
plot_time_courses(datarun000,c, 1, 1);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'tcallt1', gcf)


cell_id = 6196;
%Receptive Field from data000
plot_rf(datarun000, cell_id, 'scale', 5);
title ('Spatial Receptive Field', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'rft2', gcf)

%Time course from data000
plot_time_course(datarun000,cell_id); %'clear', false needed at varargin param?
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Temporal Receptive Field', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'tct2', gcf)


%EI from data000
plot_ei(datarun000, cell_id, 'coordinates' , 'array');
title ('Electrophysiological Image', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('microns', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('microns', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'eit2', gcf)


%ACF from data000
plot_autocorrelograms(datarun000, cell_id, 'foa', 0)
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Autocorrelation Function', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('Time Difference (s)', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('# of pairs', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'acft2', gcf)

figure;
[T R Unew Vnew] = polar_plots_one(rho(:,spatPer),theta(:,spatPer),U(:,spatPer),V(:,spatPer),num,get_cell_indices(datarun002, cell_id));
close;
h3 = figure;
all = 1:length(datarun002.stimulus.combinations);
[an indic] = sortrows(StimComb, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [6 3 2 1 4 7 8 9]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(3,3,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 0, 8, 0, 8, 'stop', 12, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,3,5)
p = polar(T{tempPer,1}, R{tempPer,1});
th = findall(gca,'Type','text');
for i = 1:length(th),
set(th(i),'FontName', 'AvantGarde', 'FontSize', 14)
end
h = findall(gca, 'type', 'line'); % find all lines in the current figure
h(h==p) = []; % remove the line you ploted from the list.
set(h, 'LineStyle', '--');
hold on;
h1 = compass(Unew(tempPer,1),Vnew(tempPer,1), 'r'); %Vector average plot
set(h1,'linewidth',2) 
hold off;
title ('Drifting Grating Rasters and Polar Plot', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0], 'Position' , [0, 4.8, 18.23]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT1T2/', 'dgt2', gcf)





%off type 3
[COEFF,SCORE] = princomp(acnormmx(:,b)');
plot(SCORE(:,1), SCORE(:,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
title ('PCs of Normalized Autocorrelation Functions for OFF Cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
xlabel('Principal Component 1', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('Principal Component 2', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT3/', 'classfplot', gcf)


%mosaics
plot_rf_summaries(datarun000, c, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT3/', 'mosaict3', gcf)


%tcs
plot_time_courses(datarun000,c, 1, 1);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT3/', 'tcallt3', gcf)



cell_id = 3286;
%Receptive Field from data000
plot_rf(datarun000, cell_id, 'scale', 5);
title ('Spatial Receptive Field', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/UNC/', 'rft5', gcf)

%Time course from data000
plot_time_course(datarun000,cell_id); %'clear', false needed at varargin param?
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Temporal Receptive Field', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/UNC/', 'tct5', gcf)


%EI from data000
plot_ei(datarun000, cell_id, 'coordinates' , 'array');
title ('Electrophysiological Image', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('microns', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('microns', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT3/', 'eit3', gcf)


%ACF from data000
plot_autocorrelograms(datarun000, cell_id, 'foa', 0)
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Autocorrelation Function', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('Time Difference (s)', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('# of pairs', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT3/', 'acft3', gcf)

figure;
[T R Unew Vnew] = polar_plots_one(rho(:,spatPer),theta(:,spatPer),U(:,spatPer),V(:,spatPer),num,get_cell_indices(datarun002, cell_id));
close;
h3 = figure;
all = 1:length(datarun002.stimulus.combinations);
[an indic] = sortrows(StimComb, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [6 3 2 1 4 7 8 9]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(3,3,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 0, 8, 0, 8, 'stop', 12, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,3,5)
p = polar(T{tempPer,1}, R{tempPer,1});
th = findall(gca,'Type','text');
for i = 1:length(th),
set(th(i),'FontName', 'AvantGarde', 'FontSize', 14)
end
h = findall(gca, 'type', 'line'); % find all lines in the current figure
h(h==p) = []; % remove the line you ploted from the list.
set(h, 'LineStyle', '--');
hold on;
h1 = compass(Unew(tempPer,1),Vnew(tempPer,1), 'r'); %Vector average plot
set(h1,'linewidth',2) 
hold off;
title ('Drifting Grating Rasters and Polar Plot', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0], 'Position' , [0, 4.8, 18.23]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/UNC/', 'dgt5', gcf)



pulse_analysis(datarun001, get_cell_indices(datarun001, cell_id), 0, path); %change raster plot axes
title ('Full Field Raster', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]); 
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
xlabel('seconds', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('trials', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/UNC/', 'pulsest5', gcf)


c = [3635 3889 7040 76 5672 407 2521 1427 1564 2206 1130 2373 4021 4097 4697 1411 843 6992 901 2732 2794 3857 4279 6139 6376 7517 226 2253 2555 3091 4732 6257 6380 1246 3695 4771 5116 6260 6886 2326 4235];

c = [182 528 857 1532 1591 3586 4231 6155 6439 6737 2328 3812 4130 6125 6363 6392 6752 649 1084 1578 2236 2716 5705 7069 424 496 2311 3152 3241 3530 3859 4145 4501 6152 6304 7186]; 



[COEFF,SCORE] = princomp(tcnormauc(:,b)');
plot(SCORE(:,1), SCORE(:,2),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
title ('Principal Components of Normalized Time Courses for ON Cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
xlabel('Principal Component 1', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('Principal Component 2', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/ONT1T2T3/', 'classfplot', gcf)

c = [62 991 1156 4234 4278 4487 5733 6286 6931];
%mosaics
plot_rf_summaries(datarun000, c, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT3/', 'rfalloff-2', gcf)


%tcs
plot_time_courses(datarun000,c, 1, 1);
h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 1.2);
title ('Temporal Receptive Fields of all cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT3/', 'tcalloff-2', gcf)



cell_id = 3452;
%2nd plot
acnormmx,tcnormmx

[COEFF1,SCORE1] = princomp(acnormmx(:,b)');
[COEFF2,SCORE2] = princomp(tcnormmx(:,b)');

plot(SCORE1(:,1)./(max(SCORE1(:,1))), SCORE2(:,1)./(max(SCORE2(:,1))),'s','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
title ('Principal Components of Normalized Time Courses for OFF Cells', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
xlabel('Principal Component 1', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('Principal Component 2', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/OFFT3/', 'classfplot2', gcf)


pulsenormsdPSTH , tcnormmx

c = [889 1684 2478 2896 3200 3648 3935 4188 4294 5645 7278];

plot_rf_summaries(datarun000, c, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/ARVO/2013/UNC/', 'mosaicon', gcf)






% set(h1,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
% 
% 
% 
%  
%   , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'YGrid'       , 'on'      , ...
%   'XColor'      , [.3 .3 .3], ...
%   'YColor'      , [.3 .3 .3], ...
%   'YTick'       , 0:500:2500, ...
%   'LineWidth'   , 1         );
% 
% 
% close;
% 
% 
% 
% 
% set(hE                            , ...
%   'LineStyle'       , 'none'      , ...
%   'Marker'          , '.'         , ...
%   'Color'           , [.3 .3 .3]  );
% set(hData                         , ...
%   'LineStyle'       , 'none'      , ...
%   'Marker'          , '.'         );
% set(hModel                        , ...
%   'LineStyle'       , '--'        , ...
%   'Color'           , 'r'         );
% set(hCI(1)                        , ...
%   'LineStyle'       , '-.'        , ...
%   'Color'           , [0 .5 0]    );
% set(hCI(2)                        , ...
%   'LineStyle'       , '-.'        , ...
%   'Color'           , [0 .5 0]    );
% 
% 
% 
% set(hFit                          , ...
%   'LineWidth'       , 2           );
% set(hE                            , ...
%   'LineWidth'       , 1           , ...
%   'Marker'          , 'o'         , ...
%   'MarkerSize'      , 6           , ...
%   'MarkerEdgeColor' , [.2 .2 .2]  , ...
%   'MarkerFaceColor' , [.7 .7 .7]  );
% set(hData                         , ...
%   'Marker'          , 'o'         , ...
%   'MarkerSize'      , 5           , ...
%   'MarkerEdgeColor' , 'none'      , ...
%   'MarkerFaceColor' , [.75 .75 1] );
% set(hModel                        , ...
%   'LineWidth'       , 1.5         );
% set(hCI(1)                        , ...
%   'LineWidth'       , 1.5         );
% set(hCI(2)                        , ...
%   'LineWidth'       , 1.5         );