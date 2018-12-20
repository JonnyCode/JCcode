function[] = allplotsinone(datarun000, datarun001, datarun002, cellids, rho,theta,U,V,num,StimComb, sp, temp, path, wh, gr) 

addpath(path);
spatPer = find(unique(StimComb(:,1))==sp);
tempPer = find(unique(StimComb(:,2))==temp);

%Code to analyze 2012-10-15-0 dataset

%Load all 3 stimulus runs (2 hour white noise, flashing full field pulses,
% %drifitng gratings)
% [datarun000] = load_dsdata();
% 
% [datarun001] = load_dsdata();
% 
% [datarun002] = load_dsdata();
% 
% %Find cell IDs that exist in all 3 datasets
% cellids = intersect(intersect(datarun000.cell_ids, datarun001.cell_ids),datarun002.cell_ids);

%Put spatial receptive field, time course, EI, polar+raster plot, pulse
%plot of those cells in a figure and save each cell as one pdf file

for aa = 1:length(cellids)
cell_id = cellids(1, aa);    
    
h3 = figure;
%Receptive Field from data000
ax1 = subplot(3,7,[1 2], 'Parent', h3);
plot_rf(datarun000, cell_id, 'foa', ax1, 'scale', 5);
title ('Spatial RF', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0])
set(ax1,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)

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
subpin = [12 5 4 3 10 17 18 19]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(3,7,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 'axis_range', [0 8 0 10], 'stop', 10, 'foa', destaxes, 'tic_color', [0 0 0]);
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
% title ('DG Rasters', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0], 'Position', [0.16, 5.3, 18.23])

% all_rasterplots_pdf(NumSpikesCell,StimComb,datarun002, rho,theta,U,V,num, get_cell_indices(datarun002, cell_id), 0); %change raster plot axes
% f_c = openfig('one.fig');
%     set(gcf, 'Visible', 'off')
%     ax2 = gca;
%     figure(h3)
%    s2 = subplot(3,3,[2 5 8], 'Parent', h3);
%     fig2 = get(ax2,'children');
%     copyobj(fig2,s2);
% 
% 
% fig2 = get(gca,'children');
% findall(gcf)
% 
% get(gca,'Children')
% allchild(gca)
% 
% copyobj(allchild(gcf),ax2);   
% 
% 
% set(gcf, 'Visible', 'off')
% s2 = gca;
% figure(h3)
% 
% fig2 = get(s2,'children');


% >> saveas(h,'AllPlot','jpg') 
% hgsave(gcf,'one.fig') 
% axes_to_be_copied = findobj(f_c,'type','axes'); 
% chilred_to_be_copied = get(axes_to_be_copied,'children'); 
% 
% copyobj(chilred_to_be_copied,ax2);
% fig2 = get(gca,'children');
% copyobj(fig2,ax2);
% copyobj(allchild(ax1),ax2);



%Pulse Raster from data001
ax3 = subplot(3,7,[6 7 13 14 20 21], 'Parent', h3);
[hp ap] = pulse_analysis(datarun001, get_cell_indices(datarun001, cell_id), 0, path, wh, gr,10,true,0.1); %change raster plot axes
fig3 = get(ap,'children');
copyobj(fig3,ax3);
close(hp);
title ('W-G-B-G', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]); 
%set(ax3, 'FontName', 'Helvetica', 'FontSize', 12);
set(ax3,'Box', 'off', 'TickDir','out', 'TickLength', [.005 .005],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
xlabel('seconds', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
%ylabel('trials', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


%Time course from data000
ax4 = subplot(3,7,[8 9], 'Parent', h3);
plot_time_course(datarun000,cell_id, 'figure', ax4, 'clear', false); %'clear', false needed at varargin param?
title ('Temporal RF', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax4,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax4, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);

% 
% % %EI from data000
% ax5 = subplot(3,7,[15 16], 'Parent', h3);
% plot_ei(datarun000, cell_id, 'foa', ax5, 'coordinates' , 'array');
% title ('Electrophysiological Image', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
% set(ax5,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
% set(ax5, 'FontName', 'Helvetica', 'FontSize', 12);
% xlabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
% ylabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);

% ACF from data000
% ax5 = subplot(3,7,[15 16], 'Parent', h3);
% plot(datarun000.autocorrelation{get_cell_indices(datarun000, cell_id),1}.bins, datarun000.autocorrelation{get_cell_indices(datarun000, cell_id),1}.probabilities);
% title ('ACF', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
% set(ax5,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
% set(ax5, 'FontName', 'Helvetica', 'FontSize', 12);
% xlabel('bins', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
% ylabel('probability', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);

%ISI from data000
ax5 = subplot(3,7,[15 16], 'Parent', h3);
plot(datarun000.interspikeinterval{get_cell_indices(datarun000, cell_id),1}.bins, datarun000.interspikeinterval{get_cell_indices(datarun000, cell_id),1}.probabilities);
title ('ACF', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax5,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax5, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('bins', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('probability', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);





save_figure_pdf(path, num2str(cellids(1,aa)), h3) %uncommnet

% pathname = ['/Users/sravi/matlab/DS cell analysis/2012-10-15-0/AllPlotsInOne-512/'  num2str(cellids(1,aa))];
% set(h3,'PaperOrientation','landscape');
% set(h3,'PaperUnits','normalized');
% set(h3,'PaperPosition', [0 0 1 1]);
% print(h3,'-dpdf', pathname);

close(h3); %uncomment
%set(gcf, 'Color', [1 1 1]);
end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 


