[datarun006] = load_dsdata('/Analysis/sravi/Mouse/2013-04-09-0/', 'data006/data006', 1, '/stimuli/s06', 0);
[datarun011] = load_dsdata('/Analysis/sravi/Mouse/2013-04-09-0/', 'data011/data011', 1, '/stimuli/s11', 0);
[datarun022] = load_dsdata('/Analysis/sravi/Mouse/2013-04-09-0/', 'data022/data022', 1, '/stimuli/s22', 0);
 
[NumSpikesCell6, StimComb6] = get_spikescellstim(datarun006, datarun006.cell_ids, 0);
[NumSpikesCell11, StimComb11] = get_spikescellstim(datarun011, datarun011.cell_ids, 0);
[NumSpikesCell22, StimComb22] = get_spikescellstim(datarun022, datarun022.cell_ids, 0);
 
 
[magin6 dsindexin6 magmaxin6 magavein6 anglein6 rhoin6 thetain6 numin6 Uin6 Vin6] = dscellanalysis(NumSpikesCell6, StimComb6);
[magin11 dsindexin11 magmaxin11 magavein11 anglein11 rhoin11 thetain11 numin11 Uin11 Vin11] = dscellanalysis(NumSpikesCell11, StimComb11);
[magin22 dsindexin22 magmaxin22 magavein22 anglein22 rhoin22 thetain22 numin22 Uin22 Vin22] = dscellanalysis(NumSpikesCell22, StimComb22);

plot(magin22{3,1}, magin22{4,1},'o','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
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
cellids = datarun022.cell_ids(1,chos);

[datarun022006ind datarun006ind] = mapwithei(datarun022, datarun006, cellids, 0.89);
[datarun022011ind datarun011ind] = mapwithei(datarun022, datarun011, datarun022.cell_ids(1, datarun022006ind), 0.82);



c6 = datarun006ind;
c11 = datarun011ind;
c22 = datarun022006ind;


% figure(1)
% UU22{1,1} = Uin22{4,1}(1,c22);
% VV22{1,1} = Vin22{4,1}(1,c22);
% magmag22{1,1} = magin22{4,1}(1,c22);
% polar_plots_all(UU22, VV22, 64, magmag22)
% 
% figure(2)
% UU11{1,1} = Uin11{4,1}(1,c11);
% VV11{1,1} = Vin11{4,1}(1,c11);
% magmag11{1,1} = magin11{4,1}(1,c11);
% polar_plots_all(UU11, VV11, 64, magmag11)
% 
% figure(3)
% UU6{1,1} = Uin6{4,1}(1,c6);
% VV6{1,1} = Vin6{4,1}(1,c6);
% magmag6{1,1} = magin6{4,1}(1,c6);
% polar_plots_all(UU6, VV6, 64, magmag6)
 












% all_rasterplots_pdf(NumSpikesCell22,StimComb22,datarun022, {rhoin22{3,1}},{thetain22{3,1}},{Uin22{3,1}},{Vin22{3,1}},numin22(3,1), c22(1,1), 0,1);
% all_rasterplots_pdf(NumSpikesCell11,StimComb11,datarun011, {rhoin11{3,1}},{thetain11{3,1}},{Uin11{3,1}},{Vin11{3,1}},numin11(3,1), c11(1,1), 0,1)
% all_rasterplots_pdf(NumSpikesCell6,StimComb6,datarun006, {rhoin6{3,1}},{thetain6{3,1}},{Uin6{3,1}},{Vin6{3,1}},numin6(3,1), c6(1,1), 0,1)
% %  raster_plot(StimComb22,datarun022,,numin22(3,1),T,R,Unew, Vnew, c22(1,1));

%Make raster plots for all cells

for k = 1:length(c22)
[T22 R22 Unew22 Vnew22] = polar_plots_one({rhoin22{6,1}},{thetain22{6,1}},{Uin22{6,1}},{Vin22{6,1}},numin22(6,1),c22(1,k));

all = 1:length(datarun022.stimulus.combinations);

[an indic] = sortrows(StimComb22, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [12 3 2 1 10 19 20 21]; %Places on subplot for each angle from 0 to 315 deg
h3 = figure;
in1 = ismember(an(:,1), 64)'; %
SC2 = an(ismember(an(:,1), 64), :);%
A1 = all(in1);%

i = 1;
in = ismember(SC2(:,2), numin22(6,1))';
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun022.stimulus.trial_list,A(j));
    destaxes=subplot(3,9,subpin(i,j), 'Parent', h3);
    spikesbytrials = get_raster(datarun022.spikes{c22(1,k),1}, datarun022.stimulus.triggers(trigpre), 0, 8, 0, 4, 'stop', 10, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,9,11, 'Parent', h3)
polar(T22{i,1}, R22{i,1});
hold on;
h1 = compass(Unew22(i,1),Vnew22(i,1), 'r'); %Vector average plot
set(h1,'linewidth',3) 
hold off;

figure;
[T11 R11 Unew11 Vnew11] = polar_plots_one({rhoin11{6,1}},{thetain11{6,1}},{Uin11{6,1}},{Vin11{6,1}},numin11(6,1),c11(1,k));
figure(2)
all = 1:length(datarun011.stimulus.combinations);

[an indic] = sortrows(StimComb11, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [15 6 5 4 13 22 23 24]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), 64)'; %
SC2 = an(ismember(an(:,1), 64), :);%
A1 = all(in1);%

i = 1;
in = ismember(SC2(:,2), numin11(6,1))';
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun011.stimulus.trial_list,A(j));
    destaxes=subplot(3,9,subpin(i,j), 'Parent', h3);
    spikesbytrials = get_raster(datarun011.spikes{c11(1,k),1}, datarun011.stimulus.triggers(trigpre), 0, 8, 0, 4, 'stop', 10, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,9,14, 'Parent', h3)
polar(T11{i,1}, R11{i,1});
hold on;
h1 = compass(Unew11(i,1),Vnew11(i,1), 'r'); %Vector average plot
set(h1,'linewidth',3) 
hold off;

figure;
[T6 R6 Unew6 Vnew6] = polar_plots_one({rhoin6{6,1}},{thetain6{6,1}},{Uin6{6,1}},{Vin6{6,1}},numin6(6,1),c6(1,k));
figure(2)
all = 1:length(datarun006.stimulus.combinations);

[an indic] = sortrows(StimComb6, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [18 9 8 7 16 25 26 27]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), 64)'; %
SC2 = an(ismember(an(:,1), 64), :);%
A1 = all(in1);%

i = 1;
in = ismember(SC2(:,2), numin6(6,1))';
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun006.stimulus.trial_list,A(j));
    destaxes=subplot(3,9,subpin(i,j), 'Parent', h3);
    spikesbytrials = get_raster(datarun006.spikes{c6(1,k),1}, datarun006.stimulus.triggers(trigpre), 0, 8, 0, 4, 'stop', 10, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,9,17, 'Parent', h3)
polar(T6{i,1}, R6{i,1});
hold on;
h1 = compass(Unew6(i,1),Vnew6(i,1), 'r'); %Vector average plot
set(h1,'linewidth',3) 
hold off;

save_figure_pdf('/Analysis/sravi/Mouse/2013-04-09-0/DS Analysis/64-256/', num2str(datarun022.cell_ids(1,(c22(1,k)))),h3);

close all;

end

%Angle and vector average tuning subplots


close all;
for i = 1:length(c22)
    xxv = 0:1:512;
    va = zeros(3, length(magin22));
    ang = zeros(3, length(magin22));
    for j = 1:length(magin22)
        va(1,j) = magin22{j,1}(1, c22(1,i));
        va(2,j) = magin11{j,1}(1, c11(1,i));
        va(3,j) = magin6{j,1}(1, c6(1,i));
        ang(1,j) = anglein22{j,1}(1, c22(1,i));
        ang(2,j) = anglein11{j,1}(1, c11(1,i));
        ang(3,j) = anglein6{j,1}(1, c6(1,i));
    end
    subplot(3,2,1)
    plot(numin22, va(1,:), 'b', 'LineWidth',1.5, 'Marker', 'o', 'MarkerFaceColor', 'w')
    hold on;
    plot(numin22, va(2,:), 'g', 'LineWidth',1.5, 'Marker', 'o',  'MarkerFaceColor', 'w')
    plot(numin22, va(3,:), 'r', 'LineWidth',1.5, 'Marker', 'o',  'MarkerFaceColor', 'w')
    legend('NDF 0', 'NDF 3', 'NDF 4');
    title(num2str(datarun022.cell_ids(1,c22(1,i))));
    xlabel('Speed')
    ylabel('Vector Average')
    hold off;
    subplot(3,2,2)
    plot(numin22, ang(1,:), 'b', 'LineWidth',1.5, 'Marker', 'o',  'MarkerFaceColor', 'w')
    hold on;
    plot(numin22, ang(2,:), 'g', 'LineWidth',1.5, 'Marker', 'o',  'MarkerFaceColor', 'w')
    plot(numin22, ang(3,:), 'r', 'LineWidth',1.5, 'Marker', 'o',  'MarkerFaceColor', 'w')
    legend('NDF 0', 'NDF 3', 'NDF 4');
    xlabel('Speed')
    ylabel('Vector Angle')
    hold off;
    
    yyv1 = interp1(numin22,va(1,:),xxv, 'pchip');
    yyv2 = interp1(numin22,va(2,:),xxv, 'pchip');
    yyv3 = interp1(numin22,va(3,:),xxv, 'pchip');
   
    subplot(3,2,3)
    plot(numin22,va(1,:),'o', xxv, yyv1, 'b');
    hold on;
    plot(numin22,va(2,:), 'o', xxv, yyv2, 'g');
    plot(numin22,va(3,:), 'o', xxv, yyv3, 'r');
    hold off;
    
    yya1 = interp1(numin22,ang(1,:),xxv, 'pchip');
    yya2 = interp1(numin22,ang(2,:),xxv, 'pchip');
    yya3 = interp1(numin22,ang(3,:),xxv, 'pchip');

    subplot(3,2,4)
    plot(numin22,ang(1,:),'o', xxv, yya1, 'b');
    hold on;
    plot(numin22,ang(2,:), 'o', xxv, yya2, 'g');
    plot(numin22,ang(3,:), 'o', xxv, yya3, 'r');
    hold off;
    
    subplot(3,2,5)
    plot(xxv, yya1.*yyv1, 'b', 'LineWidth', 1.2);
    hold on;
    plot(xxv, yya2.*yyv2, 'g', 'LineWidth', 1.2);
    plot(xxv, yya3.*yyv3, 'r', 'LineWidth', 1.2);
    hold off;
    xlabel('Speed')
    ylabel('Vector Mag * Vector Angle')
%     pause;
%     save_figure_pdf('/Analysis/sravi/Mouse/2013-04-09-0/DS Analysis/VA plots/', num2str(datarun022.cell_ids(1,c22(1,i))), gcf);
    poptunv1(i,:) = yyv1;
    poptuna1(i,:) = yya1;
    poptunv2(i,:) = yyv2;
    poptuna2(i,:) = yya2;
    poptunv3(i,:) = yyv3;
    poptuna3(i,:) = yya3;
end

subplot(2,3,1)
plot(xxv, poptunv1, 'b');
xlabel('Speed')
ylabel('Vector Average')
title('Population tuning')
subplot(2,3,2)
plot(xxv, poptunv2, 'g');
subplot(2,3,3)
plot(xxv, poptunv3, 'r');
subplot(2,3,4)
plot(xxv, poptuna1, 'b');
xlabel('Speed')
ylabel('Vector Angle')
subplot(2,3,5)
plot(xxv, poptuna2, 'g');
subplot(2,3,6)
plot(xxv, poptuna3, 'r');



'/Analysis/sravi/Mouse/2013-04-09-0/DS Analysis/64-256/', num2str(datarun022.cell_ids(1,(c22(1,k)))),h3









