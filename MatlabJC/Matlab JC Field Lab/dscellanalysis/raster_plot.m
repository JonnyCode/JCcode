function[] = raster_plot(StimComb,datarun,num,T,R,Unew, Vnew,cellnum)
%Assumes 8 directions, if more, change subpin

all = 1:length(datarun.stimulus.combinations);

[an indic] = sortrows(StimComb, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [6 3 2 1 4 7 8 9]; %Places on subplot for each angle from 0 to 315 deg
h3 = figure;
in1 = ismember(an(:,1), 64)'; %
SC2 = an(ismember(an(:,1), 64), :);%
A1 = all(in1);%

i = 1;
in = ismember(SC2(:,2), num(i))';
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun.stimulus.trial_list,A(j));
    destaxes=subplot(3,3,subpin(i,j));
    spikesbytrials = get_raster(datarun.spikes{cellnum,1}, datarun.stimulus.triggers(trigpre), 0, 8, 0, 4, 'stop', 10, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,3,5)
polar(T{i,1}, R{i,1});
hold on;
h1 = compass(Unew(i,1),Vnew(i,1), 'r'); %Vector average plot
set(h1,'linewidth',3) 
hold off;
% h=gcf;
% a = get(gca,'children');
end