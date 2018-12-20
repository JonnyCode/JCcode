function[] = raster_plots(StimComb,datarun,num,T,R,Unew, Vnew,cellnum)
%Assumes 8 directions, if more, change subpin

all = 1:length(datarun.stimulus.combinations);

[an indic] = sortrows(StimComb, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [9 3 2 1 7 13 14 15; 12 6 5 4 10 16 17 18]; %Places on subplot for each angle from 0 to 315 deg
h3 = figure;
in1 = ismember(an(:,1), 64)'; %
SC2 = an(ismember(an(:,1), 64), :);%
A1 = all(in1);%


for i = 1:length(num)
    in = ismember(SC2(:,2), num(i))';
    A = A1(in);
    for  j = 1:sum(in)
        trigpre = ismember(datarun.stimulus.trial_list,A(j));
        destaxes=subplot(3,6,subpin(i,j));
        %spikesbytrials = get_raster(datarun.spikes{cellnum,1}, datarun.stimulus.triggers(trigpre),'stop', mean(diff(datarun.stimulus.triggers)), 'foa', destaxes, 'tic_color', [0 0 0]);
                spikesbytrials = get_raster(datarun.spikes{cellnum,1}, datarun.stimulus.triggers(trigpre), 'axis_range', [0, 8, 0, 10], 'stop', 10, 'foa', destaxes, 'tic_color', [0 0 0]);
        %hold on
   %        rasterphli(datarun, datarun.cell_ids(1,cellnum), datarun.stimulus.triggers(trigpre), 'stop',8, 'histax' , destaxes, 'hist', true, 'hist_bin', 0.2)
  %get_psth(datarun.spikes{cellnum,1},  datarun.stimulus.triggers(trigpre),'stop', 8, 'plot_hist', true, 'bin_size', 1, 'foa', destaxes);
       %get_psth(datarun.spikes{cellnum,1},  datarun.stimulus.triggers(trigpre),'stop', 8, 'plot_hist', true, 'bin_size', 2, 'foa', destaxes);

    end
    if(i == 1)
          subplot(3,6,8)
          polar(T{i,1}, R{i,1});
          hold on;
          h1 = compass(Unew(i,1),Vnew(i,1), 'r'); %Vector average plot
          set(h1,'linewidth',3) 
          hold off;
    elseif(i == 2)
            subplot(3,6,11)
            polar(T{i,1}, R{i,1});
            hold on;
            h1 = compass(Unew(i,1),Vnew(i,1), 'r'); %Vector average plot
            set(h1,'linewidth',3) 
            hold off;
    end
end

end