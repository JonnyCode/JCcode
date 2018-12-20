function[h, a, spikesbytrials, sumSpTrTrig, nhist] = pulse_analysis(datarun, chos, save, path, whiteT, grayAWh, tEnd, P,binWidth)

%Prints out spikes in response to Pulse Data for cells
%Input datarun structure with pulse stimuli
%Input cells (chos) to print out spikes for
%Input 1 or 0 for save if you want to save file to pdf
%Treat each W,G,B,G as a trial

addpath(path);
spikesbytrials = cell(length(chos), 1);
sumSpTrTrig = cell(length(chos), 1);
nhist = cell(length(chos), 1);
%nhistt = cell(length(chos), 1);
%xouthist = cell(length(chos), 1);
binSize = 0:binWidth:tEnd;


% whiteT = datarun.triggers(1:4:length(datarun.triggers), 1); %For 2012-10-10-1 dataset, that is how the triggers are arranges
% grayAWh = datarun.triggers(2:4:length(datarun.triggers),1); %Might change with dataset
% grayAWh = datarun.triggers(3:4:length(datarun.triggers),1); %Might change with dataset

numtrials = min(length(whiteT), length(grayAWh));

if(length(whiteT) > numtrials)
    whiteT(length(whiteT), :) = []; %Make them have the same number of triggers / trials
end

if(length(grayAWh) > numtrials)
    grayAWh(length(grayAWh), :) = [];
end

b = numtrials+1; %Light steps: white, gray, black
g = b+0.5;
w = g+0.5; 

for i = 1:length(chos)
    pathname = [path  num2str(chos(1,i)) '-' num2str(datarun.cell_ids(chos(1,i)))];
    spikesbytrials{i,1} = get_raster(datarun.spikes{chos(1,i), 1}, whiteT, 'tic_color', [0 0 0], 'plot', P, 'axis_range', [0, tEnd, 0, w+3]); %change depending on whether you want plot
    hold on
    if (tEnd == 10)
        stairs([0 3 5 8 10],[w g b g g], 'Color', 'r', 'LineWidth',1);
    elseif (tEnd == 12)
            stairs([0 3 6 9 12],[w g b g g], 'Color', 'r', 'LineWidth',1);
    elseif (tEnd == 20)        
            stairs([0 5 10 15 20],[w g b g g], 'Color', 'r', 'LineWidth',1);
    elseif (tEnd == 24)        
            stairs([0 3 12 15 24],[w g b g g], 'Color', 'r', 'LineWidth',1);
    end

    for j = 1:length(whiteT)
        if (tEnd == 10)
            sumSpTrTrig{i,1}(j,1) = sum(spikesbytrials{i,1}{j,1} >= 0 & spikesbytrials{i,1}{j,1} < 3);
            sumSpTrTrig{i,1}(j,2) = sum(spikesbytrials{i,1}{j,1} >= 3 & spikesbytrials{i,1}{j,1} < 5);
            sumSpTrTrig{i,1}(j,3) = sum(spikesbytrials{i,1}{j,1} >= 5 & spikesbytrials{i,1}{j,1} < 8);
            sumSpTrTrig{i,1}(j,4) = sum(spikesbytrials{i,1}{j,1} >= 8 & spikesbytrials{i,1}{j,1} < 10);
        elseif (tEnd == 12)
            sumSpTrTrig{i,1}(j,1) = sum(spikesbytrials{i,1}{j,1} >= 0 & spikesbytrials{i,1}{j,1} < 3);
            sumSpTrTrig{i,1}(j,2) = sum(spikesbytrials{i,1}{j,1} >= 3 & spikesbytrials{i,1}{j,1} < 6);
            sumSpTrTrig{i,1}(j,3) = sum(spikesbytrials{i,1}{j,1} >= 6 & spikesbytrials{i,1}{j,1} < 9);
            sumSpTrTrig{i,1}(j,4) = sum(spikesbytrials{i,1}{j,1} >= 9 & spikesbytrials{i,1}{j,1} < 12);
        elseif (tEnd == 20)
            sumSpTrTrig{i,1}(j,1) = sum(spikesbytrials{i,1}{j,1} >= 0 & spikesbytrials{i,1}{j,1} < 5);
            sumSpTrTrig{i,1}(j,2) = sum(spikesbytrials{i,1}{j,1} >= 5 & spikesbytrials{i,1}{j,1} < 10);
            sumSpTrTrig{i,1}(j,3) = sum(spikesbytrials{i,1}{j,1} >= 10 & spikesbytrials{i,1}{j,1} < 15);
            sumSpTrTrig{i,1}(j,4) = sum(spikesbytrials{i,1}{j,1} >= 15 & spikesbytrials{i,1}{j,1} < 20);
         elseif (tEnd == 24)
            sumSpTrTrig{i,1}(j,1) = sum(spikesbytrials{i,1}{j,1} >= 0 & spikesbytrials{i,1}{j,1} < 3);
            sumSpTrTrig{i,1}(j,2) = sum(spikesbytrials{i,1}{j,1} >= 3 & spikesbytrials{i,1}{j,1} < 12);
            sumSpTrTrig{i,1}(j,3) = sum(spikesbytrials{i,1}{j,1} >= 12 & spikesbytrials{i,1}{j,1} < 15);
            sumSpTrTrig{i,1}(j,4) = sum(spikesbytrials{i,1}{j,1} >= 15 & spikesbytrials{i,1}{j,1} < 24);
        end
    for k = 1:(length(binSize)-1)
        nhist{i,1}(j,k) = sum(spikesbytrials{i,1}{j,1} >= binSize(1,k) & spikesbytrials{i,1}{j,1} < binSize(1,k+1));
    end
    
    %[nhistt{i,1}(j,:) xouthist{i,1}(j,:)] = hist(spikesbytrials{i,1}{j,1},100); %gives different bins every time
    end
    psth = sum(nhist{i,1})/numtrials;
    psth = psth./max(psth);
    psth = psth.*2;
    stairs(binWidth:binWidth:tEnd,psth+w+1);
        hold off
        h = gcf;
        a = gca;
        %close;
     if(save)
        set(h,'PaperOrientation','landscape');
        set(h,'PaperUnits','normalized');
        set(h,'PaperPosition', [0 0 1 1]);
        print(h,'-dpdf', pathname);
    end

end

%Spikes by trials has spikes during each trial for each cell (eg - 27
%cells, each having 50 trials)




%close;
end


