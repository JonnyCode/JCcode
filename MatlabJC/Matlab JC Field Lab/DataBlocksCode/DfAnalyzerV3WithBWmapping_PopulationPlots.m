
% This script will combine analyses from many DB from 'DfAnalyzerV2WithBWmapping'
% JC 2017

%% populations

DBsetArray{1} = [3,23,24] ; % 1 m rescue
PlotColorArray{1} = 'r' ;

DBsetArray{2} = [19,20,21,30,31] ; % 3m KO
PlotColorArray{2} = 'b' ;

DBsetArray{3} = [26,28,29] ; % 1m KO 
PlotColorArray{3} = 'c' ;

DBsetArray{4} = [12,25,27] ; % CNG deltaCam
PlotColorArray{4} = 'k' ;

% DBsetArray{5} = [32,33,34] ; % ELFN1 mice
% PlotColorArray{5} = 'g' ;

%% Within DBsetArray stats

% % make histograms cumulative, normalize, and calculate fraction no threshold (BAD CODING: this should all be done in DfAnalyzerV3WithBWmapping)
% for a=1:length(DBsetArray) ; % make arrays and find pop statistics
%     DBset = DBsetArray{a} ;
%     
%     for DB = DBset ;
%         % calculate fraction no threshold
%         ForMatlab(DB).Threshold_FracNoThresh = ForMatlab(DB).Threshold_HistValues(end)/sum(ForMatlab(DB).Threshold_HistValues) ;
%         if ~isnan(ForMatlab(DB).ExampleCellTypei) % cell type average threshold 
%             ForMatlab(DB).Threshold_byCellType_FracNoThresh{ForMatlab(DB).ExampleCellTypei} = ...
%                 ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei}(end)/...
%                 sum(ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei}) ;
%         end
%     end
% 
%         % make cumulative and get rid of no thresh cells
%     for DB = DBset ;
%         ForMatlab(DB).Threshold_HistValues = cumsum(ForMatlab(DB).Threshold_HistValues(1:end-1)) ;
%         ForMatlab(DB).Bg_rate_HistValues = cumsum(ForMatlab(DB).Bg_rate_HistValues(1:end-1)) ;
%         if ~isnan(ForMatlab(DB).ExampleCellTypei) % cell type average threshold 
%             ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei} = cumsum(ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei}(1:end-1)) ;
%         end
%         ForMatlab(DB).Threshold_HistBins = ForMatlab(DB).Threshold_HistBins(1:end-1) ;
%     end
% 
%         % normalize to max
%     for DB = DBset ;
%         ForMatlab(DB).Threshold_HistValues = ForMatlab(DB).Threshold_HistValues/ForMatlab(DB).Threshold_HistValues(end) ;
%         ForMatlab(DB).Bg_rate_HistValues = ForMatlab(DB).Bg_rate_HistValues/ForMatlab(DB).Bg_rate_HistValues(end) ;
%         if ~isnan(ForMatlab(DB).ExampleCellTypei) % cell type average threshold 
%             ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei} =...
%                 ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei}/...
%             ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei}(end) ;
%         end
%     end
% end

for a=1:length(DBsetArray) ; % make arrays and find pop statistics
    DBset = DBsetArray{a} ;
    st = 0 ;
    Threshold_HistValues_array= [];
    Threshold_FracNoThresh_array= [];
    Bg_rate_HistValues_array= [] ;
    Threshold_byCellType_HistValues_array = [] ;
    
    for DB = DBset ;
        st=st+1 ;
        Threshold_HistValues_array(st,:) = ForMatlab(DB).Threshold_HistValues ;
        Threshold_FracNoThresh_array(st,:) = ForMatlab(DB).Threshold_FracNoThresh ;
        Bg_rate_HistValues_array(st,:) = ForMatlab(DB).Bg_rate_HistValues ;
%         if ~isnan(ForMatlab(DB).ExampleCellTypei) % cell type average threshold 
%             Threshold_byCellType_HistValues_array(st,:) = ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei} ;
%         else
%             Threshold_byCellType_HistValues_array(st,:) = nan(size(Threshold_HistValues_array(st,:))) ;
%         end
    end
    
    Threshold_HistValues_AllCellAverage(a,:) = nanmean(Threshold_HistValues_array,1) ;
    Threshold_FracNoThresh_AllCellAverage(a,:) = nanmean(Threshold_FracNoThresh_array,1) ;
    Bg_rate_HistValues_AllCellAverage(a,:) = nanmean(Bg_rate_HistValues_array,1) ;
    %Threshold_byCellType_HistValues_AllCellAverage(a,:) = nanmean(Threshold_byCellType_HistValues_array,1) ;
    
    Threshold_HistValues_AllCellSem(a,:) = nanstd(Threshold_HistValues_array,[],1)/sqrt(length(DBset)) ;
    Threshold_FracNoThresh_AllCellSem(a,:) = nanstd(Threshold_FracNoThresh_array,[],1)/sqrt(length(DBset)) ;
    Bg_rate_HistValues_AllCellSem(a,:) = nanstd(Bg_rate_HistValues_array,[],1)/sqrt(length(DBset)) ;
    %Threshold_byCellType_HistValues_AllCellSem(a,:) = nanstd(Threshold_byCellType_HistValues_array,[],1)/sqrt(sum(~isnan(Threshold_byCellType_HistValues_array(:,1)))) ;
end
  
% number of cells
for a=1:length(DBsetArray) ; % make arrays and find pop statistics
    DBset = DBsetArray{a} ;
    N{a} = [] ;
    for DB = DBset ;
        N{a} = [N{a},ForMatlab(DB).HistCellNumber] ;
    end 
end
        
    
    

%% Between DBsetArray stats

% KS test between all cell threshold means
for a=1:length(DBsetArray) ; % each DBsetArray
    for a2=1:length(DBsetArray) ; % each other DBsetArray
        [temp,KS(a,a2)] = kstest2(Threshold_HistValues_AllCellAverage(a,:),Threshold_HistValues_AllCellAverage(a2,:)) ;
    end
end

% Ttest between fraction of cells with no threshold
Temp = cell(1,length(DBsetArray)) ; % prep array

for a=1:length(DBsetArray) ; % each DBsetArray
    for DB = DBsetArray{a} ; 
        Temp{a} = [Temp{a},ForMatlab(DB).Threshold_FracNoThresh] ;
    end
end
 
for a=1:length(DBsetArray) ; % each DBsetArray
    for a2=1:length(DBsetArray) ; % each other DBsetArray
        [temp,NoThreshFractTest(a,a2)] = ttest2(Temp{a},Temp{a2}) ;
    end
end
        
% difference in means of individual threshold distributions
Temp = cell(1,length(DBsetArray)) ; % prep array

for a=1:length(DBsetArray) ; % each DBsetArray
    for DB = DBsetArray{a} ; 
        Temp{a} = [Temp{a},sum([0,diff(ForMatlab(DB).Threshold_HistValues)].*ForMatlab(DB).Threshold_HistBins)]; % mean
    end
    ThresholdMeanAllCellAv(a) = 10^mean(Temp{a}) ; % means
    ThresholdMeanAllCellSem(a) = 10^(std(Temp{a})/sqrt(length(DBsetArray{a}))) ;
end
        
for a=1:length(DBsetArray) ; % each DBsetArray
    for a2=1:length(DBsetArray) ; % each DBsetArray 
        [temp,meantTest(a,a2)] = ttest2(Temp{a},Temp{a2});
    end  
end

% difference in medians of individual threshold distributions
Temp = cell(1,length(DBsetArray)) ; % prep array

for a=1:length(DBsetArray) ; % each DBsetArray
    for DB = DBsetArray{a} ; 
        i1=find(ForMatlab(DB).Threshold_HistValues<0.5,1,'last'); % interpolate for median
        Temp{a} = [Temp{a},(ForMatlab(DB).Threshold_HistBins(i1)+ForMatlab(DB).Threshold_HistBins(i1+1))/2] ;
    end
    ThresholdMedianAllCellAv(a) = 10^mean(Temp{a}) ; % means
    ThresholdMedianAllCellSem(a) = 10^(std(Temp{a})/sqrt(length(DBsetArray{a}))) ;
end
        
for a=1:length(DBsetArray) ; % each DBsetArray
    for a2=1:length(DBsetArray) ; % each DBsetArray 
        [temp,mediantTest(a,a2)] = ttest2(Temp{a},Temp{a2});
    end  
end

%% for Igor
ForIgor = struct() ;

for a=1:length(DBsetArray) ;
    VecName = ['ThreshHistValuesMean','DbSet',num2str(a)] ;
    ForIgor = setfield(ForIgor,VecName,Threshold_HistValues_AllCellAverage(a,:)) ; 

    VecName = ['ThreshHistValuesSem','DbSet',num2str(a)] ;
    ForIgor = setfield(ForIgor,VecName,Threshold_HistValues_AllCellSem(a,:)) ;
    
    VecName = ['FracNoThreshMean','DbSet',num2str(a)] ;
    ForIgor = setfield(ForIgor,VecName,Threshold_FracNoThresh_AllCellAverage(a,:)) ; 

    VecName = ['FracNoThreshSem','DbSet',num2str(a)] ;
    ForIgor = setfield(ForIgor,VecName,Threshold_FracNoThresh_AllCellSem(a,:)) ;
end

%% figures

figure % all cell histograms
for a=1:length(DBsetArray) ; % plots
    DBset = DBsetArray{a} ;
    PlotColor = PlotColorArray{a}
    
    for DB = DBset ; 
        plot(ForMatlab(DB).Threshold_HistBins,ForMatlab(DB).Threshold_HistValues,[PlotColor,'-']) ;
        hold on
        %errorbar(ForMatlab(DB).Threshold_HistBins,Threshold_HistValues_AllCellAverage(a,:),...
            %Threshold_HistValues_AllCellSem(a,:),[PlotColor,'-'],'LineWidth',1) ;
    end
    ylabel('fraction of cells')
    xlabel('theshold')
end

figure % all cell histograms (no threshold cells)
for a=1:length(DBsetArray) ; % plots
    DBset = DBsetArray{a} ;
    PlotColor = PlotColorArray{a} ;
    
    for DB = DBset ; 
        plot(a,ForMatlab(DB).Threshold_FracNoThresh,[PlotColor,'*'])
        hold on
        plot(a,Threshold_FracNoThresh_AllCellAverage(a),[PlotColor,'+']) 
    end
    ylabel('fraction cells no threshold')
end
xlim([0, a+1])    

figure(3) % all cell histograms bg rates
for a=1:length(DBsetArray) ; % plots
    DBset = DBsetArray{a} ;
    PlotColor = PlotColorArray{a} ;

    for DB = DBset ; 
        DBset = DBsetArray{a} ;
        PlotColor = PlotColorArray{a} ;

        plot(ForMatlab(DB).Bg_rate_HistBins(1:end-1),...
            cumsum(ForMatlab(DB).Bg_rate_HistValues(1:end-1))/...
            sum(ForMatlab(DB).Bg_rate_HistValues(1:end-1)),[PlotColor,'-']) ;
        hold on
        plot(ForMatlab(DB).Bg_rate_HistBins(1:end-1),...
            cumsum(BgRate_HistValues_AllCellAverage(1:end-1))/...
            sum(BgRate_HistValues_AllCellAverage(1:end-1)),[PlotColor,'-'],'LineWidth',5) ;
    end
    ylabel('fraction of cells')
    xlabel('Bg Rate')
end
    



figure % example cell type histograms
for a=1:length(DBsetArray) ;
    DBset = DBsetArray{a} ;
    PlotColor = PlotColorArray{a}
    
    for DB = DBset ; 
        if ~isnan(ForMatlab(DB).ExampleCellTypei)
            plot(ForMatlab(DB).Threshold_HistBins,...
                ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei},...
                [PlotColor,'-']) ;
            hold on
        %     plot(ForMatlab(DB).Threshold_HistBins(1:end-1),...
        %         cumsum(Threshold_HistValues_byCellTypeAverage(1:end-1))/...
        %         sum(Threshold_HistValues_byCellTypeAverage(1:end-1)),...
        %         [PlotColor,'-'],'LineWidth',2) ;
        end
    end
end
    
figure % example cell type histograms (no threshold cells)
for a=1:length(DBsetArray) ;
DBset = DBsetArray{a} ;
PlotColor = PlotColorArray{a}
    
    for DB = DBset ; 
        if ~isnan(ForMatlab(DB).ExampleCellTypei)
            plot(ForMatlab(DB).Threshold_HistBins(1),...
                ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei}(end)/...
                sum(ForMatlab(DB).Threshold_byCellType_HistValues{ForMatlab(DB).ExampleCellTypei}(1:end)),...
                [PlotColor,'*'])
            hold on
        %     plot(ForMatlab(DB).Threshold_HistBins(1),...
        %         Threshold_HistValues_byCellTypeAverage(end),...
        %         [PlotColor,'+']) 
        end
    end
end

figure % example cell tuning curve
for a=1:length(DBsetArray) ;
DBset = DBsetArray{a} ;
PlotColor = PlotColorArray{a}
    for DB = DBset ; 
        errorbar(ForMatlab(DB).Irel_unique,...
            ForMatlab(DB).psth_peak_delta_trls_mean(:,ForMatlab(DB).ExampleCelli),...
            ForMatlab(DB).psth_peak_delta_trls_std(:,ForMatlab(DB).ExampleCelli),...
            PlotColor) 
        hold on
    end
    set(gca,'Xscale','log')
end

 % example cell psth
for a=1:length(DBsetArray) ;
DBset = DBsetArray{a} ;
PlotColor = PlotColorArray{a}
    for DB = DBset ;
        figure
        for ts = 1:length(ForMatlab(DB).psth_delta) ; % for each trigger set
            plot(ForMatlab(DB).psthTime,...
                ForMatlab(DB).psth_delta{ts}(ForMatlab(DB).ExampleCelli,:),...
                PlotColor) 
            hold on
        end
    end
end





