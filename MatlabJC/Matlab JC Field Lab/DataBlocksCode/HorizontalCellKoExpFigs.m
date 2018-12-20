
% DataBlocks (DB)

% DataBlock(11).DgConcat = [RootAnalysisDir,'2015-11-12-0/data003-7-11/data003-7-11'] ; % wash 100% contrast
% DataBlock(12).DgConcat = [RootAnalysisDir,'2015-11-17-0/data000-2-5_RS/data000-2-5_RS'] ; % wash 12,24,48% contrast
% DataBlock(13).DgConcat = [RootAnalysisDir,'2016-01-25-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
% DataBlock(14).DgConcat = [RootAnalysisDir,'2016-02-03-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
% DataBlock(16).DgConcat = [RootAnalysisDir,'2016-03-09-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
% DataBlock(17).DgConcat = [RootAnalysisDir,'2016-03-10-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
% DataBlock(18).DgConcat = [RootAnalysisDir,'2016-03-14-0/data000-2-4/data000-2-4'] ; % wash 12,24,48% contrast
% DataBlock(19).DgConcat = [RootAnalysisDir,'2016-03-16-0/data002-4-6/data002-4-6'] ; % 12,24,48% contrast
% DataBlock(20).DgConcat = [RootAnalysisDir,'2016-03-18-0/data002-5-7/data002-5-7'] ; % 12,24,48% contrast
% DataBlock(21).DgConcat = [RootAnalysisDir,'2016-03-21-0/data002-4-6/data002-4-6'] ; % 12,24,48% contrast
% DataBlock(23).DgConcat = [RootAnalysisDir,'2016-04-20-0/data002-4-6/data002-4-6'] ; % 12,24,48% contrast

% (19) 2016-03-16 - no gfp seen on array (lost tissue in IHC), shows similar effect to November 
% (12) 2015-11-17 - sparse gfp most clear effects
% (21) 2016-03-21 - very weak/sparse gfp shows some effects similar to November
% (16) 2016-03-09 - good transfection, no clear effect seen yet
% (17) 2016-03-10 - good transfection, no clear effect seen yet
% (20) 2016-03-18 - ok transfection, no clear effect seen yet
% (13) 2016-01-25 - CBA control
% (14) 2016-02-03 - CBA control
% (18) 2016-03-14 - very good transfection, no clear effect yet
% (15) 2016-02-05 - CBA control
% (11) 2016-11-12 - sparse gfp, no clear effect (long experiment) 

MatFilePath = '/Volumes/lab/temp/cx57_psam_analysis/MatFiles' ; % path for mat files
DB = 13 ;

load([MatFilePath,'/KoSpatialV4',num2str(DB)])

% FIGURES
   
% psths of each cell in control, +psem, wash
cells=1 ;
while cells <=l_cells ; % for each cell
    figure(1)
    set(gcf,'name',[num2str(cells),' LfMse W ',num2str(ceil(MseRelWashLfMean(cells)*100)/100),...
                    ' D ',num2str(ceil(MseRelDrugLfMean(cells)*100)/100),...
                    ' D/W ',num2str(ceil(MseRelDrug_DivMseRelWashLfMean(cells)*100)/100)])
    clf
    for dset=1:3 ;
        for sp=1:l_sp ; % for each spatial period 
            for cntrst=1:l_cntrst ; % for each contrast
                subplot(l_sp,l_cntrst,l_cntrst*(sp-1)+cntrst) ;
                plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cells,:),'color',Color_list{dset})
                axis tight
                hold on
                title(['Mse W ',num2str(ceil(MseRelWash{cells}(cntrst,sp)*100)/100),...
                    ' D ',num2str(ceil(MseRelDrug{cells}(cntrst,sp)*100)/100),...
                    ' D/W ',num2str(ceil(MseRelDrug_DivMseRelWash{cells}(cntrst,sp)*100)/100)]) ;
            end
        end
    end
    
    figure(2)
    clf
    % master cell indicy
    masteri=[] ;
    for clm=1:length(cell_list_map) ;
        if cell_list_map{clm}==dataRun_cell_ids(cells) ;
            masteri=clm ;
        end
    end
    if ~isempty(masteri) && ~isempty(srf{masteri}) ; % if there is a mapped master cell that has a rf
        subplot(8,l_cntrst,1) % spatial receptive field
        temp_rf = srf{masteri} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
        title(num2str(dataRunMaster_cell_ids(masteri)))
        
        subplot(8,l_cntrst,2) % temporal receptive field
        plot(fliplr(trf_time),trf{masteri}) ;
        title(cell_type_list{masteri})
    end
    for dset=1:3 ;
        for cntrst=1:l_cntrst ; % for each contrast

            subplot(8,l_cntrst,l_cntrst+cntrst) % psth variance (spatial frequency)
            plot(log(Spatial_frequency),psth_var{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            %errorbar(log(Spatial_frequency),psthEsem_var_mean{dset}{1}{cntrst}(cells,:),psthEsem_var_sem{dset}{1}{cntrst}(cells,:),psthEsem_var_sem{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth variance')
            title(num2str(dataRun_cell_ids(cells)))
            hold on
            
            subplot(8,l_cntrst,l_cntrst*2+cntrst) % psth norm variance (spatial frequency)
            plot(log(Spatial_frequency),psth_var_norm{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('norm psth variance')
            %set(gca,'yscale','log')
            hold on
            
            subplot(8,l_cntrst,l_cntrst*3+cntrst) % psth range (spatial frequency)
            plot(log(Spatial_frequency),psth_range{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth range')
            hold on

            subplot(8,l_cntrst,l_cntrst*4+cntrst) % psth mean (spatial frequency)
            plot(log(Spatial_frequency),psth_mean{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth mean')
            hold on

            subplot(8,l_cntrst,l_cntrst*5+cntrst) % psth peak/mean  (spatial frquency)
            plot(log(Spatial_frequency),psth_peakDivMean{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth peak/mean')
            hold on
            
            subplot(8,l_cntrst,l_cntrst*6+cntrst) % psth duty  (spatial frquency)
            plot(log(Spatial_frequency),psth_duty{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth duty')
            hold on
        end
    end
    
    nxtval = input('next (cell num, 0=back, return=forward)') ;
    if isempty(nxtval) ;
        cells=cells+1 ;
    elseif nxtval == 0 ;
        cells=cells-1 ;
    elseif nxtval>0 ;
        cells=nxtval ;
    end
end

% psth [control, +PSEM, wash] organized by unique cell types for a single spatial freq and contrast
sp = 7 ; % spatial frequency default
cntrst = 1 ; % contrast default
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ; 
            plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        title(num2str(cell_i{uc}(cells)))
    end
end

% mse (how does drug change of psth compare to wash change of psth?)
figure
subplot(2,2,1)
for cells=1:l_cells ; % for each cell of this type 
    plot(MseRelWash{cells}(:),MseRelDrug{cells}(:),'k*') ;
    hold on
    plot([min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],[min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],'g')
    plot([MseRelWashThreshold,MseRelWashThreshold],[min(MseRelDrug{cells}(:)),max(MseRelDrug{cells}(:))],'r')
    plot([min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],MseRelDrugDivWashThreshold*[min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],'r')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end
axis tight
xlabel('Rel Mse Wash')
ylabel('Rel Mse Drug')
title('all responses')

subplot(2,2,2)
plot(MseRelWashLfMean,MseRelDrugLfMean,'k*') ;
hold on
plot([min(MseRelWashLfMean),max(MseRelWashLfMean)],[min(MseRelWashLfMean),max(MseRelWashLfMean)],'g')
plot([MseRelWashThreshold,MseRelWashThreshold],[min(MseRelDrugLfMean),max(MseRelDrugLfMean)],'r')
plot([min(MseRelWashLfMean),max(MseRelWashLfMean)],MseRelDrugDivWashThreshold*[min(MseRelWashLfMean),max(MseRelWashLfMean)],'r')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis tight
xlabel('Rel Mse Wash')
ylabel('Rel Mse Drug')
title('cell mean lf')
 
subplot(2,2,3)
for cells=1:l_cells ; % for each cell of this type 
    plotyy(MseRelHistX,MseRelDrug_DivMseRelWash_Hist,MseRelHistX,MseRelDrug_DivMseRelWash_Hist_cumsum) ;
    text(.5,.9,['%aboveThresh=',num2str(100*MseRelDrug_DivMseRelWash_FractAboveThresh)],'units','norm')
end
xlabel('Rel Mse drug/wash')
ylabel('num obs')
title('all responses')

subplot(2,2,4)
for cells=1:l_cells ; % for each cell of this type 
    plotyy(MseRelHistX,MseRelDrug_DivMseRelWashLfMean_Hist,MseRelHistX,MseRelDrug_DivMseRelWashLfMean_Hist_cumsum) ;
    text(.5,.9,['%aboveThresh=',num2str(100*MseRel_FractAboveDrugDivWashThresh)],'units','norm')
end
xlabel('Rel Mse drug/wash')
ylabel('num obs')
title('cell mean lf')

% mse (does magintude of relative drug change depend on array location or cell type?)
figure

subplot(2,1,1)
for c = 1:length(dataRunMaster_cell_ids) ; % for each with a BW mapped ref
    if ~isempty(cell_list_map{c}) ; % if there is a mapped cell
        lw = MseRelDrug_DivMseRelWashLfMean(slave_c(c))/max(MseRelDrug_DivMseRelWashLfMean) ; 
        
        if ~isnan(lw) ;
            [X,Y] = drawEllipse([ctr{c} rad{c} angle{c}]) ;
            if ~any(isnan([X,Y])) ;
                [X,Y] = tformfwd(coord_tform, X, Y) ;
                plot(X,Y,'k','color',[1-lw,1-lw,1-lw],'linewidth',lw*3)
                hold on
            end
        end
    end
end   

subplot(2,1,2)
for uc = 1:length(UniqueCellTypes) ;
    tempMean = 0 ;
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        plot(uc,MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells)),'ko')
        if ~isnan(MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))) 
            tempMean = tempMean+MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))/lc ;
        end
        hold on
    end
    plot(uc,tempMean,'r+','MarkerSize',30)
end

% LfMean parameter comparisons for cells with (black) and without (cyan) relatively stable cells
figure
for uc = 1:length(UniqueCellTypes) ;
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        if MseRelWashLfMean(cell_i{uc}(cells))<MseRelWashThreshold ; % if this cell had a drug effect
            PlotColor = 'k' ;
        else
            PlotColor = 'c' ;
        end
        
        if strcmp('ON',UniqueCellTypes{uc}(1:2))
            SignPnt = 'o' ;
        elseif strcmp('OF',UniqueCellTypes{uc}(1:2))
            SignPnt = '*' ;
        else
            SignPnt = '+' ;
        end
        
        subplot(6,2,1)
        plot(uc,psth_var_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('Variance drug/cntrl')

        subplot(6,2,3)
        plot(uc,psth_range_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('Range drug/cntrl')

        subplot(6,2,5)
        plot(uc,psth_mean_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('Mean drug/cntrl')

        subplot(6,2,7)
        plot(uc,psth_duty_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('Duty drug/cntrl')

        subplot(6,2,9)
        plot(uc,psth_peakDivMean_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('peak/mean drug/cntrl')
    end
end
subplot(6,2,1); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,3); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,5); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,7); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,9); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,11); plot([1,length(UniqueCellTypes)],[0,0])

subplot(6,2,2)
plot(psth_var_divCntrl_LfMean_histX,psth_var_DrugdivCntrl_LfMean_hist)
xlabel('var drug/cntrl')
ylabel('# obs')

subplot(6,2,4)
plot(psth_range_divCntrl_LfMean_histX,psth_range_DrugdivCntrl_LfMean_hist)
xlabel('range drug/cntrl')
ylabel('# obs')

subplot(6,2,6)
plot(psth_mean_divCntrl_LfMean_histX,psth_mean_DrugdivCntrl_LfMean_hist)
xlabel('mean drug/cntrl')
ylabel('# obs')

subplot(6,2,8)
plot(psth_duty_divCntrl_LfMean_histX,psth_duty_DrugdivCntrl_LfMean_hist)
xlabel('duty drug/cntrl')
ylabel('# obs')

subplot(6,2,10)
plot(psth_peakDivMean_divCntrl_LfMean_histX,psth_peakDivMean_DrugdivCntrl_LfMean_hist)
xlabel('peak/mean drug/cntrl')
ylabel('# obs')

% contrast gain functions for each cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ; 
            plot(Spatial_frequency,psth_var_cgain{dset}{1}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        set(gca,'Xscale','log')
        title(num2str(cell_i{uc}(cells)))
    end
end
    
% spatial tuning functions averaged across contrasts for each cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ; 
            plot(Spatial_frequency,psth_var_cMean{dset}{1}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        set(gca,'Xscale','log')
        title(num2str(cell_i{uc}(cells)))
    end
end

% parameter changes as function of sp averaged across contrast
% organized by cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        PlotColor = 'c' ;
        if MseRelWashLfMean(cell_i{uc}(cells))<MseRelWashThreshold ; % if this cell was stable
            PlotColor = 'b' ;
        end
        if MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))>MseRelDrugDivWashThreshold ; % if this cell had a drug change > wash change
            PlotColor = 'r' ;
        end
        if MseRelWashLfMean(cell_i{uc}(cells))<MseRelWashThreshold & MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))>MseRelDrugDivWashThreshold ; % if this cell was stable and had a bigger drug effect
            PlotColor = 'k' ;
        end
        
        subplot(5,1,1) 
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on
        
        subplot(5,1,2) 
        plot(Spatial_frequency,psth_var_cgain_DivCntrl{2}(cell_i{uc}(cells),:), PlotColor)
        hold on
        
        subplot(5,1,3) 
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on

        subplot(5,1,4) 
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on
    end
    
    subplot(5,1,1) 
    %plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('var drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,2) 
    %plot(Spatial_frequency,psth_var_cgain_DrugDivCntrl_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('var gain drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,3) 
    %plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('mean drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,4) 
    %plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('duty drug/cntrl')
    xlabel('Spatial Frequency (c/um)')
end

% parameter changes as function of sp averaged across contrast
% comparing averages across cell types
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(6,1,1) 
    plot(Spatial_frequency,nanmean(psth_var_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,2) 
    plot(Spatial_frequency,nanmean(psth_range_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,3) 
    plot(Spatial_frequency,nanmean(psth_mean_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,4) 
    plot(Spatial_frequency,nanmean(psth_duty_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,5) 
    plot(Spatial_frequency,nanmean(psth_peakDivMean_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight
end

% average stats for all cell types for 'stable cells' 
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs 
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr{Master_i{uc}(c)} rad{Master_i{uc}(c)} angle{Master_i{uc}(c)}]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i{uc}(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'r')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end
        
% cells with large Drug mse /wash mse 
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr{Master_i{uc}(c)} rad{Master_i{uc}(c)} angle{Master_i{uc}(c)}]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i{uc}(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'r')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end

% cells that are stable and have large drug/wash mse
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr(Master_i{uc}(c)) rad(Master_i{uc}(c)) angle(Master_i{uc}(c))]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'r')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end
    
        
    
    


