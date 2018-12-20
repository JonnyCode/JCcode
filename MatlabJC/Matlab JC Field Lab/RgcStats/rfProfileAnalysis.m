function ForIgor = rfProfileAnalysis(DataBlock, DB, Params) 

% get some stats from srf and trf

% JC sept 2016

% load data
dataRun = load_data(DataBlock(DB).BwPath{1}) ;
dataRun = load_neurons(dataRun) ;

marks_params.thresh = 3 ;
dataRun = load_sta(dataRun) ;
dataRun = load_params(dataRun) ;
dataRun = get_sta_summaries(dataRun, 'all','marks_params', marks_params) ;

dataRun.stimulus.monitor_refresh = 60 ; % correct refresh rate

numCells = length(dataRun.cell_ids) ;

% srf params
for cl=1:numCells ; 
    [srfProfileX{cl}, srfProfile{cl}] = get_rf_profiles(dataRun,dataRun.cell_ids(cl),'normalize','none') ;
    srfProfile_smooth{cl} = smooth(srfProfile{cl},5) ;
    Xi_10 = find(srfProfileX{cl}>10,1,'first') ; % find first point after 10 pix
    Xi_20 = find(srfProfileX{cl}>20,1,'first') ; % find first point after 20 pix
    srf_surroundFact(cl) = min(srfProfile_smooth{cl}(1:Xi_10))/std(srfProfile_smooth{cl}(Xi_10:Xi_20)) ;
end

[srf_surroundFact_hist,srf_surroundFact_histX] = hist(srf_surroundFact,[-10:.1:0]) ; % histogram of surround factors
srf_surroundFact_chist = cumsum(srf_surroundFact_hist) ;

% trf params
TrfTime = fliplr(-[0:dataRun.stas.depth-1]*dataRun.stimulus.interval/dataRun.stimulus.monitor_refresh) ; 

for cl=1:numCells ; 
    Trf(cl,:) = dataRun.stas.time_courses{cl} ;
    Trf_norm(cl,:) = Trf(cl,:)/norm(Trf(cl,:)) ;
    [tcFit(cl,:), final_params(cl,:)] = fit_time_course(Trf(cl,:)', 'verbose', false) ;

    tcParamsTemp= tc_params_finder(fliplr(tcFit(cl,:)),fliplr(TrfTime)) ;
    
    Trf_ZeroTime(cl) = tcParamsTemp.zeroCrosst ;
    Trf_DoT(cl) = tcParamsTemp.DoT ;
end

[trf_zeroTime_hist, trf_zeroTime_histX] = hist(Trf_ZeroTime,[-1:.01:0]) ; % histogram 
trf_zeroTime_chist = cumsum(trf_zeroTime_hist) ;

[trf_DoT_hist, trf_DoT_histX] = hist(Trf_DoT,[0:.1:2]) ; % histogram 
trf_DoT_chist = cumsum(trf_DoT_hist) ;

% for Igor
Field = ['SurroundFactHistDb',num2str(DB)]
ForIgor.(Field) = [srf_surroundFact_histX;srf_surroundFact_chist] ;

Field = ['ZeroTimeHistDb',num2str(DB)]
ForIgor.(Field) = [trf_zeroTime_histX;trf_zeroTime_chist] ;

Field = ['DotHistDb',num2str(DB)]
ForIgor.(Field) = [trf_DoT_histX; trf_DoT_chist] ;

% % figures to assess qualitative success of "surroundFact"
% 
% % load receptive fields and filter
% filt_params.radius = 0.75;
% dataRunMaster = get_rfs_filtered(dataRun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');
% 
% % srf
% for cl=1:numCells ; 
%     srf{cl} = dataRunMaster.stas.filt_rfs{cl}*dataRunMaster.stas.polarities{cl} ; % give the correct pix polarity 
% end
%    
% for cl=1:numCells ; 
%     
%     subplot(1,3,1); 
%     temp_rf = srf{cl} ;
%     norm_rf = norm_image(temp_rf);
%     imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
%     colormap(brewermap([],'RdBu'))
%     caxis([0,1]) 
%     set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
%     
%     subplot(1,3,2) ; 
%     plot(srfProfileX{cl}, srfProfile{cl})
%     hold on
%     plot(srfProfileX{cl},srfProfile_smooth{cl},'r')
%     hold off
%     
%     subplot(1,3,3);
%     plot(srf_surroundFact_histX,srf_surroundFact_hist)
%     hold on
%     plot([srf_surroundFact(cl),srf_surroundFact(cl)],[0,max(srf_surroundFact_hist)],'r')
%     hold off
%     
%     pause
% end
%     
    
    