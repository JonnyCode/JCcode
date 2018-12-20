% DSPopCoding figures with multiple datasets
% JC 2018-12-07

saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/2018Paper/ImJitterAnalysisV4_sta/Figs/'] ;

DB=19 ;
RunId1 = '19538' ;
RunId2 = '43564' ;
% figure 1
figure
colorList = {'k','r'} ;
for a=1:2 ;
    load([saveFigPath,'ImJitterAnalysisV4_sta_DB',num2str(DB),'_',eval(['RunId',num2str(a)])]) ; % save mat file
    for DsType=1:length(DsTypeName) ; % for each DS cell type
        for c = 1:length(cell_i{DsType}) ; % for each DS cell
            subplot(1,3,1)
            plot(StaTime(StaWindowPnts),StaDirJitterVector{DsType}(c,StaWindowPnts)*180/pi,colorList{a} )
            hold on

            subplot(1,3,2)
            plot(StaTime(StaWindowPnts),StaMagJitterVector{DsType}(c,StaWindowPnts)*180/pi,colorList{a})
            hold on
     
        end
    end
    
    subplot(1,3,3)
    plot(StaTime(StaWindowPnts),StaMagJitterVector_mean(StaWindowPnts)*180/pi,colorList{a})
    hold on
    plot(StaTime(StaWindowPnts),StaMagJitterVector_mean(StaWindowPnts)*180/pi...
    +StaMagJitterVector_sem(StaWindowPnts)*180/pi,[colorList{a},':'])
    plot(StaTime(StaWindowPnts),StaMagJitterVector_mean(StaWindowPnts)*180/pi...
    -StaMagJitterVector_sem(StaWindowPnts)*180/pi,[colorList{a},':'])
end

saveas(gcf,[saveFigPath,'StaSet_DB',num2str(DB),'Fig_',RunId1,'_',RunId2])
print(gcf, '-dpdf',[saveFigPath,'StaSet_DB',num2str(DB),'Fig_',RunId1,'_',RunId2])

