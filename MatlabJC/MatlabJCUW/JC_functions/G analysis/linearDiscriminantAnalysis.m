function ForIgor = linearDiscriminantAnalysis(ForIgor,id)

% this function will look at the population data assembeled after
% 'DCdcPairAnalyzer3.m' create bubble plots for igor and assess the percent
% correct of each bubble plot for an ideal linear discriminant function.

% JC 10/19/11

for a=1:3 ; % for each bar

    % mean matrices
    identifier = ['sn1SimMean',id,'Bar',num2str(a),'AllMean'] ;
    MeanPc(a,1) = ForIgor.(identifier) ; 

    identifier = ['sn2SimMean',id,'Bar',num2str(a),'AllMean'] ;
    MeanPc(a,2) = ForIgor.(identifier) ;  

    identifier = ['sn1ShuffMean',id,'Bar',num2str(a),'AllMean'] ;
    MeanMc(a,1) = ForIgor.(identifier) ; 

    identifier = ['sn2ShuffMean',id,'Bar',num2str(a),'AllMean'] ;
    MeanMc(a,2) = ForIgor.(identifier) ; 

    % varaiance in covariance matrices
    identifier = ['sn1SimVar',id,'Bar',num2str(a),'AllMean'] ;
    CovPcPp{a}(1,1) = ForIgor.(identifier) ; 
    CovPcMp{a}(1,1) = ForIgor.(identifier) ;

    identifier = ['sn2SimVar',id,'Bar',num2str(a),'AllMean'] ;
    CovPcPp{a}(2,2) = ForIgor.(identifier) ; 
    CovPcMp{a}(2,2) = ForIgor.(identifier) ;

    identifier = ['sn1ShuffVar',id,'Bar',num2str(a),'AllMean'] ;
    CovMcPp{a}(1,1) = ForIgor.(identifier) ; 
    CovMcMp{a}(1,1) = ForIgor.(identifier) ;

    identifier = ['sn2ShuffVar',id,'Bar',num2str(a),'AllMean'] ;
    CovMcPp{a}(2,2) = ForIgor.(identifier) ; 
    CovMcMp{a}(2,2) = ForIgor.(identifier) ;

    % covariance in covariance matrices
    identifier = ['snCovsPcPp',id,'Bar',num2str(a),'AllMean'] ;
    CovPcPp{a}(1,2) = ForIgor.(identifier) ; 
    CovPcPp{a}(2,1) = ForIgor.(identifier) ; 

    identifier = ['snCovsMcPp',id,'Bar',num2str(a),'AllMean'] ;
    CovMcPp{a}(1,2) = ForIgor.(identifier) ; 
    CovMcPp{a}(2,1) = ForIgor.(identifier) ;

    identifier = ['snCovsPcMp',id,'Bar',num2str(a),'AllMean'] ;
    CovPcMp{a}(1,2) = ForIgor.(identifier) ; 
    CovPcMp{a}(2,1) = ForIgor.(identifier) ;

    identifier = ['snCovsMcMp',id,'Bar',num2str(a),'AllMean'] ;
    CovMcMp{a}(1,2) = ForIgor.(identifier) ; 
    CovMcMp{a}(2,1) = ForIgor.(identifier) ;

    % bubble plots
    figure
    
    h=error_ellipse(CovPcPp{a},MeanPc(a,:),'conf',.95) ; % confidence elipse
    sn_PcPp_ellipseX(a,:) = get(h,'XData') ;
    sn_PcPp_ellipseY(a,:) = get(h,'YData') ;

    h=error_ellipse(CovMcPp{a},MeanMc(a,:),'conf',.95) ; % confidence elipse
    sn_McPp_ellipseX(a,:) = get(h,'XData') ;
    sn_McPp_ellipseY(a,:) = get(h,'YData') ;
    
    h=error_ellipse(CovMcMp{a},MeanMc(a,:),'conf',.95) ; % confidence elipse
    sn_McMp_ellipseX(a,:) = get(h,'XData') ;
    sn_McMp_ellipseY(a,:) = get(h,'YData') ;
    
    close
end

% percent correct
snRange = [0:50] ;
pCorrectPcPp = linearDiscriminant(CovPcPp,MeanPc,snRange,snRange ) ;

pCorrectMcPp = linearDiscriminant(CovMcPp,MeanMc,snRange,snRange ) ;

pCorrectMcMp = linearDiscriminant(CovMcMp,MeanMc,snRange,snRange ) ;

% figure

figure % bubble plots

subplot(1,3,1)
plot(sn_PcPp_ellipseX',sn_PcPp_ellipseY','k')
axis([0 20 0 20])
text(.1,.9,num2str(pCorrectPcPp),'Units','normalized')
title('control')

subplot(1,3,2)
plot(sn_McPp_ellipseX',sn_McPp_ellipseY','g')
axis([0 20 0 20])
text(.1,.9,num2str(pCorrectMcPp),'Units','normalized')
title('minus EI')

subplot(1,3,3)
plot(sn_McMp_ellipseX',sn_McMp_ellipseY','y')
axis([0 20 0 20])
text(.1,.9,num2str(pCorrectMcMp),'Units','normalized')
title('minus all')

% ForIgors
for a=1:3 ; % for each bar 
    identifier = ['snPcPpEllipseX',id,'Bar',num2str(a)] ;
    ForIgor.(identifier) = sn_PcPp_ellipseX(a,:) ; 

    identifier = ['snPcPpEllipseY',id,'Bar',num2str(a)] ;
    ForIgor.(identifier) = sn_PcPp_ellipseY(a,:) ; 
    
    
    identifier = ['snMcPpEllipseX',id,'Bar',num2str(a)] ;
    ForIgor.(identifier) = sn_McPp_ellipseX(a,:) ; 
    
    identifier = ['snMcPpEllipseY',id,'Bar',num2str(a)] ;
    ForIgor.(identifier) = sn_McPp_ellipseY(a,:) ; 
    
    
    identifier = ['snMcMpEllipseX',id,'Bar',num2str(a)] ;
    ForIgor.(identifier) = sn_McMp_ellipseX(a,:) ; 
    
    identifier = ['snMcMpEllipseY',id,'Bar',num2str(a)] ;
    ForIgor.(identifier) = sn_McMp_ellipseY(a,:) ; 

end

identifier = ['ILDpCorrectPcPp',id] ;
ForIgor.(identifier) = pCorrectPcPp ; 

identifier = ['ILDpCorrectMcPp',id] ;
ForIgor.(identifier) = pCorrectMcPp ; 

identifier = ['ILDpCorrectMcMp',id] ;
ForIgor.(identifier) = pCorrectMcMp ; 



