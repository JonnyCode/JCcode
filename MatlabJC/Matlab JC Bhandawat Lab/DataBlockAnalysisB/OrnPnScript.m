BgConcentrationRange = [0,10.^[-7:-4]] ;
ConcentrationRange = 10.^[-7:-2] ;
ConcentrationRangeFit = 10.^[-7:.1:-2] ;

for a=2:length(BgConcentrationRange) ;
    ConcentrationRangeDivBg(a,:) = ConcentrationRange./BgConcentrationRange(a) ;
end

Conc(1,:) = [.7,.6,.5,.4,0,0,0,0] ;
Conc(2,:) = [1,1,1,1,.7,.6,.5,.4] ;
Conc(3,:) = [.7,.6,.5,.4,.3,.2,.1,0] ;
for a=1:size(ForIgor.PnDr,1) ; % each concentration is a matrix and each background is a row within that matrix
    colorMat{a} = [Conc(3,a),Conc(3,a),Conc(3,a); Conc(2,a),Conc(1,a),Conc(1,a); Conc(1,a),Conc(1,a),Conc(2,a);... 
        Conc(1,a),Conc(2,a),Conc(1,a);Conc(2,a),Conc(2,a),Conc(1,a);...
        Conc(2,a),Conc(1,a),Conc(2,a);Conc(1,a),Conc(2,a),Conc(2,a);...
        Conc(3,a),Conc(3,a),Conc(1,a)] ;
end

% glomerular transformation
for a=1:size(ForIgor.PnDr,1) ;
    figure(1)
    plot(log10(ConcentrationRange),ForIgor.PnDr(a,:))
    hold on 
    plot(log10(ConcentrationRange),ForIgor.OrnDr(a,:),'--')

    errorbarXY(ForIgor.OrnDr(a,:),ForIgor.PnDr(a,:),ForIgor.OrnDr_sem(a,:),ForIgor.PnDr_sem(a,:),'o',colorMat{5}(a,:),2) 

    FitCoefs{a} = nlinfit(ForIgor.OrnDr(a,:),ForIgor.PnDr(a,:),@SatFun,[max(ForIgor.PnDr(a,:)),50]) ;
    OrnDrInterp = [0:1:max(ForIgor.OrnDr(:))] ;
    SatFit(a,:) = SatFun(FitCoefs{a},OrnDrInterp) ;
    figure(2)
    plot(OrnDrInterp,SatFit(a,:),'color',colorMat{5}(a,:),'LineWidth',1.5)
end

% comparing adaptation in ORN and PN
for a = 2:size(ForIgor.OrnDivCon,1) ;
    errorbarXY(ForIgor.OrnDivCon(a,:),ForIgor.PnDivCon(a,:),ForIgor.OrnDivCon_sem(a,:),ForIgor.PnDivCon_sem(a,:),'*',colorMat{5}(a,:),3)
    hold on
    plot([0,2],[0,2],'k')
end

figure
for a = 2:length(BgConcentrationRange) ; % for each background (but 'wash')
    plot(log10(ConcentrationRange),ForIgor.OrnDivCon(a,:),'*','Color',colorMat{5}(a,:))
    hold on
    plot(log10(ConcentrationRange),ForIgor.PnDivCon(a,:),'o','Color',colorMat{5}(a,:))
end
set(gca,'xlim',[-9,-1])
legend('ORN','PN')

% adaptation models (fixed or variable ORN or glomerular transform)
% transforms
for a = 1:length(BgConcentrationRange) ; % for each background (but 'wash')
    PNmodel_adaptORN_adaptTrans(a,:) = interp1(OrnDrInterp,SatFit(a,:),ForIgor.OrnDr(a,:),'linear','extrap') ;
    PNmodel_adaptORN_fixTrans(a,:) = interp1(OrnDrInterp,SatFit(1,:),ForIgor.OrnDr(a,:),'linear','extrap') ;
    PNmodel_fixORN_adaptTrans(a,:) = interp1(OrnDrInterp,SatFit(a,:),ForIgor.OrnDr(1,:),'linear','extrap') ;
end
    
figure
for a = 1:length(BgConcentrationRange) ; % for each background 
    plot(log10(ConcentrationRange),ForIgor.PnDr(a,:),'-o','Color',colorMat{5}(a,:),'LineWidth',1.5) ; % real data
    hold on
    plot(log10(ConcentrationRange),PNmodel_adaptORN_adaptTrans(a,:),'--*','Color',colorMat{5}(a,:),'LineWidth',1.5)
end
legend('data','model')

figure
for a = 1:length(BgConcentrationRange) ; % for each background 
    plot(log10(ConcentrationRange),PNmodel_adaptORN_adaptTrans(a,:),'-*','Color',colorMat{5}(a,:))
    hold on
    plot(log10(ConcentrationRange),PNmodel_adaptORN_fixTrans(a,:),'--*','Color',colorMat{5}(a,:))
    plot(log10(ConcentrationRange),PNmodel_fixORN_adaptTrans(a,:),'-.*','Color',colorMat{5}(a,:))
end

figure
for a = 2:length(BgConcentrationRange) ; % for each background
    subplot(1,3,1)
    plot(log10(ConcentrationRangeDivBg(a,:)),PNmodel_adaptORN_adaptTrans(a,:),'-*','Color',colorMat{5}(a,:),'LineWidth',1.5)
    hold on
    
    subplot(1,3,2)
    plot(log10(ConcentrationRangeDivBg(a,:)),PNmodel_adaptORN_fixTrans(a,:),'-*','Color',colorMat{5}(a,:),'LineWidth',1.5)
    hold on
    
    subplot(1,3,3)
    plot(log10(ConcentrationRangeDivBg(a,:)),PNmodel_fixORN_adaptTrans(a,:),'-*','Color',colorMat{5}(a,:),'LineWidth',1.5)
    hold on
end
    
% adapation factors from model
for a = 1:length(BgConcentrationRange) ; % for each background 
    PNmodel_adaptORN_adaptTrans_interp(a,:) = interp1(log10(ConcentrationRange),PNmodel_adaptORN_adaptTrans(a,:),log10(ConcentrationRangeFit)) ;
    PNmodel_adaptORN_fixTrans_interp(a,:) = interp1(log10(ConcentrationRange),PNmodel_adaptORN_fixTrans(a,:),log10(ConcentrationRangeFit)) ;
    PNmodel_fixORN_adaptTrans_interp(a,:) = interp1(log10(ConcentrationRange),PNmodel_fixORN_adaptTrans(a,:),log10(ConcentrationRangeFit)) ;
end
    
for a = 2:length(BgConcentrationRange) ; % for each background
    Model_adaptORN_fixTrans_FracCont(a,:) = (PNmodel_adaptORN_adaptTrans_interp(1,:)-PNmodel_adaptORN_fixTrans_interp(a,:))./(PNmodel_adaptORN_adaptTrans_interp(1,:)-PNmodel_adaptORN_adaptTrans_interp(a,:)) ;
    Model_fixORN_adaptTrans_FracCont(a,:) = (PNmodel_adaptORN_adaptTrans_interp(1,:)-PNmodel_fixORN_adaptTrans_interp(a,:))./(PNmodel_adaptORN_adaptTrans_interp(1,:)-PNmodel_adaptORN_adaptTrans_interp(a,:)) ; 
end

figure
imagesc(Model_adaptORN_fixTrans_FracCont,[0,1])
figure
imagesc(Model_fixORN_adaptTrans_FracCont,[0,1])

% models of 

