% script to calculate some extra vectors from MbAnalyzer data
% JC 3/25/1016

for DB=1:5 ;
    for MbType = 1:4 ; % for each moving bar group
        % load mat file of all data vectors
        load([DataBlock(DB).DSmatFilePath,'Db',num2str(DB),'Mbtype',num2str(MbType),'MbAnalysisMat.mat']) ;

        first_phase_width_cHist{DB}(MbType,:) = cumsum(first_phase_width_Hist)/sum(first_phase_width_Hist) ; % phase width
        second_phase_width_cHist{DB}(MbType,:) = cumsum(second_phase_width_Hist)/sum(second_phase_width_Hist) ; % phase width
        first_phase_dsi_cHist{DB}(MbType,:) = cumsum(first_phase_dsi_Hist)/sum(first_phase_dsi_Hist) ; % phase vector dsi
        second_phase_dsi_cHist{DB}(MbType,:) = cumsum(second_phase_dsi_Hist)/sum(second_phase_dsi_Hist) ; %

        PhaseDirDiff_cHist{DB}(MbType,:) = cumsum(PhaseDirDiff_Hist)/sum(PhaseDirDiff_Hist) ; 
        
        N(DB,MbType) = sum(first_phase_width_Hist) ;
        
        clearvars -except DataBlock DB N first_phase_width_cHist second_phase_width_cHist first_phase_dsi_cHist second_phase_dsi_cHist PhaseDirDiff_cHist phase_width_HistX phase_dsi_HistX phase_angle_HistX
    end
end


for MbType = 1:4 ; % for each moving bar group
    first_phase_width_cHist_Ko(MbType,:) = mean([first_phase_width_cHist{1}(MbType,:);first_phase_width_cHist{2}(MbType,:)]) ;
    second_phase_width_cHist_Ko(MbType,:) = mean([second_phase_width_cHist{1}(MbType,:);second_phase_width_cHist{2}(MbType,:)]) ;

    first_phase_width_cHist_Het(MbType,:) = mean([first_phase_width_cHist{3}(MbType,:);first_phase_width_cHist{4}(MbType,:)]) ;
    second_phase_width_cHist_Het(MbType,:) = mean([second_phase_width_cHist{3}(MbType,:);second_phase_width_cHist{4}(MbType,:)]) ;
    
    first_phase_dsi_cHist_Ko(MbType,:) = mean([first_phase_dsi_cHist{1}(MbType,:);first_phase_dsi_cHist{2}(MbType,:)]) ;
    second_phase_dsi_cHist_Ko(MbType,:) = mean([second_phase_dsi_cHist{1}(MbType,:);second_phase_dsi_cHist{2}(MbType,:)]) ;

    first_phase_dsi_cHist_Het(MbType,:) = mean([first_phase_dsi_cHist{3}(MbType,:);first_phase_dsi_cHist{4}(MbType,:)]) ;
    second_phase_dsi_cHist_Het(MbType,:) = mean([second_phase_dsi_cHist{3}(MbType,:);second_phase_dsi_cHist{4}(MbType,:)]) ;
    
    PhaseDirDiff_cHist_Ko(MbType,:) = mean([PhaseDirDiff_cHist{1}(MbType,:);PhaseDirDiff_cHist{2}(MbType,:)]) ; 
    PhaseDirDiff_cHist_Het(MbType,:) = mean([PhaseDirDiff_cHist{3}(MbType,:);PhaseDirDiff_cHist{4}(MbType,:)]) ;
end

figure
subplot(2,2,1)
plot(phase_dsi_HistX,first_phase_dsi_cHist_Ko(2,:),'r') 
hold on
plot(phase_dsi_HistX,first_phase_dsi_cHist_Het(2,:),'k') 
xlabel('DS index')
ylabel('Fraction of cells')
title ('On low contrast')

subplot(2,2,2)
plot(phase_dsi_HistX,second_phase_dsi_cHist_Ko(2,:),'r') 
hold on
plot(phase_dsi_HistX,second_phase_dsi_cHist_Het(2,:),'k') 
xlabel('DS index')
ylabel('Fraction of cells')
title ('Off low contrast')

subplot(2,2,3)
plot(phase_angle_HistX,PhaseDirDiff_cHist_Ko(2,:),'r') 
hold on
plot(phase_angle_HistX,PhaseDirDiff_cHist_Het(2,:),'k') 
xlabel('angle diff. (degrees)')
ylabel('Fraction of cells')
title ('low contrast')



figure
subplot(2,2,1)
plot(phase_width_HistX,first_phase_width_cHist_Ko(1,:),'r') 
hold on
plot(phase_width_HistX,first_phase_width_cHist_Het(1,:),'k') 
xlabel('Circular std (degrees)')
ylabel('Fraction of cells')
title ('On high contrast')

subplot(2,2,2)
plot(phase_width_HistX,first_phase_width_cHist_Ko(3,:),'r') 
hold on
plot(phase_width_HistX,first_phase_width_cHist_Het(3,:),'k') 
xlabel('Circular std (degrees)')
ylabel('Fraction of cells')
title ('Off high contrast')

subplot(2,2,3)
plot(phase_width_HistX,first_phase_width_cHist_Ko(2,:),'r') 
hold on
plot(phase_width_HistX,first_phase_width_cHist_Het(2,:),'k') 
xlabel('Circular std (degrees)')
ylabel('Fraction of cells')
title ('On low contrast')

subplot(2,2,4)
plot(phase_width_HistX,first_phase_width_cHist_Ko(4,:),'r') 
hold on
plot(phase_width_HistX,first_phase_width_cHist_Het(4,:),'k') 
xlabel('Circular std (degrees)')
ylabel('Fraction of cells')
title ('Off low contrast')

figure
subplot(2,2,1)
plot(phase_dsi_HistX,first_phase_dsi_cHist_Ko(1,:),'r') 
hold on
plot(phase_dsi_HistX,first_phase_dsi_cHist_Het(1,:),'k') 
xlabel('DS index')
ylabel('Fraction of cells')
title ('On high contrast')

subplot(2,2,2)
plot(phase_dsi_HistX,first_phase_dsi_cHist_Ko(3,:),'r') 
hold on
plot(phase_dsi_HistX,first_phase_dsi_cHist_Het(3,:),'k') 
xlabel('DS index')
ylabel('Fraction of cells')
title ('Off high contrast')


subplot(2,2,3)
plot(phase_dsi_HistX,first_phase_dsi_cHist_Ko(2,:),'r') 
hold on
plot(phase_dsi_HistX,first_phase_dsi_cHist_Het(2,:),'k') 
xlabel('DS index')
ylabel('Fraction of cells')
title ('On low contrast')


subplot(2,2,4)
plot(phase_dsi_HistX,first_phase_dsi_cHist_Ko(4,:),'r') 
hold on
plot(phase_dsi_HistX,first_phase_dsi_cHist_Het(4,:),'k') 
xlabel('DS index')
ylabel('Fraction of cells')
title ('Off low contrast')


