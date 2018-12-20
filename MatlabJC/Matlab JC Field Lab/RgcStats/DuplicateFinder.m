function Duplicate_cell_ids = DuplicateFinder(dataRunMaster,dataRunSlave,varargin)

% this function will take a master and slave datarun and look for
% duplicated cells using the srf and trf fits.  

% If looking for duplicatesin a single data set, set dataRunSlave= [].
% Duplicate_cell_ids{master_i} = {[master_cell_id],[slave_cell_ids]}

% J Cafaro 08/18/2016

verify_duplicates=1 ; % verify duplicates in figure

% threshold parameters
MeanDist_threshold = .5 ; % (units)
SdDist_threshold = .2 ; % (fraction) minimum difference between cells major or minor axis, also minum size of cell to check angle
AglDist_threshold = 90 ; % (degrees) minimum anlge difference between eliptical cells
SdDiff_threshold = .8 ; % (fraction min/maj axes) min value to consider cell not round
TrfDist_threshold = .75 ; % (dot product from norm vectors) min similar trf

% set up input
if isempty(dataRunSlave) ; % if there is no slave
    dataRunSlave = dataRunMaster ; % look for duplicates in itself
    SelfSearchFlag = 1 ; 
end

% number of cells
numCellsMaster =length(dataRunMaster.stas.fits) ; % number cells in master
numCellsSlave =length(dataRunSlave.stas.fits) ; % number cells in slave

% get important parameters from dataRuns
for cl=1:numCellsMaster ; % for each cell in master
    Master.mean(cl,:) = dataRunMaster.stas.fits{cl}.mean ;
    Master.sd(cl,:) = dataRunMaster.stas.fits{cl}.sd ;
    Master.angle(cl) = dataRunMaster.stas.fits{cl}.angle ;
    Master.tc(cl,:) = dataRunMaster.stas.time_courses{cl} ;
end

for cl=1:numCellsSlave ; % for each cell in slave
    Slave.mean(cl,:) = dataRunSlave.stas.fits{cl}.mean ;
    Slave.sd(cl,:) = dataRunSlave.stas.fits{cl}.sd ;
    Slave.angle(cl) = dataRunSlave.stas.fits{cl}.angle ;
    Slave.tc(cl,:) = dataRunSlave.stas.time_courses{cl} ;
end

% prep parameter matrices
TempMat = nans(numCellsMaster, numCellsSlave) ; 
Mdist = TempMat ;
MjrSdDist = TempMat ;
MnrSdDist = TempMat ;
TrfDist = TempMat ;
Master_SdDiff = TempMat ;
Slave_SdDiff = TempMat ;
AglDist = TempMat ; 

% find putative duplicate cells
MasterSlaveGroups = cell(numCellsMaster) ;
for Mcl=1:numCellsMaster ; % for each cell in master
    for Scl=1:numCellsSlave ; % for each cell in slave
        if ~(Mcl==Scl && SelfSearchFlag == 1)  ; % if not searching for duplicates in same data set
        
            Mdist(Mcl,Scl)= pdist([Master.mean(Mcl,:);Slave.mean(Scl,:)]) ; % distance of centers
            if Mdist(Mcl,Scl)<MeanDist_threshold ; % if the centers are close enough...

                MjrSdDist(Mcl,Scl) = abs(max(Master.sd(Mcl,:))- max(Slave.sd(Scl,:)))/max(max(Master.sd(Mcl,:)), max(Slave.sd(Scl,:))) ; % fractional difference in major axis
                if MjrSdDist(Mcl,Scl)<SdDist_threshold ; % if the cells are close enough in size...

                    TrfDist(Mcl,Scl) = (Master.tc(Mcl,:)/norm(Master.tc(Mcl,:)))*(Slave.tc(Scl,:)'/norm(Slave.tc(Scl,:))) ; % dot product of norms
                    if TrfDist(Mcl,Scl)>TrfDist_threshold ; % if the time courses are similar enough... 

                        if max(Master.sd(Mcl,:))>SdDist_threshold ; % if the cells are big enough...  

                            MnrSdDist(Mcl,Scl) = abs(min(Master.sd(Mcl,:))- min(Slave.sd(Scl,:)))/max(min(Master.sd(Mcl,:)), min(Slave.sd(Scl,:))) ; % fractional difference in minor axis                        
                            if MnrSdDist(Mcl,Scl)<SdDist_threshold ; % if the cells are close enough in size...                  

                                Master_SdDiff(Mcl,Scl) = min(Master.sd(Mcl,:))/max(Master.sd(Mcl,:)) ; % fraction of minor/major (extent circularity)
                                Slave_SdDiff(Mcl,Scl) = min(Slave.sd(Scl,:))/max(Slave.sd(Scl,:)) ;

                                if Master_SdDiff(Mcl,Scl)<SdDiff_threshold && Slave_SdDiff(Mcl,Scl)<SdDiff_threshold ; % if the cells are not almost round
                                    AglDist(Mcl,Scl) = acuteAngle(Master.angle(Mcl)*180/pi,Slave.angle(Scl)*180/pi) ; % difference in angle
                                    if AglDist(Mcl,Scl)<AglDist_threshold ; % if their angles are similar...
                                        MasterSlaveGroups{Mcl} = [MasterSlaveGroups{Mcl},Scl] ; % putative duplicate
                                    end
                                else % if the cells are basically round don't bother checking elipse angle
                                    MasterSlaveGroups{Mcl} = [MasterSlaveGroups{Mcl},Scl] ; % putative duplicate
                                end
                            end                          
                        else % if the cells are too small don't bother checking minor axis and angle
                            MasterSlaveGroups{Mcl} = [MasterSlaveGroups{Mcl},Scl] ; % putative duplicate
                        end
                    end
                end
            end
        end
    end
end
              
% inspect putative duplicates
if verify_duplicates==1 ; % if you want to verify the duplicates, check plots
    putDupFig = figure ;

    for Mcl=1:numCellsMaster ; % for each cell in master
        if ~isempty(MasterSlaveGroups{Mcl}) ; % if there is a putative duplicate
            figure(putDupFig) ; clf
            subplot(1,3,1)

            [X,Y] = drawEllipse([Master.mean(Mcl,:), Master.sd(Mcl,:), Master.angle(Mcl)]) ;
            plot(X,Y,'k')
            hold on

            subplot(1,3,2)
            plot(Master.tc(Mcl,:),'k')
            hold on

            Scl_out=[] ;
            for cl=1:length((MasterSlaveGroups{Mcl})) ; % for each putative duplicate
                Scl = MasterSlaveGroups{Mcl}(cl) ;

                [X,Y] = drawEllipse([Slave.mean(Scl,:), Slave.sd(Scl,:), Slave.angle(Scl)]) ;

                subplot(1,3,1)
                plot(X,Y,'r')
                title(['Master: ',num2str(Mcl)])

                tx1=text(.8,.9,{['Mdist: ',num2str(Mdist(Mcl,Scl))],['MjrSdDist: ',num2str(MjrSdDist(Mcl,Scl))],...
                    ['MnrSdDist: ',num2str(MnrSdDist(Mcl,Scl))],['MstSdDiff: ',num2str(Master_SdDiff(Mcl,Scl))],...
                    ['SlvSdDiff: ',num2str(Slave_SdDiff(Mcl,Scl))],['AglDist: ',num2str(AglDist(Mcl,Scl))]},...
                    'Units','normalize') ;

                subplot(1,3,2)
                plot(Slave.tc(Scl,:),'r')
                title(['Slave: ',num2str(Scl)])

                tx2=text(.1,.9,['TrfDist: ',num2str(TrfDist(Mcl,Scl))], 'Units','normalize') ;

                subplot(1,3,3)
                plot_ei_(dataRunMaster.ei.eis{Mcl},dataRunMaster.ei.position,0,'neg_color','k','pos_color','k')
                hold on
                plot_ei_(dataRunMaster.ei.eis{Scl},dataRunMaster.ei.position,0,'neg_color','r','pos_color','r')

                DupString = input('duplicate (n,y):','s') ;
                if strcmp(DupString,'n'); % if its not a duplicate
                    Scl_out=[Scl_out,Scl] ;
                end

                subplot(1,3,1)
                plot(X,Y,'g')

                subplot(1,3,2)
                plot(Slave.tc(Scl,:),'g')

                subplot(1,3,3)
                cla

                delete(tx1) % delete text
                delete(tx2)          
            end
            MasterSlaveGroups{Mcl} = setdiff(MasterSlaveGroups{Mcl},Scl_out) ; % clear cells that user indicated are not duplicates
        end
    end
end
 
% make output {[master_cell_id],[slave_cell_id]}
for Mcl=1:numCellsMaster ; % for each cell in master
    Duplicate_cell_ids{Mcl} = {[dataRunMaster.cell_ids(Mcl)],...
        [dataRunSlave.cell_ids(MasterSlaveGroups{Mcl})]} ;
end
        
    


    


