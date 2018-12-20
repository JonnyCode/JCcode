RootAnalysisDir = '/Volumes/lab/Experiments/Array/Analysis/' ; % Brahms server
RootMovieDir = '/Volumes/lab/acquisition/movie-xml/' ; % Brahms server xml movie root 

Hbin = 0.1 ; % (sec)

% load BW data
dataRunMaster = load_data([RootAnalysisDir,'2016-06-22-0/data000/data000']) ;
dataRunMaster = load_neurons(dataRunMaster) ;
dataRunMaster = load_ei(dataRunMaster, 'all') ;
dataRunMaster = load_params(dataRunMaster,'cell_type_depth', 5) ;
dataRunMaster = load_sta(dataRunMaster) ; % only necessary to get trf_time

% Spatial receptive fields
for c = 1:length(dataRunMaster.cell_ids) ; % for each with a BW mapped ref
    ctr(c,:) = dataRunMaster.stas.fits{c}.mean ;
    rad(c,:) = dataRunMaster.stas.fits{c}.sd ;
    angle(c) = dataRunMaster.stas.fits{c}.angle ;
end

ctr_dist = squareform(pdist(ctr)) ; % distance between srf centers

for c = 1:length(dataRunMaster.cell_ids) ; % for each cell (point cell)
    
    % cells that likely have a lot of receptive field overlap
    CellSet = find(ctr_dist(c,:)< min(rad(c,:))) ; 
    
    % make hisograms
    for csi = 1:length(CellSet) ;
        H(csi,:) = hist(dataRunMaster.spikes{CellSet(csi)},[0:Hbin:dataRunMaster.duration]) ;
    end

    % identify clear responses (1) and non-responses (2), anything in the noise is (nan)
    B = nans(size(H)) ;
    for csi = 1:length(CellSet) ;
        UpperThresh = mean(H(csi,:))+2*std(H(csi,:)) ;
        B(csi,H(csi,:)>=UpperThresh)= 1 ;
        
        LowerThresh = max((mean(H(csi,:))-2*std(H(csi,:))),min(H(csi,:))) ;
        B(csi,H(csi,:)<=LowerThresh)=0 ;
    end
    
    % find clear stim responses in point cell 
    r = find(B(1,:)==1) ;
    
    % make respones pairs (1=+/-, 0=-/- or +/+, -1=-/+ pairs)
    for csi = 1:length(CellSet) ;
        for a=1:length(r) ;
            R{csi} = B(csi,r(a)) - B(csi,:) ; % difference in stim response logicals  
        end
    end
    
    % find clear response/non-response (+/-) stim response pairs in point cell (R=1)
    L{1} = R{1}==1 ;
    
    % find -/+, -/-, and +/+ stim response pairs in the non-point cells (i.e. stim
    % responses that are cell different than the point cell)
    for csi = 2:length(CellSet) ;
        L{csi} = R{csi}==0|-1 ;
    end
    
    % find unique stim pair times in point cells
    P=L{1} ;
    for csi = 2:length(CellSet) ;
        P = P*L{csi} ;
    end
    
    % pull movies leading to that those stim pairs (UNFINISHED)
    for a=1:size(P,2) ; % for each positive stim in point cell
        for b=1:sum(P(a,:)) ; % for each unique pair
            PosStim{a} = movie ;
            NegStim{a}{b} = movie ;
        end
    end     
end
        
figure
for c = 1:length(dataRunMaster.cell_ids) ; % for each with a BW mapped ref
    clf
    
    % cells that likely have a lot of receptive field overlap
    CellSet = find(ctr_dist(c,:)< min(rad(c,:))) ;
    
    for csi = 1:length(CellSet) ;
        [X,Y] = drawEllipse([ctr(CellSet(csi),:) rad(CellSet(csi),:) angle(CellSet(csi))]) ;
        if ~any(isnan([X,Y])) ;
            subplot(2,1,1)
            if CellSet(csi)==c ;
                plot(X,Y,'r')
                hold on
                plot(ctr(CellSet(csi),1),ctr(CellSet(csi),2),'r*')
            else
                plot(X,Y,'k')
                plot(ctr(CellSet(csi),1),ctr(CellSet(csi),2),'g*')
            end
        end
        title(num2str(c))
        
        subplot(2,1,2)
        plot((H(csi,:)-mean(H(csi,:)))/std(H(csi,:)))
        hold on
    end
    pause
end


         