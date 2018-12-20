%This script loads up and parses the data
%then runs the OLE and OQE script to do the analysis
%it makes plots of 1) error vs time for estimator;
%2) weights vs DSI for OLE and the linear OQE weight;
%3) and of interaction weighting vs distance for OQE

% load the data set
addpath('/Users/joelz/Dropbox/DS_2017_CafaroField/Ds Pop Coding');

load ForJoelSlide4;
DB = 11 ; 
ImPathNum = 3 ; %3 image data analyzed
DsPathNum = 2 ; % drifting grating data analyzed

% ForJoel(DB,ImPathNum,DsPathNum).stimulus = struct(stimulus) ;
% ForJoel(DB,ImPathNum,DsPathNum).psth = psth ; % {cell}{stimulus}(trial,:)
% ForJoel(DB,ImPathNum,DsPathNum).psth_mean = psth_mean ; % {cell}(stimulus,:)
% ForJoel(DB,ImPathNum,DsPathNum).psth_timeBins = psthTime ; % 
% ForJoel(DB,ImPathNum,DsPathNum).EiCenters = EiCtr ; % (cells)
% ForJoel(DB,ImPathNum,DsPathNum).DsCellsIndices = DsCelli ; % {DsType}(cell)

%extract the data srcuture of responses
dat_struct_img = ForJoel(DB,ImPathNum,DsPathNum).psth; % {cell}{stimulus}(trial,:)
dat_struct_grating = ForJoel(DB,ImPathNum,DsPathNum).psthDg;
DS_cell_IDs = ForJoel(DB,ImPathNum,DsPathNum).DsCellsIndices;  %{DsType}(cell)
DSlist = [DS_cell_IDs{1} DS_cell_IDs{2} DS_cell_IDs{3} DS_cell_IDs{4}];
Ncells = 392;
NONDS = setdiff((1:Ncells),DSlist);
EIcenters = ForJoel(DB,ImPathNum,DsPathNum).EiCenters;

%pick an array of cells: celllist
%celllist = DSlist;   
celllist = find(EIcenters(:,1)<0 & EIcenters(:,2)<0); %here, I took one quadrant (all cells in that quadrant)
%DSlister = zeros(Ncells,1); DSlister(DSlist) = 1; celllist = find(EIcenters(:,1)<0 & EIcenters(:,2)<0 & DSlister ==1);



%shape the data and labels...
[ntrials, ntimepts] = size(dat_struct_img{1}{1});
for stim = 1:length(dat_struct_img{1})
        for cellnum = 1:length(dat_struct_img)
            dat_img(cellnum,stim,:,:) = dat_struct_img{cellnum}{stim};
        end
        dirn_img(stim,:,:) = (pi/4)*stim*ones(ntrials,ntimepts);
        timepts(stim,:,:) = repmat(1:ntimepts,ntrials,1); %use this later for error vs time...
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('data reshaped and ready')


%now take the responses that I want for the OLE and OQE
%and put them in an array that's Nobs x Ncells
%also make "angles" array that's Ntrials x 1

responses = reshape(dat_img,cellnum,8*ntrials*ntimepts)';
responses = responses(:,celllist); %just do this for the specified cells
angles = reshape(dirn_img,1,8*ntrials*ntimepts)';
timevals = reshape(timepts,1,8*ntrials*ntimepts)';

%now run the OLE and OQE
OLE_OQE_Angles_Regularized

%make a few plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first, OLE and OQE Errors vs Time
%in 20 time bins
time_binner = linspace(1,max(timevals),21);
times_test = timevals(test_set); %test_set is an array from the OLE_OQE script, identifying the time bins / trials used for testing

for tt = 1:20
    
    OLE_mean_err(tt) = mean(angle_deviations_OLE((times_test>time_binner(tt)) & (times_test<time_binner(tt+1))));
    OLE_std_err(tt) = std(angle_deviations_OLE((times_test>time_binner(tt)) & (times_test<time_binner(tt+1)))) / sqrt( sum( (times_test>time_binner(tt)) & (times_test<time_binner(tt+1))  ) );

    OQE_mean_err(tt) = mean(angle_deviations_OQE((times_test>time_binner(tt)) & (times_test<time_binner(tt+1))));
    OQE_std_err(tt) = std(angle_deviations_OQE((times_test>time_binner(tt)) & (times_test<time_binner(tt+1))))/sqrt( sum( (times_test>time_binner(tt)) & (times_test<time_binner(tt+1))  ) );
    
    timecenters(tt) = mean(times_test((times_test>time_binner(tt)) & (times_test<time_binner(tt+1))));
    
end
figure(1)
hold on
title('OLE (R) and OQE (B) Error vs Time; IMAGES')
errorbar(timecenters,OLE_mean_err*180/pi,OLE_std_err*180/pi,'r')
errorbar(timecenters,OQE_mean_err*180/pi,OQE_std_err*180/pi,'b')
xlabel('time after stim onset (34 ms bins)')
ylabel('mean estimator error (Deg.)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now some descriptive stats
%mean error for OLE and OQE (all times)
av_error_OLE = mean(angle_deviations_OLE)*180/pi
av_error_OQE = mean(angle_deviations_OQE)*180/pi
[h p] = ttest(angle_deviations_OLE,angle_deviations_OQE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
