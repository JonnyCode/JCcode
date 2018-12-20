function ForIgor = FirstOrderLNlight2G(Input,Parameters,id,A) ;

% This function will find the first order linear and nonlinear filters that convert
% light into conductance

% get data
try
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
catch
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);
end

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [stmOrig(a,:), error] = ITCReadEpochStm(epochs(a), 0, fp);  % get the light stimulus
    [data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp);    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SR(a) = 1/(SI(a) * 1e-6); % Sampling rate in Hz
    [seed(a),error] = ITCGetSeed(epochs(a),0,0,fp) ;
end

UniqueSeed = unique(seed) ; % unique ff light stim
if length(UniqueSeed)<2 ;
    errormessage('only one light seed used')
end

[prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

% rectifier, zero stim, and get rid of pre/post points
stm = stmOrig ;
negstmPnts = find(stm<0) ;  % find indicies which would be getting a negative stim voltage
stm(negstmPnts) = 0 ;       % make these points zero
stm = stm - repmat(mean(stm(:,prePnts:end-postPnts),2),1,size(stm,2)) ; % subtract off mean of stimulus during time varying stimulus
stm = stm(:,prePnts:end-postPnts) ;

% change I to conductance and get rid of pre/post points
if strcmp(id,'Exc')==1 ;
    G = data/-62 ;
else
    G = data/62 ;
end
G = G - repmat(mean(G,2),1,size(G,2)) ; % subtract mean 
G = G(:,prePnts:end-postPnts) ;

time = [1:length(G)]/SR(1) ;

% first order LN model
% linear filter from all data but one trial
% prep stim and response for lin filter
withheldi = find(seed == UniqueSeed(1)) ; % indices with same seed as first seed
notwithheldi =  find(seed ~= UniqueSeed(1)) ; % indices not being withheld

cutfactor = rem(length(stm)-Parameters.STAPnts,Parameters.ChopPnts) ; % how many pnts remain when we attempt to evenly divide the signal 
reshfactor = floor((length(stm)-Parameters.STAPnts)/Parameters.ChopPnts)*length(notwithheldi) ; % number of rows will we have when we divide up the signal further 

signal = stm(notwithheldi,Parameters.STAPnts+1:end-cutfactor) ; %leave out trials with withheld stim leave off any possible remainder before reshaping rows 
response = G(notwithheldi,Parameters.STAPnts+1:end-cutfactor) ;

signal = reshape(signal',Parameters.ChopPnts,reshfactor)' ; % cut stimuli in rows into more rows
response = reshape(response',Parameters.ChopPnts,reshfactor)' ;

[LinearFilter] = LinFilterFinder(signal,response, SR(1), 60) ; 

for trial = 1:size(data,1); % linear prediction
    LinPred(trial,:) = conv(LinearFilter,stm(trial,:)) ;
end
LinPred = LinPred(:,1:length(G)) ;

[InputBins,Nonlinearity,Nonlinearity_sem] = NonLinFilterFinder(LinPred(notwithheldi,:),G(notwithheldi,:),1) ; % static nonlinearity from all but withheld trials (input,oouput,binsize (pA))

% find repeats and get std and mean of g to see how well lnpred fits
for a=1:length(UniqueSeed) ;
    samestmi = find(seed==UniqueSeed(a)) ; 
    
    LNPred(a,:) = interp1(InputBins,Nonlinearity,LinPred(samestmi(1),:),'linear','extrap') ; % predict ln
    
    if length(samestmi)>1 ; 
        G_std(a,:) = std(G([samestmi],:),[],1) ;
        G_mean(a,:) = mean(G([samestmi],:),1) ;
    else
        G_std(a,:) = nans(1,length(G)) ;
        G_mean(a,:) = G_mean(samestmi,:) ;
    end
 
end

time_LF = [0:length(LinearFilter)-1]/SR(1)*-1 ;

% figure

% figure % linear filter
% plot(time_LF,LinearFilter)
% 
% figure % nonlinear filter
% plot(LinPred(notwithheldi,:),G(notwithheldi,:),'k.')
% hold on
% errorbar(InputBins,Nonlinearity,Nonlinearity_sem,Nonlinearity_sem,'r')
% 
% figure % withheld tested
% plot(time,G(withheldi,:))
% hold on
% plot(time,LNPred(1,:),'k--')

figure % all tested
plot_timestep = .5 ; % sec 
for a=1:length(UniqueSeed) ;
    plot(time,G_mean(a,:))
    hold on
    plot(time,G_mean(a,:)+(2*G_std(a,:)),'b--')
    plot(time,G_mean(a,:)-(2*G_std(a,:)),'b--')
    plot(time,LNPred(a,:),'r')
    hold off
    pause
    for b=1:ceil(max(time)/plot_timestep) ;
        xlim([b*plot_timestep-plot_timestep,b*plot_timestep])
        pause
    end    
end


% For Igor

identifier = ['LinearFilter',id,num2str(A)] ;
ForIgor.(identifier) = LinearFilter ;

identifier = ['timeLF',id,num2str(A)] ;
ForIgor.(identifier) = time_LF ;

identifier = ['InputBins',id,num2str(A)] ;
ForIgor.(identifier) = InputBins ;

identifier = ['NonLinearity',id,num2str(A)] ;
ForIgor.(identifier) = Nonlinearity ;

identifier = ['GmeanWH',id,num2str(A)] ;
ForIgor.(identifier) = G_mean(1,:) ;

identifier = ['GmeanPlus2stdWH',id,num2str(A)] ;
ForIgor.(identifier) = G_mean(1,:)+2*G_std(1,:) ;

identifier = ['GmeanMinus2stdWH',id,num2str(A)] ;
ForIgor.(identifier) = G_mean(1,:)-2*G_std(1,:) ;

identifier = ['LNPredWH',id,num2str(A)] ;
ForIgor.(identifier) = LNPred(1,:) ;

identifier = ['time',id,num2str(A)] ;
ForIgor.(identifier) = time ;



