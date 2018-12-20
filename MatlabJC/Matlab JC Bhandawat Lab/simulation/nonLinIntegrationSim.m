% simulation for linear vs. nonlinear summation

% params
numTrials = 1000 ;

numChannels = 600 ;

OdorRsp_mean = 10 ;
OdorRsp_std = 3 ;
SpontAct_mean = 5 ;
SpontAct_std = 3 ;
RspRange = [0:.1:100] ;

% number of channels
OdorRsp_numChannels = floor(sqrt(numChannels)*0.75) ;
SpontAct_numChannels = numChannels - OdorRsp_numChannels ;

% distributions
OdorRsp_Dist = gmdistribution(OdorRsp_mean,OdorRsp_std^2) ;
OdorRsp_Dist = pdf(OdorRsp_Dist,RspRange') ;

SpontAct_Dist = gmdistribution(SpontAct_mean,SpontAct_std^2) ;
SpontAct_Dist = pdf(SpontAct_Dist,RspRange') ;

% optimal nonlinearity
OptNonLin = OdorRsp_Dist./(OdorRsp_Dist+SpontAct_Dist) ;

% integration
for a = 1:numTrials ;
    % odor present
    OdorRsp = normrnd(OdorRsp_mean,OdorRsp_std,1,OdorRsp_numChannels) ;
    SpontAct = normrnd(SpontAct_mean,SpontAct_std,1,SpontAct_numChannels) ;
    Input = [OdorRsp,SpontAct] ;
    
    OdorRsp_LinOutput(a) = mean(Input) ;
    
    for b=1:length(Input) ;
        [m,mi] = min(abs(Input(b)-RspRange)) ;
        NonLinInput(b) = OptNonLin(mi)*Input(b) ;
    end
    OdorRsp_NonLinOutput(a) = mean(NonLinInput) ;
    
    % no odor present
    SpontAct = normrnd(SpontAct_mean,SpontAct_std,1,numChannels) ;
    Input = SpontAct ;
    
    SpontAct_LinOutput(a) = mean(Input) ;
    
    for b=1:length(Input) ;
        [m,mi] = min(abs(Input(b)-RspRange)) ;
        NonLinInput(b) = OptNonLin(mi)*Input(b) ;
    end
    SpontAct_NonLinOutput(a) = mean(NonLinInput) ;
    
end

% output distributions
[OdorRspLin_hist,OdorRspLin_histX] = hist(OdorRsp_LinOutput,RspRange) ;
[OdorRspNonLin_hist,OdorRspNonLin_histX] = hist(OdorRsp_NonLinOutput,RspRange) ;

[SpontActLin_hist,SpontActLin_histX] = hist(SpontAct_LinOutput,RspRange) ;
[SpontActNonLin_hist,SpontActNonLin_histX] = hist(SpontAct_NonLinOutput,RspRange) ;

% discriminability
sum(OdorRspLin_hist(OdorRspLin_hist>SpontActLin_hist))/sum([OdorRspLin_hist,SpontActLin_hist]) ;
  
sum(OdorRspNonLin_hist(OdorRspNonLin_hist>SpontActNonLin_hist))/sum([OdorRspNonLin_hist,SpontActNonLin_hist]) ;

% figure
figure
plot(RspRange,SpontAct_Dist)
hold on
plot(RspRange,OdorRsp_Dist,'r')
plot(RspRange,OptNonLin,'k')

figure
plot(OdorRspLin_histX,OdorRspLin_hist,'r')
hold on
plot(OdorRspNonLin_histX,OdorRspNonLin_hist,'r--')
plot(SpontActLin_histX,SpontActLin_hist,'b')
plot(SpontActNonLin_histX,SpontActNonLin_hist,'b--')


