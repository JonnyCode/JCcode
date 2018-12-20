% ORN-PN simulation

% ORN noise and signal params
ORNfr = [0:50] ; % firing rate of single ORN at several concentrations of odor
ORNfrStd = 10 ; % std of signal 

ORNfrSpont = 4 ; % average firing rate of spont activity
ORNfrSpontStd = 4 ; % std of firing rate of spont activty

ORNfrToAssess = 1 ; % ability to discriminate this ORN response from spont. noise

% ORN-PN nonlinearity params
Rmax = 30 ; % maximum firing rate of PN
hmax = 2 ; % half max of ORN-PN function
m = 100 ; % scale factor of supression for this ORN
ORNSfr = [300:350] ; % firing rate of all ORNS to given odor

% make ORN sig and noise dist
ORNsigDist = gmdistribution(ORNfrToAssess,ORNfrStd) ;
ORNsigDist = pdf(ORNsigDist,ORNfr') ;

ORNspontDist = gmdistribution(ORNfrSpont,ORNfrSpontStd) ;
ORNspontDist = pdf(ORNspontDist,ORNfr') ;

% ORN-PN nonlinearity
PNfr = Rmax*(ORNfr.^1.5./(ORNfr.^1.5+m*(ORNSfr./190)+hmax)) ; % PN firing rate

% figure
figure % ORN signal and noise dist
plot(ORNfr,ORNsigDist)
hold on
plot(ORNfr,ORNspontDist,'r')

figure % ORN-PN function
plot(ORNfr,PNfr) 