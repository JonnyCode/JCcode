% script will run a set of simulated filters to test the ability of LN
% model to encode movie frames in history independant way

% JC 2018/1/8

%% simulation

% trf params
Trf_tpeak = 0.05 ; % (sec)
Trf_peakRise = 0.01 ; % (sec) 
Trf_ttrough = 0.1 ; % (sec)
Trf_troughDecay = 0.06 ; % (sec)
Trf_peakAmp = 1 ; % (sec)
Trf_troughAmp = 0.2 ; % (sec)

% load image slide movie
load('/Volumes/lab/Documents/Movies/squirrel_video/squirrel_mean117_sd62_0to255.mat')

% reverse movie
movFlip = flip(mov,3) ;

% length params
numf = size(mov,3) ;
numX = size(mov,1) ;
numY = size(mov,2) ;

% time
Time = [1:numf]*1/60 ;

% construct temporal recptive field (trf)

Trf_peak = exp(-((Time-Trf_tpeak).^2)/(2*Trf_peakRise^2)) ;  % make gaussian for peak
Trf_trough = exp(-((Time-Trf_ttrough).^2)/(2*Trf_troughDecay^2)) ;  % make gausian for trough
Trf_exp = exp((-Trf_tpeak ./Time).^3) ;                 % make sat exponential to avoid tough before peak
Trf = (Trf_peakAmp*Trf_peak)-((Trf_troughAmp*Trf_trough).*Trf_exp) ; % add them up with exp weighting on trough 
Trf = Trf/sum(Trf(:)) ; % normalize

% non-DS filter
Params.centers = [numX,numY]/2 ; % (pix) X,Y positions of filter centers (Nx2), 
Params.direction_preferences = 90 ; % (deg) direction preference of each filter (Nx1)
Params.speed = 8 ; % (pix/frame) speed preference of filters
Params.SrfFlag = true ;
Params.radi = 150/4 ;
Params.Trf = Trf(1:20) ;

% make filter and generate linear prediction
Lp_mov = LinDirFilters(mov, Params) ;
Lp_movFlip = LinDirFilters(movFlip, Params) ;

cc = xcov(Lp_mov,fliplr(Lp_movFlip),'coef') ;

% square Lp
Lp2_mov = Lp_mov.^2 ;
Lp2_movFlip = Lp_movFlip.^2 ;

% rectify Lp
Thresh = 0 ; % fraction of max below which = 0 
Lp3_mov = Lp2_mov ;
Lp3_movFlip = Lp2_movFlip ;

ri = Lp2_mov<max(Lp2_mov)*Thresh ;
Lp3_mov(ri) = 0;

ri = Lp2_movFlip<max(Lp2_movFlip)*Thresh ;
Lp3_movFlip(ri) = 0;


