% Testing NLP model with STA and FBP

% JC 1/14/15

% parameters
Fbp_flag = true ;
TimeMax = 300 ; % (sec) time of experiment
TimeStep = .001 ; % (sec) time step

Spatial_N = 100 ; % number of spatial pnts
RefreshRate = .01 ; % (sec) time of each display frame

Wn_NperStixel = 5 ; % number of spatial pnts per whitenoise stixel
Wn_Ivalues = [0,128,256] ; % intensity values stixels can assume

Fbp_width = 10 ; % bar width
Fbp_angleValues = [0:30:150] ; % an
Fbp_offsetValues = [-45:45] ; % offset range
Fbp_Ivalue = 128 ; % intensity

Srf_type = 'SumOfGaussian' ; % function generating spatial receptive field
Trf_type = 'SumOfGaussian' ; % function generating temporal receptive field
Nl_type = 'linearThresh' ; % function generating temporal receptive field

Srf_center_center = [50,50] ; % receptive field center
Srf_center_cov = [50,0;0,50] ; % covariance matrix
Srf_surround_cov = [100,0;0,100] ; % covariance of surround
Srf_surround_strength = .09 ; % integral of surround relative to center

Trf_tpeak = 0.05 ; % (sec)
Trf_peakRise = 0.01 ; % (sec) 
Trf_ttrough = 0.1 ; % (sec)
Trf_troughDecay = 0.06 ; % (sec)
Trf_peakAmp = 1 ; % (sec)
Trf_troughAmp = 0.2 ; % (sec)

Nl_minSpikeRate = 0 ; % hz
Nl_maxSpikeRate = 500 ; % hz
Nl_minGenerator = 0 ; % generator
Nl_maxGenerator = 256 ; % generator

refractory_amp = 1000 ; % (Hz)
refractory_tau = 0.0001 ; % (sec)

SpikeRateNoise_std = 0 ; % (Hz)

StaTime = .1 ; % (sec)
StaSigStd = 3 ; % std of intensities that are significant
StaShift = 0 ; % (sec) time before spike sta starts

% from parameters
Time = [TimeStep:TimeStep:TimeMax] ;
TimePnts = length(Time) ; % number of time points
Generator = nan(1,TimePnts) ; % preallocate

%% model

% contruct spatial receptive field
if strcmp(Srf_type,'SumOfGaussian') ;
    [grid_x,grid_y] = meshgrid([1:Spatial_N],[1:Spatial_N]) ; 
    
    temp = mvnpdf([grid_x(:),grid_y(:)], Srf_center_center, Srf_center_cov) ;
    Srf_center = reshape(temp, Spatial_N, Spatial_N) ;
    
    temp = mvnpdf([grid_x(:),grid_y(:)], Srf_center_center, Srf_surround_cov) ;
    Srf_surround = -reshape(temp, Spatial_N, Spatial_N)*Srf_surround_strength ;
    
    Srf = Srf_center + Srf_surround ;
    
    Srf = Srf/sum(Srf(:)) ; % normalize
end
    
% construct temporal recptive field
if strcmp(Trf_type,'SumOfGaussian') ;
    Trf_peak = exp(-((Time-Trf_tpeak).^2)/(2*Trf_peakRise^2)) ;  % make gaussian for peak
    Trf_trough = exp(-((Time-Trf_ttrough).^2)/(2*Trf_troughDecay^2)) ;  % make gausian for trough
    Trf_exp = exp((-Trf_tpeak ./Time).^3) ;                 % make sat exponential to avoid tough before peak
    Trf = (Trf_peakAmp*Trf_peak)-((Trf_troughAmp*Trf_trough).*Trf_exp) ; % add them up with exp weighting on trough 
    
    Trf = Trf/sum(Trf(:)) ; % normalize
end

% construct nonlinearity
Nlx = [Nl_minGenerator-1,Nl_minGenerator,Nl_maxGenerator,Nl_maxGenerator+1] ;
Nly = [Nl_minSpikeRate,Nl_minSpikeRate,Nl_maxSpikeRate,Nl_maxSpikeRate] ;

%% stimuli
% FOR MEMORY AND SPEED THINK ABOUT: "matfile"- read and write only chuck of
% movie at a time, switch to all singles, try single frame random number
% generation, switch computers, ....

if ~Fbp_flag ; % if white noise
    % make white noise stimulus
    Wn_TperFrame = RefreshRate/TimeStep ;
    FrameStep = [0:Wn_TperFrame:TimePnts] ; % frame steps
    StixelStep = [0:Wn_NperStixel:Spatial_N] ; % stixel steps
    Wn = nan(Spatial_N,Spatial_N,TimePnts) ; % preallocate

    for t = 1:length(FrameStep)-1 ; % for each unique movie frame
        for a = 1:length(StixelStep)-1 ; % for each stixel X
            for b = 1:length(StixelStep)-1 ; % for each stixel Y
                Wn([StixelStep(a)+1:StixelStep(a+1)],[StixelStep(b)+1:StixelStep(b+1)],[FrameStep(t)+1:FrameStep(t+1)]) = Wn_Ivalues(randi(length(Wn_Ivalues))) ; % 
            end
        end
    end
else % if bars
    % construct stimulus for FBP
    Wn_TperFrame = RefreshRate/TimeStep ;
    FrameStep = [0:Wn_TperFrame:TimePnts] ; % frame steps
    Wn = nan(Spatial_N,Spatial_N,TimePnts) ; % preallocate

    for t = 1:length(FrameStep)-1 ; % for each unique movie frame
        if rem(t,2)~=0 ; % if odd   
            temp = nan(Spatial_N) ;
            angle(t) = Fbp_angleValues(randi(length(Fbp_angleValues))) ; % select random angle
            offset(t) = Fbp_offsetValues(randi(length(Fbp_offsetValues))) ; % select random offset
            for a = 1:Spatial_N ; % for each pixel X
                for b = 1:Spatial_N ; % for each pixel Y
                    temp(a,b) = (a-Spatial_N/2)*cosd(angle(t)) + (b-Spatial_N/2)*sind(angle(t)) - offset(t)  ;
                end
            end
            bar = (temp>=-Fbp_width/2 & temp<=Fbp_width/2)*Fbp_Ivalue ;
            for t2 = 1:Wn_TperFrame ;
                Wn(:,:,FrameStep(t)+t2) = bar ;
            end
        else
            angle(t) = angle(t-1) ; % 
            offset(t) = offset(t-1) ; % 
            for t2 = 1:Wn_TperFrame ;
                Wn(:,:,FrameStep(t)+t2) = zeros(Spatial_N) ;
            end
        end
    end
end

%% run stim on model

% run model with white noise stimulus
for t = 1:TimePnts ; % for each time point
    Generator(t) = sum(sum(Wn(:,:,t).*Srf)) ; % dot product with Srf
end
    
% convolve with Trf
Generator = conv(Generator,Trf) ;
Generator = Generator(1:TimePnts) ;

% nonlinearity
SpikeRate = interp1(Nlx, Nly,Generator, 'linear', 'extrap') ;

% additive noise in spike generation
SpikeRateNoise = normrnd(0, SpikeRateNoise_std, 1, TimePnts) ;
SpikeRate = SpikeRate + SpikeRateNoise ;

% poisson spike generation
Refractory_curve = -refractory_amp * exp(-Time/refractory_tau) ;
SpikeTrain = zeros(1,length(SpikeRate)) ;
SpikeRate_mod = SpikeRate ; % spike rate that is modified for refractory period

for a = 1:TimePnts ;
    if SpikeRate_mod(a)*TimeStep>random('unif',0,1) ; % 
        SpikeTrain(a) = 1 ;
        SpikeRate_mod(a:end) = SpikeRate_mod(a:end)+ Refractory_curve(1:TimePnts-a+1) ; %
    end
end
    
SpikePnts = find(SpikeTrain==1) ;

%% estimate model

if ~Fbp_flag ;
    % Estimate model using STA
    StaShiftPnts = StaShift/TimeStep ;
    for t=1:StaTime/TimeStep ;
        StaSpikePnts = SpikePnts-t+1+StaShiftPnts ;
        Sta(:,:,t) = mean(Wn(:,:,StaSpikePnts(StaSpikePnts>0 & StaSpikePnts<TimePnts)),3) ; 
    end

    % receptive field estimate
    Sta_high = mean(Sta(:)) + StaSigStd * std(Sta(:)) ; % highest Sta values
    Sta_low = mean(Sta(:)) - StaSigStd * std(Sta(:)) ; % lowest Sta values

    for t=1:StaTime/TimeStep ;
        Sta_sig_i = Sta(:,:,t)>Sta_high | Sta(:,:,t)<Sta_low ;
        Temp = Sta(:,:,t).*Sta_sig_i ;
        Trf_StaEst(t) = sum(Temp(:))/sum(Sta_sig_i(:)) ; % temporal receptive field estimate
    end

    Trf_weights = Trf_StaEst/nansum(Trf_StaEst) ;
    for t=1:StaTime/TimeStep ;
        Srf_Weighted(:,:,t) = Sta(:,:,t)*Trf_weights(t) ;
    end
    Srf_StaEst = nanmean(Srf_Weighted,3) ;
    
else
    % Estimate model using FBP
    RadonTransform_On = zeros(length(Fbp_offsetValues),length(Fbp_angleValues)) ;
    RadonTransform_Off = zeros(length(Fbp_offsetValues),length(Fbp_angleValues)) ;

    for t = 1:length(FrameStep)-1 ; % for each unique movie frame
        offset_i = find(Fbp_offsetValues==offset(t)) ;
        angle_i = find(Fbp_angleValues==angle(t)) ;

        if rem(t,2)~=0 ; % if odd  
            OnResponse(t) = sum(SpikeTrain(FrameStep(t)+1:FrameStep(t+1))) ;
            RadonTransform_On(offset_i,angle_i) = RadonTransform_On(offset_i,angle_i) + OnResponse(t)/sum([angle==angle(t)].*[offset==offset(t)]) ;
        else
            OffResponse(t) = sum(SpikeTrain(FrameStep(t)+1:FrameStep(t+1))) ;
            RadonTransform_Off(offset_i,angle_i) = RadonTransform_Off(offset_i,angle_i) + OffResponse(t)/sum([angle==angle(t)].*[offset==offset(t)]) ;
        end

    end

    for ang = 1:length(Fbp_angleValues) ; % for each angle
        Xvals = find(RadonTransform_On(:,ang)~=0) ; 
        RadonTransform_On_Interp(:,ang) = interp1(Xvals,RadonTransform_On(Xvals,ang),[1:length(RadonTransform_On)],'linear') ;

        Xvals = find(RadonTransform_Off(:,ang)~=0) ; 
        RadonTransform_Off_Interp(:,ang) = interp1(Xvals,RadonTransform_Off(Xvals,ang),[1:length(RadonTransform_Off)],'linear') ;
    end


    for a = 1:Spatial_N ; % for each pixel X
        for b = 1:Spatial_N ; % for each pixel Y
            for ang = 1:length(Fbp_angleValues) ; % for each angle
                z = (a-Spatial_N/2)*cosd(Fbp_angleValues(ang)) + (b-Spatial_N/2)*sind(Fbp_angleValues(ang)) ;
                offset_i = find(Fbp_offsetValues==round(z)) ;
                if ~isempty(offset_i) ;
                    %check = [check;[a,b,offset_i,ang]] ;
                    Srf_est_On_z(a,b,ang) = RadonTransform_On_Interp(offset_i,ang) ;
                    Srf_est_Off_z(a,b,ang) = RadonTransform_Off_Interp(offset_i,ang) ;
                else
                    Srf_est_On_z(a,b,ang) = 0 ;
                    Srf_est_Off_z(a,b,ang) = 0 ;
                end

            end
            Srf_est_On(a,b) = sum(Srf_est_On_z(a,b,:),3) ;
            Srf_est_Off(a,b) = sum(Srf_est_Off_z(a,b,:),3) ;
        end
    end
end



    

