% parameters
signal_fraction = [.01:.01:.1,.2:.1:1] ;

signal_mean = 8 ;
signal_std = 4 ;
noise_mean = 0 ;
noise_std = 4 ;

numTrials = 100 ;
pnts = 10000 ;

hist_range = -100:100 ;

% create signal
signal = normrnd(signal_mean,signal_std,1,pnts) ;

for round = 1:length(signal_fraction) ;
    for trial = 1:numTrials ;

        % create noise
        noise = normrnd(noise_mean,noise_std,1,pnts) ;

        % histograms
        hist_signal = hist(signal,hist_range) ;
        hist_noise = hist(noise,hist_range) ;

        % nonlinearity
        nonlinearity = zeros(1,length(hist_range)) ;
        hist_sum = hist_signal+hist_noise ;
        nonlinearity(hist_sum~=0)  = hist_signal(hist_sum~=0) ./hist_sum(hist_sum~=0)  ; %Field and Rieke 2002
        nonlinearity = nonlinearity/max(nonlinearity) ;

        if trial<=signal_fraction(round)*numTrials ;
            input(trial,:) = noise + signal ;
        else
            input(trial,:) = noise ;
        end

        % implement nonlinear weighting
        for a=1:numel(input(trial,:)) ;
            output(trial,a) = input(trial,a)*nonlinearity(find(hist_range>=input(trial,a),1)) ;
        end

    end

    
    hist_output = hist(output(:),hist_range) ;
    hist_input = hist(input(:),hist_range) ;

    % average trials 
    output_mean = mean(output) ;
    input_mean = mean(input) ;

    % correlation of signal and input (no nonlinear filtering)
    r = corrcoef(input_mean',signal') ;
    corrInputSignal(round) = r(2,1) ;

    % correlation of signal and output (with filtering)
    r = corrcoef(output_mean',signal') ;
    corrOutputSignal(round) = r(2,1) ;
    
    % mse of signal and input (no nonlinear filtering)
    mseInputSignal(round) = mean((input_mean-signal).^2) ;
    
    % mse of signal and output (with filtering)
    mseOutputSignal(round) = mean((output_mean-signal).^2) ;

    

%     figure
%     plot(hist_range,hist_input/max(hist_input),'c')
%     hold on
%     plot(hist_range,hist_signal/max(hist_signal),'b')
%     plot(hist_range,hist_noise/max(hist_noise),'k')
%     plot(hist_range,hist_output/max(hist_output),'g')
%     plot(hist_range,nonlinearity,'r')
%     legend('input','signal','noise','output','nonlinearity')

end

% theoretical corr of signal and noise without non-linearity
%theoreticalCorr = (((signal_fraction*numTrials).^2)*signal_std^2)./((((signal_fraction*numTrials).^2)*signal_std^2)+((1-signal_fraction)*numTrials*noise_std^2)) ;
theoreticalCorr = ((signal_fraction.^2)*signal_std^2)./(((signal_fraction.^2)*signal_std^2)+((1-signal_fraction)*noise_std^2)/numTrials) ;
theoreticalCorr = ((signal_fraction.^2)*var(signal))./(((signal_fraction.^2)*var(signal))+(((1-signal_fraction)*var(noise))/numTrials)) ;

theoreticalSNR = (signal_fraction*numTrials)./sqrt((1-signal_fraction)*numTrials) ;

theoreticalMSE = (1-signal_fraction).^2*signal_std^2+(noise_std^2/numTrials) ;

figure
plot(signal_fraction,corrInputSignal,'b-*')
hold on
plot(signal_fraction,corrOutputSignal,'r-*')
plot(signal_fraction,theoreticalCorr,'b--o')

figure
plot(signal_fraction,corrOutputSignal./corrInputSignal,'k-*')

figure
plot(signal_fraction,mseInputSignal,'b-*')
hold on
plot(signal_fraction,mseOutputSignal,'r-*')
plot(signal_fraction,theoreticalMSE,'b--o')







