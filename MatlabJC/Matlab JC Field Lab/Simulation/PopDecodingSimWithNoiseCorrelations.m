% model testing how SNR, population size, and covariance impacts pop error 
% for correlated population when independence is assumed (or not)

% JC 2018-10-30

nt = 10000 ; % time points

% stim values
stimSet = [0,.5,1] ; % possible stim values 
stim = stimSet(randi([1,length(stimSet)],nt,1)) ; % stim - uniform dist of possible values

% cell params
g1 = 1 ; % stim gain cell 1
g2set = [0:.2:2] ; % stim gain cell 2 (all cell 2 will be the same)
C2NumSet = [1,2,10] ; % number of cell 2

Var = .5 ; % var for all cells and all stim (additive noise)
CorrSet = [0,.5,.99] ; % noise covariance

% decoder params
AssumeIndFlag = true ;

% loop set params and run decoder
clear mse
for C2NumSet_loop = 1:length(C2NumSet) ; % for each 'cell 2' population size to test
    for CorrSet_loop = 1:length(CorrSet) ; % for each noise corr value to test
        for g2set_loop = 1:length(g2set) ; % for each 'cell 2' snr to test 

            C2Num = C2NumSet(C2NumSet_loop) ; % number of "cell 2"
            Corr = CorrSet(CorrSet_loop) ; % noise correlations 
            g2 = g2set(g2set_loop) ; % gain of cell 2

            % mean response
            C1mean = stim'*g1 ; % high snr
            C2mean = stim'*g2 ; % low snr

            % noise
            C1noise = normrnd(0,sqrt(Var),nt,1) ; % noise cell 1
            for c=1:C2Num ; % for each cell 2  
                C2noise(:,c) = sqrt(Corr^2)*C1noise + sqrt(1-Corr^2)*normrnd(0,sqrt(Var),nt,1) ; % noise in cell 2
            end

            % response
            C1rsp = C1mean + C1noise ;
            for c=1:C2Num ; % for each cell 2 
                C2rsp(:,c) = C2mean + C2noise(:,c) ;
            end

            % bayes decoding

            % cov mat
            if AssumeIndFlag ;
                sigma = cov([C1noise,C2noise]).*eye(1+C2Num) ; % assume nc=0
            else
                sigma = cov([C1noise,C2noise]) ; % use correct cov
            end

            for s = 1:length(stimSet);% for each stim
                mu = [stimSet(s)*g1,ones(1,C2Num)*stimSet(s)*g2] ; % mean
                prgs(:,s) = mvnpdf([C1rsp,C2rsp],mu,sigma) ; % prob rsp given stim 1
            end

            pr = sum(prgs,2) ; % probability of the response
            
            Sest = (prgs*stimSet')./pr ; % stim estimate

            % error
            mse(C2NumSet_loop,CorrSet_loop,g2set_loop) = mean((stim'-Sest).^2) ;
            
            clear prgs C1noise C2noise C1rsp C2rsp
        end
    end
end

% figure

figure % example
plot(stim)
hold on
plot(Sest,'r--')

figure % error
subplot(3,1,1) % error as function of cell 2 gain
for C2NumSet_loop = 1:length(C2NumSet) ; % for each 'cell 2' population size to test
    for CorrSet_loop = 1:length(CorrSet) ; % for each noise corr value to test
        plot(g2set,squeeze(mse(C2NumSet_loop,CorrSet_loop,:)))
        hold on
    end
end

subplot(3,1,2) % error as function of cell 2 population size
for CorrSet_loop = 1:length(CorrSet) ; % for each noise corr value to test
    for g2set_loop = 1:length(g2set) ; % for each 'cell 2' snr to test
        plot(C2NumSet,squeeze(mse(:,CorrSet_loop,g2set_loop)))
        hold on
    end
end

subplot(3,1,3) % error as function of correlation
for C2NumSet_loop = 1:length(C2NumSet) ; % for each 'cell 2' population size to test
    for g2set_loop = 1:length(g2set) ; % for each 'cell 2' snr to test
        plot(CorrSet,squeeze(mse(C2NumSet_loop,:,g2set_loop)))
        hold on
    end
end


