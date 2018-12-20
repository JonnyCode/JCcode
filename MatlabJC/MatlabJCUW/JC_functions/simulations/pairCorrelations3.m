% simple simulation of pairwise correlations (modified from
% pairCorreatlions2.m): 

% JC 4/28/11
% parameters
Nt_var = 100 ; % total noise variance per g per cell  
FracNc = [0:.1:1] ; % fraction of noise that orignates before gain
FracNpair = [0:.1:1] ; % fraction of common noise that is correlated accross pairs (i.e. how much overlap)

gain_E1 = 10 ;
gain_I1 = 1 ;
gain_E2 = 10 ;
gain_I2 = 1 ;

numPnts = 1000000 ;

for a = 1:length(FracNpair) ;
    for b = 1:length(FracNc) ;

        % variances
        Npair_var = Nt_var*FracNpair(a)*FracNc(b) ;
        NindCell_var = Nt_var*FracNc(b)*(1 - FracNpair(a)) ;
        NindG_var = Nt_var*(1 - FracNpair(a)*FracNc(b) - FracNc(b)*(1 - FracNpair(a))) ;

        % create independant vectors
        N1 = normrnd(0,sqrt(NindCell_var),1,numPnts) ;
        N2 = normrnd(0,sqrt(NindCell_var),1,numPnts) ;

        Np = normrnd(0,sqrt(Npair_var),1,numPnts) ;

        Ne1 = normrnd(0,sqrt(NindG_var),1,numPnts) ;
        Ne2 = normrnd(0,sqrt(NindG_var),1,numPnts) ;
        Ni1 = normrnd(0,sqrt(NindG_var),1,numPnts) ;
        Ni2 = normrnd(0,sqrt(NindG_var),1,numPnts) ;

        % create dependant vectors
        Nc1 = N1 + Np ;
        Nc2 = N2 + Np ;

        Et1 = Nc1*gain_E1 + Ne1 ;
        Et2 = Nc2*gain_E2 + Ne2 ;

        It1 = Nc1*gain_I1 + Ni1 ;
        It2 = Nc2*gain_I2 + Ni2 ;

        % total synaptic current
        Current_1 = Et1-It1 ;
        Current_2 = Et2-It2 ;

        % correlations
        corrEt1Et2(a,b) = mean(Et1.*Et2) ;
        corrIt1It2(a,b) = mean(It1.*It2) ;

        corrEt1It1(a,b) = mean(Et1.*It1) ;
        corrEt2It2(a,b) = mean(Et2.*It2) ;

        corrEt1It2(a,b) = mean(Et1.*It2) ;
        corrEt2It1(a,b) = mean(Et2.*It1) ;

        corrCurrent(a,b) = mean(Current_1.*Current_2) ;

        % corr coef 
        corrCoefEt1Et2(a,b) = corrEt1Et2(a,b)/sqrt(mean(Et1.^2)*mean(Et2.^2)) ;
        corrCoefIt1It2(a,b) = corrIt1It2(a,b)/sqrt(mean(It1.^2)*mean(It2.^2)) ;

        corrCoefEt1It1(a,b) = corrEt1It1(a,b)/sqrt(mean(Et1.^2)*mean(It1.^2)) ; 
        corrCoefEt2It2(a,b) = corrEt2It2(a,b)/sqrt(mean(Et2.^2)*mean(It2.^2)) ; 

        corrCoefEt1It2(a,b) = corrEt1It2(a,b)/sqrt(mean(Et1.^2)*mean(It2.^2)) ; 
        corrCoefEt2It1(a,b) = corrEt2It1(a,b)/sqrt(mean(Et2.^2)*mean(It1.^2)) ; 

        corrCoefCurrent(a,b) = corrCurrent(a,b)/sqrt(mean(Current_1.^2)*mean(Current_2.^2)) ;

        % prediction of corrCurrent
        corrCurrent_Pred = corrEt1Et2(a,b) + corrIt1It2(a,b) - corrEt1It2(a,b) - corrEt2It1(a,b) ;
        Denominator = sqrt((mean(Et1.^2) + mean(It1.^2) - 2*corrEt1It1(a,b))*(mean(Et2.^2) + mean(It2.^2) - 2*corrEt2It2(a,b))) ;
        corrCoefCurrent_Pred(a,b) = corrCurrent_Pred/Denominator ;
    end 
end

figure
plot(FracNc,corrCoefEt1Et2(5,:),'b')
hold on
plot(FracNc,corrCoefIt1It2(5,:),'r')
plot(FracNc,corrCoefEt1It2(5,:),'c--')
plot(FracNc,corrCoefEt2It1(5,:),'m:')
plot(FracNc,corrCoefCurrent(5,:),'k')
plot(FracNc,corrCoefCurrent_Pred(5,:),'y--')
plot(FracNc,corrCoefEt1It1(5,:),'c:')
plot(FracNc,corrCoefEt2It2(5,:),'m--')
xlabel('Fraction of total noise originating before gain')
ylabel('corr coef')
legend('Pee','Pii','Pei','Pie','Pss simulation','Pss analytical','Cei1','Cei2')

figure
plot(FracNpair ,corrCoefEt1Et2(:,5),'b')
hold on
plot(FracNpair ,corrCoefIt1It2(:,5),'r')
plot(FracNpair ,corrCoefEt1It2(:,5),'c--')
plot(FracNpair ,corrCoefEt2It1(:,5),'m:')
plot(FracNpair ,corrCoefCurrent(:,5),'k')
plot(FracNpair ,corrCoefCurrent_Pred(:,5),'y--')
plot(FracNpair ,corrCoefEt1It1(:,5),'c:')
plot(FracNpair ,corrCoefEt2It2(:,5),'m--')
xlabel('Fraction of original noise that is shared between cells')
ylabel('corr coef')
legend('Pee','Pii','Pei','Pie','Pss simulation','Pss analytical','Cei1','Cei2')

figure
pcolor(FracNpair,FracNc,corrCoefCurrent)
colorbar
xlabel('Fraction of total noise originating before gain')
ylabel('Fraction of original noise that is shared between cells')
