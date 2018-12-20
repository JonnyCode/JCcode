% simple simulation of pairwise correlations (modified from
% pairCorreatlions.m): 

% JC 3/28/11

numPnts = 1000000 ;

Noise_std_E1 = 50 ;
%Noise_std_EI1 = zeros(1,9) ; 
Noise_std_EI1 = [0:25:200] ;
Noise_std_I1 = 50 ;

Noise_std_E2 = 50 ;
%Noise_std_EI2 = zeros(1,9) ; 
Noise_std_EI2 = [0:25:200] ;
Noise_std_I2 = 50 ;

Noise_std_Ec = 50 ;
Noise_std_EIc = [0:25:200] ;
Noise_std_Ic = 50 ;

gain_E1 = .3 ;
gain_I1 = 1 ;
gain_E2 = 1 ;
gain_I2 = .3 ;

for a=1:length(Noise_std_EIc) ;
    
    % create independant vectors
    E1 = normrnd(0,Noise_std_E1,1,numPnts) ;
    EI1 = normrnd(0,Noise_std_EI1(a),1,numPnts) ;
    I1 = normrnd(0,Noise_std_I1,1,numPnts) ;

    E2 = normrnd(0,Noise_std_E2,1,numPnts) ;
    EI2 = normrnd(0,Noise_std_EI2(a),1,numPnts) ;
    I2 = normrnd(0,Noise_std_I2,1,numPnts) ;

    Ec = normrnd(0,Noise_std_Ec,1,numPnts) ;
    EIc = normrnd(0,Noise_std_EIc(a),1,numPnts) ;
    Ic = normrnd(0,Noise_std_Ic,1,numPnts) ;

    % create exc and inh to cell pair
    Et1 = E1+EI1+Ec+EIc ;
    It1 = I1+EI1+Ic+EIc ;

    Et2 = E2+EI2+Ec+EIc ;
    It2 = I2+EI2+Ic+EIc ;

    % gain: see this as stim dependant nonlinearity
    Et1 = gain_E1*Et1 ;
    It1 = gain_I1*It1 ;
    Et2 = gain_E2*Et2 ;
    It2 = gain_I2*It2 ;

    % total synaptic current
    Current_1 = Et1-It1 ;
    Current_2 = Et2-It2 ;

    % correlations
    corrEt1Et2(a) = mean(Et1.*Et2) ;
    corrIt1It2(a) = mean(It1.*It2) ;

    corrEt1It1(a) = mean(Et1.*It1) ;
    corrEt2It2(a) = mean(Et2.*It2) ;

    corrEt1It2(a) = mean(Et1.*It2) ;
    corrEt2It1(a) = mean(Et2.*It1) ;

    corrCurrent(a) = mean(Current_1.*Current_2) ;

    % corr coef 
    corrCoefEt1Et2(a) = corrEt1Et2(a)/sqrt(mean(Et1.^2)*mean(Et2.^2)) ;
    corrCoefIt1It2(a) = corrIt1It2(a)/sqrt(mean(It1.^2)*mean(It2.^2)) ;

    corrCoefEt1It1(a) = corrEt1It1(a)/sqrt(mean(Et1.^2)*mean(It1.^2)) ; 
    corrCoefEt2It2(a) = corrEt2It2(a)/sqrt(mean(Et2.^2)*mean(It2.^2)) ; 

    corrCoefEt1It2(a) = corrEt1It2(a)/sqrt(mean(Et1.^2)*mean(It2.^2)) ; 
    corrCoefEt2It1(a) = corrEt2It1(a)/sqrt(mean(Et2.^2)*mean(It1.^2)) ; 
    
    corrCoefCurrent(a) = corrCurrent(a)/sqrt(mean(Current_1.^2)*mean(Current_2.^2)) ;

    
    % prediction of corrCurrent
    corrCurrent_Pred(a) = corrEt1Et2(a) + corrIt1It2(a) - corrEt1It2(a) - corrEt2It1(a) ;
    Denominator = sqrt((mean(Et1.^2) + mean(It1.^2) - 2*corrEt1It1(a))*(mean(Et2.^2) + mean(It2.^2) - 2*corrEt2It2(a))) ;
    corrCoefCurrent_Pred(a) = corrCurrent_Pred(a)/Denominator ;
    
end

figure
plot(Noise_std_EIc,corrCoefEt1Et2,'b')
hold on
plot(Noise_std_EIc,corrCoefIt1It2,'r')
plot(Noise_std_EIc,corrCoefEt1It2,'c')
plot(Noise_std_EIc,corrCoefEt1It1,'m')
plot(Noise_std_EIc,corrCoefCurrent,'k')
plot(Noise_std_EIc,corrCoefCurrent_Pred,'y--')
xlabel('std EIc')
ylabel('corr coef')
legend('Pee','Pii','Pei','Pie','Pss simulation','Pss analytical')

% figure
% plot(Noise_std_EIc,corrE1E2)
% 
% hold on
% plot(Noise_std_EIc,corrI1I2,'r')
% 
% figure
% plot(Noise_std_EIc,corrE1I1)
% 
% figure
% plot(Noise_std_EIc,corrE2I2)
% 
% figure
% plot(Noise_std_EIc,corrE1I2)
% 
% figure
% plot(Noise_std_EIc,corrE2I1)
% 
% figure
% plot(Noise_std_EIc,corrCurrent)
% 
% hold on
% plot(Noise_std_EIc,corrCurrent_Pred,'y--')
% 
% 
% figure
% plot(Noise_std_EIc,corrCoefE1E2)
% 
% hold on
% plot(Noise_std_EIc,corrCoefI1I2,'r')
% 
% figure
% plot(Noise_std_EIc,corrCoefE1I1)
% 
% figure
% plot(Noise_std_EIc,corrCoefE2I2)
% 
% figure
% plot(Noise_std_EIc,corrCoefE1I2)
% 
% figure
% plot(Noise_std_EIc,corrCoefE2I1)

% figure
% plot(Noise_std_EIc,corrCoefCurrent)
% 
% hold on
% plot(Noise_std_EIc,corrCoefCurrent_Pred,'y--')



