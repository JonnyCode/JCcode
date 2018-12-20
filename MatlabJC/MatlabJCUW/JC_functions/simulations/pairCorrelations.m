% simple simulation of pairwise correlations: Can I extract variances below

numPnts = 1000000 ;

Noise_std_E1 = 15 ;
Noise_std_EI1 = 1 ;
Noise_std_I1 = 15 ;

Noise_std_E2 = 1 ;
Noise_std_EI2 = 0.1 ;
Noise_std_I2 = 1 ;

Noise_std_Ec = 35 ;
Noise_std_EIc = 40 ;
Noise_std_Ic = 45 ;

Thresh_E1 = 3 ;
Thresh_I1 = -10000 ;
Thresh_E2 = -10000 ;
Thresh_I2 = -0 ;


% create independant vectors
E1 = normrnd(0,Noise_std_E1,1,numPnts) ;
EI1 = normrnd(0,Noise_std_EI1,1,numPnts) ;
I1 = normrnd(0,Noise_std_I1,1,numPnts) ;

E2 = normrnd(0,Noise_std_E2,1,numPnts) ;
EI2 = normrnd(0,Noise_std_EI2,1,numPnts) ;
I2 = normrnd(0,Noise_std_I2,1,numPnts) ;

Ec = normrnd(0,Noise_std_Ec,1,numPnts) ;
EIc = normrnd(0,Noise_std_EIc,1,numPnts) ;
Ic = normrnd(0,Noise_std_Ic,1,numPnts) ;

% create exc and inh to cell pair
Et1 = E1+EI1+Ec+EIc ;
It1 = I1+EI1+Ic+EIc ;

Et2 = E2+EI2+Ec+EIc ;
It2 = I2+EI2+Ic+EIc ;

% nonlinearity
Et1(find(Et1<Thresh_E1)) = Thresh_E1 ;
It1(find(It1<Thresh_I1)) = Thresh_I1 ;
Et2(find(Et2<Thresh_E2)) = Thresh_E2 ;
It2(find(It2<Thresh_I2)) = Thresh_I2 ;


% correlations
corrE1E2 = mean(Et1.*Et2) ;
corrI1I2 = mean(It1.*It2) ;

corrE1I1 = mean(Et1.*It1) ;
corrE2I2 = mean(Et2.*It2) ;

corrE1I2 = mean(Et1.*It2) ;
corrE2I1 = mean(Et2.*It1) ;

figure
minCorrs = min([corrE1E2,corrI1I2,corrE1I1,corrE2I2,corrE1I2,corrE2I1]) ;
maxCorrs = max([corrE1E2,corrI1I2,corrE1I1,corrE2I2,corrE1I2,corrE2I1]) ;
minVar = min([Noise_std_Ec^2+Noise_std_EIc^2,corrI1I2,Noise_std_Ic^2+Noise_std_EIc^2,...
    Noise_std_EI1^2+Noise_std_EIc^2,Noise_std_EI2^2+Noise_std_EIc^2,Noise_std_EIc^2,Noise_std_EIc^2]) ;
maxVar =max([Noise_std_Ec^2+Noise_std_EIc^2,corrI1I2,Noise_std_Ic^2+Noise_std_EIc^2,...
    Noise_std_EI1^2+Noise_std_EIc^2,Noise_std_EI2^2+Noise_std_EIc^2,Noise_std_EIc^2,Noise_std_EIc^2]) ;


plot(Noise_std_Ec^2+Noise_std_EIc^2,corrE1E2,'b*')
hold on
plot(Noise_std_Ic^2+Noise_std_EIc^2,corrI1I2,'r*')

plot(Noise_std_EI1^2+Noise_std_EIc^2,corrE1I1,'y*')
plot(Noise_std_EI2^2+Noise_std_EIc^2,corrE2I2,'k*')

plot(Noise_std_EIc^2,corrE1I2,'c+')
plot(Noise_std_EIc^2,corrE2I1,'co')

line([minVar,maxVar],[minVar,maxVar])
legend('Ec','Ic','EI1','EI2','EI','IE')

xlabel('actual')
ylabel('estimated')


