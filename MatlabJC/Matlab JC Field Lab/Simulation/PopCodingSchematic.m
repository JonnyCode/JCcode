% simple simulation to illustrate pop coding perpspective

time = [0:.001:20] ; % time
lt = length(time) ;

filter = pdf('norm',time,10,.02) ;

% imagine two unlcorrelated neurons with stimuli responses defined by blurred pulses

c0p = normrnd(0,1,1,lt) ;
c0p(c0p<3) = 0 ;
c0p(c0p>0) = 1 ;

c1p = abs(c0p.*normrnd(1,5,1,lt)) ;
c2p = abs(c0p.*normrnd(1,5,1,lt)) ;

c1 = conv(c1p,filter,'same') ;
c2 = conv(c2p,filter,'same') ;
c0 = conv(c0p,filter,'same') ;

% figure
figure
subplot(4,1,1)
plot(time,c1p,'k')
hold on
plot(time,c2p,'b')
plot(time,c0p,'r')

subplot(4,1,2)
plot(c1p,c2p,'b*')
hold on
plot(c1p,c0p,'r*')

subplot(4,1,3)
plot(time,c1,'k')
hold on
plot(time,c2,'b')
plot(time,c0,'r')

subplot(4,1,4)
plot(xcov(c1,c2,'coef'),'b')
hold on
plot(xcov(c1,c0,'coef'),'r')

