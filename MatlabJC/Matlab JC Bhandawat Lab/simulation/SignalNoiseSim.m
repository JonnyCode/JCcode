% signal and noise simulation

Stim = [1:.1:10] ;
RspHalfMax = 5 ;
RspStd = 1 ;
RspMax = 10 ;

BgStd = 2 ;
NumReceptors = 40 ;
gainConstant = 5 ;

RspMean =  cumsum(exp(-((Stim-RspHalfMax).^2)/(2*RspStd^2))) ;
RspMean = RspMax*RspMean/max(RspMean) ;

BgStd2 = BgStd/sqrt(NumReceptors) ; % averaging

BgStd3 = BgStd2 + BgStd ; % added noise

BgStd4 = BgStd2*gainConstant ; % multiplication

BgStd5 = BgStd4 + BgStd ; % added noise

RspMean2 = RspMean*gainConstant ; % multiplication


figure
subplot(2,2,1)
plot(Stim,RspMean)
hold on
plot(Stim,Stim*0+BgStd,'--')
plot(Stim,Stim*0+BgStd2,'r--')

subplot(2,2,2)
plot(Stim,RspMean)
hold on
plot(Stim,Stim*0+BgStd3,'r--')

subplot(2,2,3)
plot(Stim,RspMean2,'g')
hold on
plot(Stim,Stim*0+BgStd4,'g--')

subplot(2,2,4)
plot(Stim,RspMean2,'g')
hold on
plot(Stim,Stim*0+BgStd5,'g--')
plot(Stim,RspMean)
plot(Stim,Stim*0+BgStd,'--')
