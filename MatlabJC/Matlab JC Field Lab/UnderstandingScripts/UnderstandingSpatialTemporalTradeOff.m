% simple simulation to understand temporal-spatial integration tradeoffs

% JC 8-7-2018 

% simulation #1 - benefits of space

% params
sampleRate = 1000 ; % time points/sec
Duration = 10 ; % sec
nu = 50 ; % number of units
nsp = 50 ; % number of smooth points
Sfl = 10 ; % number of points to smooth S
TintRange = [1:2:100] ; % smooth only works on odd numbers

% from params
ntp = sampleRate*Duration ; % number time points

% stim (time)
Swhite = normrnd(0,10,1,ntp) ; % global signal (white)
S = smooth(Swhite,Sfl)' ; % smooth
Sstd = std(S) ;

% noise (space, time)
for i = 1:nu ; % for each unit
    n = normrnd(0,Sstd,nu,ntp) ; % white noise
end

% power spectrum
Srs = reshape(S',[sampleRate,Duration])' ;
nrs = n(:,1:sampleRate) ;
[px,ps_S] = PowerSpectrumFinder(Srs,sampleRate) ;
[px,ps_n] = PowerSpectrumFinder(nrs,sampleRate) ;
ps_Snr = ps_S./ps_n ;

% temporal integration
for Tint = 1:length(TintRange) ;
    Stint(Tint,:) = smooth(S,TintRange(Tint)) ; 
    ntint(Tint,:) = smooth(n(1,:),TintRange(Tint)) ;
    StintReshape = reshape(Stint(Tint,:)',[sampleRate,Duration])' ;
    [px,ps_Stint(Tint,:)] = PowerSpectrumFinder(StintReshape,sampleRate) ;
    mseTint(Tint) = mean((S-(Stint(Tint,:)+ntint(Tint,:))).^2) ; % mean square error
end

% spatial integration
for Sint = [1:nu] ;
    nsint(Sint,:) = mean(n(1:Sint,:),1) ;
    mseSint(Sint) = mean((S-(S+nsint(Sint,:))).^2) ; % mean square error
end

% mean squared error match
[mseTintMin,mseTintMini] = min([mseTint]) ; % minumum mse for temporal integration
[~,mseSintMatchi] = min(abs([mseSint - mseTint(mseTintMini)])) ; % most similar mse for spatial integration

% figure
figure 
subplot(1,3,1)
plot([1:nu],mseSint,'b')
hold on
plot(TintRange,mseTint,'r')
plot(mseSintMatchi,mseSint(mseSintMatchi),'bo')
plot(TintRange(mseTintMini),mseTintMin,'ro')

subplot(1,3,2)
plot(S+nsint(mseSintMatchi,:),'b')
hold on
plot(Stint(mseTintMini,:)+ntint(mseTintMini,:),'r')
plot(S,'k')

subplot(1,3,3)
loglog(px,ps_S,'b')
hold on
loglog(px,ps_Stint(mseTintMini,:),'r')

hgsave('GenericSimulationFig')
print(gcf,'-dpdf','GenericSimulationFig')


% simulation #2 - space and time

X = 100 ;
Y = 100 ;
T = 100 ;
StdSpatialBlur = 20 ; % spatial points
StdTemporalBlur = 20 ; % time points
SfRange = [0:1:10] ;
TfRange = [0:1:10] ;

S = normrnd(0,10,X,Y,T) ; % movie

Sblur = MovBlurring(S,StdSpatialBlur,StdTemporalBlur) ; % spatial/temporal blur
Sblur = Sblur - mean(Sblur(:)) ; % zero mean
Sblur = Sblur/std(Sblur(:)) ; % make std = 1

n = normrnd(0,1,X,Y,T) ; % make noise (std=1)

R = Sblur+n ; % add signal and noise

% integrate spatially and temporally
Mse = nan(length(SfRange),length(TfRange)) ;
for Sf = 1:length(SfRange) ;
    for Tf = 1:length(TfRange) ;        
        RblurTemp = MovBlurring(R,SfRange(Sf),TfRange(Tf)) ; % spatial/temporal blur (attempting to integrate out noise)
        RblurTemp = RblurTemp - mean(RblurTemp(:)) ; % zero mean
        %Rblur{Sf,Tf} = RblurTemp/std(RblurTemp(:)) ; % make std =1 (if you want to keep Rblur for plotting)
        Rblur= RblurTemp/std(RblurTemp(:)) ; % make std =1
        %Temp = (Sblur - Rblur{Sf,Tf}).^2 ; % error (if you want to keep Rblur for plotting)
        Temp = (Sblur - Rblur).^2 ; % error
        Mse(Sf,Tf) = mean(Temp(:)) ;
        
        % power spectrum of signal and noise - temporal
%         [TimePsx,TimePs] = PowerSpectrumFinder(reshape(Sblur-Rblur,[X*Y,T]),1) ;
%         figure(1)
%         subplot(1,length(TfRange),Tf)
%         loglog(TimePsx,TimePs)
%         hold on
        
        % nonlinearity
%         [nlX,nlY,Nonlinearity_sem] = NonLinFilterFinder(Rblur,Sblur,.1)  ; % nl
%         
%         figure(1)
%         subplot(1,length(SfRange),Sf)
%         plot(nlX,nlY)
%         hold on
%         
%         figure(2)
%         subplot(1,length(TfRange),Tf)
%         plot(nlX,nlY)
%         hold on
        
        disp([Sf,Tf])
        %pause
    end
end

% figure
figure
subplot(3,1,1)
imagesc(SfRange,TfRange,flipud(log(Mse)))
colorbar
xlabel('Tf std')
ylabel('Sf std')

subplot(3,1,2)
for Sf = 1:20:length(SfRange) ; 
    plot(TfRange,log(Mse(Sf,:)),'k','LineWidth',5/(log10(Sf)+1))
    hold on
end
xlabel('Tf std')
ylabel('mse')
xlim([.5,10])

subplot(3,1,3)
for Tf = 1:20:length(TfRange) ;
    plot(SfRange,log(Mse(:,Tf)),'k','LineWidth',5/(log10(Tf)+1))
    hold on
end
xlabel('Sf std')
ylabel('mse')
xlim([.5,10])

saveas(gcf,'SimpleSim2Fig1')
print(gcf, '-dpdf','SimpleSim2Fig1')

figure % signal frames
egFrames = [1:20:100] ;
FilterEg = [20,20] ; 
for a = 1:length(egFrames) ;
    subplot(3,length(egFrames),a)
    imagesc(Sblur(:,:,egFrames(a)))
    colormap gray
    
    subplot(3,length(egFrames),a+length(egFrames))
    imagesc(R(:,:,egFrames(a)))
    colormap gray
    
    subplot(3,length(egFrames),a+length(egFrames)*2)
    temp = Rblur{FilterEg} ;
    imagesc(temp(:,:,egFrames(a)))
    colormap gray
end
saveas(gcf,'SimpleSim2Fig2')
print(gcf, '-dpdf','SimpleSim2Fig2')


figure % filter
egFilterPulse = zeros(20,20,20) ;
egFilterPulse(10,10,10) = 1 ;
egFilter = MovBlurring(egFilterPulse,5,2) ;
subplot(2,1,1) ; % spatial filter
imagesc(egFilter(:,:,10))
%colormap(brewermap([],'RdBu'))

subplot(2,1,2) ; % temporal filter
imagesc(squeeze(egFilter(10,:,:))')   
%colormap(brewermap([],'RdBu'))
saveas(gcf,'SimpleSim2Fig3')
print(gcf, '-dpdf','SimpleSim2Fig3')

figure ; % gaussian
plot(pdf('norm',[-10:.1:10],0,3),'k')
print(gcf, '-dpdf','SimpleSim2Fig4')

% 

