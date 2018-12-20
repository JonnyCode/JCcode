% this code will extract some basic parameters from isolated current steps
% JC 4/28/08

% Set cell and trace parameters
CellInfo_str = '052908Bc1_a' ;

Epochs_str_exc = {'[326:8:405]','[327:8:405]','[328:8:405]','[329:8:405]','[330:8:405]','[331:8:405]','[332:8:405]','[333:8:405]'} ;
Epochs_str_inh = {'[426:8:505]','[427:8:505]','[428:8:505]','[429:8:505]','[430:8:505]','[431:8:505]','[432:8:505]','[433:8:505]'} ;
Epochs_str = [Epochs_str_exc, Epochs_str_inh] ;

EpochCondNum_exc = [ones(1,8)*2] ;
EpochCondNum_inh = [ones(1,8)*2] ;
EpochCondNum = [EpochCondNum_exc,EpochCondNum_inh] ;

Epochs_lightAmp_exc = [.125,-.125,.25,-.25,.5,-.5,1,-1] ;
Epochs_lightAmp_inh = [.125,-.125,.25,-.25,.5,-.5,1,-1] ;
Epochs_lightAmp = [Epochs_lightAmp_exc, Epochs_lightAmp_inh] ; 

StartPnt = 5000 ; % of light step
EndPnt = 10000 ;

% Get appropriate cell structure format
cd ~/data_analysis/Index ;
load (CellInfo_str) ;
CellInfo = LoadSCIData(CellInfo,1) ;        
EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ; 

for a = 1:length(Epochs_str) ;          % for each set exc epoch

Epochs{a} = str2num(Epochs_str{a}) ; % change string to number

for b=1:length(Epochs{a}) ; % for each epoch in that set
Idata{a}(b)= find(EpochCondition(EpochCondNum(b)).EpochNumbers == Epochs{a}(b)) ; % find the index of those epochs
end
clear b

I{a} = EpochCondition(EpochCondNum(a)).EpochData.Data(Idata{a},:) ; % the individual epochs

Imean(a,:) = mean(I{a},1) ; % mean current (pA)

if a>length(Epochs_str_exc) ; % if the epoch is an inhibitory epoch 
G{a} = I{a}./60 ;     % divide by this driving force
else 
G{a} = I{a}./-60 ;    % otherwise the epoch is exc so use this driving force
end

Gmean(a,:) = mean(G{a},1) ;
Gvar(a,:) = var(G{a},1) ; % time dependent variance 
Gstepvar_mean(a) = mean(Gvar(a,StartPnt:EndPnt)) ; % mean variance during the light step

% find the peak of each trace (either the max if its primarily an increase
% or the min if its primarily a decrease)
if a<=length(Epochs_str_exc) & Epochs_lightAmp(a)<0 ; % if the epoch is an excitatory decrement epoch set... 
[Gpeak{a},Gpeakpnt{a}] = min(G{a}(:,StartPnt:EndPnt),[],2) ; % find the lowest point of each individual conductance
[Gmean_peak(a),Gmean_peakpnt(a)] = min(Gmean(a,StartPnt:EndPnt),[],2) ; % find the lowest point of the mean conductance
else % if its any other trace...    
[Gpeak{a},Gpeakpnt{a}] = max(G{a}(:,StartPnt:EndPnt),[],2) ; % find the highest point of each individual conductance
[Gmean_peak(a),Gmean_peakpnt(a)] = max(Gmean(a,StartPnt:EndPnt),[],2) ; % find the lowest point of the mean conductance
end

Gpeakpnt{a} = Gpeakpnt{a} + StartPnt -1 ; % this needs to happen in order to keep everything on the same trace bc we've narrowed the window above 
Gmean_peakpnt(a) = Gmean_peakpnt(a) + StartPnt -1 ; 

% find the points on each trace that span 20-80% of the peak  
if a<=length(Epochs_str_exc) & Epochs_lightAmp(a)<0  ; % if the epoch set is an excitatory decrement epoch set... 
for b = 1:length(Epochs{a}) ; % for each individual epoch in the set...
GuppercrossPnt{a}(b) = max(find(G{a}(b,1:Gpeakpnt{a}(b))>.8*Gpeak{a}(b),1,'last'), StartPnt) ;
GlowercrossPnt{a}(b) = max(find(G{a}(b,1:Gpeakpnt{a}(b))>.2*Gpeak{a}(b),1,'last'), StartPnt) ;  
GreturnPnt{a}(b) = min(find(G{a}(b,Gpeakpnt{a}(b):end)>.2*Gpeak{a}(b),1,'first'),length(Gmean)) + Gpeakpnt{a}(b)-1 ;
end
clear b
Gmean_uppercrossPnt(a) = max(find(Gmean(a,1:Gmean_peakpnt(a))>.8*Gmean_peak(a),1,'last'), StartPnt) ;
Gmean_lowercrossPnt(a) = max(find(Gmean(a,1:Gmean_peakpnt(a))>.2*Gmean_peak(a),1,'last'), StartPnt) ;
Gmean_returnPnt(a) = min(find(Gmean(a,Gmean_peakpnt(a):end)>.2*Gmean_peak(a),1,'first'),length(Gmean)) + Gmean_peakpnt(a)-1 ;

else
for b = 1:length(Epochs{a}) ; % for each individual epoch in the set...
GuppercrossPnt{a}(b) = max(find(G{a}(b,1:Gpeakpnt{a}(b))<.8*Gpeak{a}(b),1,'last'), StartPnt) ; % find the last point that is between the peak and crossing of the upper threshold or use the start point 
GlowercrossPnt{a}(b) = max(find(G{a}(b,1:Gpeakpnt{a}(b))<.2*Gpeak{a}(b),1,'last'), StartPnt) ; % find the last point that is between the peak and crossing the lower threshold or use the start point
GreturnPnt{a}(b) = min(find(G{a}(b,Gpeakpnt{a}(b):end)<.2*Gpeak{a}(b),1,'first'),length(Gmean)) + Gpeakpnt{a}(b)-1 ;
end
clear b
Gmean_uppercrossPnt(a) = max(find(Gmean(a,1:Gmean_peakpnt(a))<.8*Gmean_peak(a),1,'last'), StartPnt) ;
Gmean_lowercrossPnt(a) = max(find(Gmean(a,1:Gmean_peakpnt(a))<.2*Gmean_peak(a),1,'last'), StartPnt) ;
Gmean_returnPnt(a) = min(find(Gmean(a,Gmean_peakpnt(a):end)<.2*Gmean_peak(a),1,'first'),length(Gmean)) + Gmean_peakpnt(a)-1 ;
end

GriseTime{a} = GuppercrossPnt{a} - GlowercrossPnt{a} ;
Gmean_riseTime(a) = Gmean_uppercrossPnt(a) - Gmean_lowercrossPnt(a) ;

Gduration{a} = GreturnPnt{a} - GlowercrossPnt{a} ;
Gmean_duration(a) = Gmean_returnPnt(a) - Gmean_lowercrossPnt(a) ;

Gpeak_mean(a) = mean(Gpeak{a}) ;
Gpeak_var(a) = var(Gpeak{a}) ;

Gpeakpnt_mean(a) = mean(Gpeakpnt{a}) ;
Gpeakpnt_var(a) = var(Gpeakpnt{a}) ;

GlowercrossPnt_mean(a) = mean(GlowercrossPnt{a}) ;
GlowercrossPnt_var(a) = var(GlowercrossPnt{a}) ;

GriseTime_mean(a) = mean(GriseTime{a}) ;
GriseTime_var(a) = var(GriseTime{a}) ;

Gduration_mean(a) = mean(Gduration{a}) ;
Gduration_var(a) = var(Gduration{a}) ;

% find the sum of the conductance during the light step
Gsum{a} = sum(G{a}(:,StartPnt:EndPnt),2) ;
Gmean_sum(a) = sum(Gmean(a,StartPnt:EndPnt)) ;

Gsum_mean(a) =  mean(Gsum{a}) ;
Gsum_var(a) = var(Gsum{a}) ;

end % end epoch set loop (a loop)
clear a


% FIGURES
figure
title(CellInfo_str)
text('units','normalized') %this sets the text in units where (1,1) is upper right
text(.25, .5,['exc=' Epochs_str_exc],'units','normalized')
text(.1, .5,[num2str(Epochs_lightAmp_exc')],'units','normalized')
text(.6, .5,['inh=' Epochs_str_inh],'units','normalized')
text(.75, .5,[num2str(Epochs_lightAmp_inh')],'units','normalized')

% plot mean conductances with peak, upper. lower, and return cross points marked
figure 
for a = 1:length(Epochs_str) 
if a<=length(Epochs_str_exc) & Epochs_lightAmp(a)<0  ;      % if the epoch set is an excitatory decrement epoch set... 
    subplot(2,2,1)
    title('Gmean exc decs')
elseif a<=length(Epochs_str_exc) & Epochs_lightAmp(a)>0 ;   % exc inc 
    subplot(2,2,2)
     title('Gmean exc incs')
elseif a>length(Epochs_str_exc) & Epochs_lightAmp(a)<0 ;    % inh dec
    subplot(2,2,3)
     title('Gmean inh decs')
else                                                        % inh inc
    subplot(2,2,4)
     title('Gmean inh incs')
end
    plot(Gmean(a,:))
    hold on
    plot(Gmean_peakpnt(a), Gmean_peak(a), 'k*')
    plot(Gmean_lowercrossPnt(a), Gmean(a,Gmean_lowercrossPnt(a)),'g*')
    plot(Gmean_uppercrossPnt(a), Gmean(a,Gmean_uppercrossPnt(a)),'r*')
    plot(Gmean_returnPnt(a),Gmean(a,Gmean_returnPnt(a)),'y*')
end 
clear a

% plot figure comparing OFF and ON inhibtion measurements from incs and
% decs
decInhEpochs = find(Epochs_lightAmp_inh<0) + length(Epochs_lightAmp_exc) ;
incInhEpochs = find(Epochs_lightAmp_inh>0) + length(Epochs_lightAmp_exc) ;

figure
subplot(3,2,1) % Gmean peaks
plot(abs(Epochs_lightAmp(decInhEpochs)), Gmean_peak(decInhEpochs),'b*')
hold on
plot(Epochs_lightAmp(incInhEpochs), Gmean_peak(incInhEpochs),'r*')
title('Gmean inh peaks')
xlabel('Abs contrast')
ylabel('Gmean peak (nS)')
legend('decs','incs')

subplot(3,2,2) % ratio of Gmean peaks
plot(Epochs_lightAmp(incInhEpochs), Gmean_peak(decInhEpochs)./Gmean_peak(incInhEpochs),'b-*') ;
title('Ratio Gmean inh peaks')
xlabel('Abs contrast')
ylabel('Decs/incs Gmean peak (nS)')

subplot(3,2,3) % individual peaks
for a = length(Epochs_str_exc)+1:length(Epochs_str) % for each set of inh epochs
if sum(ismember(decInhEpochs,a)) % if the epoch set is an inh decrement  
    plot(repmat(abs(Epochs_lightAmp(a)),size(Gpeak{a})), Gpeak{a},'b*')
    hold on
    plot(abs(Epochs_lightAmp(a)),Gpeak_mean(a),'b+')
else
    plot(repmat(abs(Epochs_lightAmp(a)),size(Gpeak{a})), Gpeak{a},'ro')
    hold on
    plot(abs(Epochs_lightAmp(a)),Gpeak_mean(a),'r+')
end
end
title('G inh peaks')
xlabel('Abs contrast')
ylabel('G peak (nS)')

subplot(3,2,4) % ratio of mean of individual peaks
plot(Epochs_lightAmp(incInhEpochs), Gpeak_mean(decInhEpochs)./Gpeak_mean(incInhEpochs),'b-*') ;
title('Ratio mean G inh peaks')
xlabel('Abs contrast')
ylabel('Decs/incs mean G peak (nS)')

subplot(3,2,5) % mean and SEM of individual peaks
for a = length(Epochs_str_exc)+1:length(Epochs_str) % for each set of inh epochs
if sum(ismember(decInhEpochs,a)) % if the epoch set is an inh decrement  
    errorbar(abs(Epochs_lightAmp(a)),Gpeak_mean(a),sqrt(Gpeak_var(a))/length(Gpeak{a}),'b*')    
else
    errorbar(abs(Epochs_lightAmp(a)),Gpeak_mean(a),sqrt(Gpeak_var(a))/length(Gpeak{a}),'ro') 
end
hold on
end
title('G inh peaks')
xlabel('Abs contrast')
ylabel('G peak (nS)')

subplot(3,2,6) % mean vs. variance of indiivdual peaks
plot(Gpeak_mean(decInhEpochs), Gpeak_var(decInhEpochs),'b*')
hold on
plot(Gpeak_mean(incInhEpochs), Gpeak_var(incInhEpochs),'r*')
title('G peaks mean vs. var')
xlabel('mean individ G peak (nS)')
ylabel('var individ G peak (nS^2)')

% plot as figure above but use sum of conductance instead of peak
figure
subplot(4,2,1) % Gmean sums
plot(abs(Epochs_lightAmp(decInhEpochs)), Gmean_sum(decInhEpochs),'b*')
hold on
plot(Epochs_lightAmp(incInhEpochs), Gmean_sum(incInhEpochs),'r*')
title('Gmean inh sum')
xlabel('Abs contrast')
ylabel('Gmean sum (nS)')
legend('decs','incs')

subplot(4,2,2) % ratio of Gmean sums
plot(Epochs_lightAmp(incInhEpochs), Gmean_sum(decInhEpochs)./Gmean_sum(incInhEpochs),'b-*') ;
title('Ratio Gmean inh sum')
xlabel('Abs contrast')
ylabel('Decs/incs Gmean sum (nS)')

subplot(4,2,3) % individual sums
for a = length(Epochs_str_exc)+1:length(Epochs_str) % for each set of inh epochs
if sum(ismember(decInhEpochs,a)) % if the epoch set is an inh decrement  
    plot(repmat(abs(Epochs_lightAmp(a)),size(Gpeak{a})), Gsum{a},'b*')
    hold on
    plot(abs(Epochs_lightAmp(a)),Gsum_mean(a),'b+')
else
    plot(repmat(abs(Epochs_lightAmp(a)),size(Gsum{a})), Gsum{a},'ro')
    hold on
    plot(abs(Epochs_lightAmp(a)),Gsum_mean(a),'r+')
end
end
title('individ G inh sum')
xlabel('Abs contrast')
ylabel('G sum (nS)')

subplot(4,2,4) % ratio of mean of individual sums
plot(Epochs_lightAmp(incInhEpochs), Gsum_mean(decInhEpochs)./Gsum_mean(incInhEpochs),'b-*') ;
title('Ratio mean G inh sum')
xlabel('Abs contrast')
ylabel('Decs/incs mean G sum (nS)')

subplot(4,2,5) % mean and SEM of individual sums
for a = length(Epochs_str_exc)+1:length(Epochs_str) % for each set of inh epochs
if sum(ismember(decInhEpochs,a)) % if the epoch set is an inh decrement  
    errorbar(abs(Epochs_lightAmp(a)),Gsum_mean(a),sqrt(Gsum_var(a))/length(Gsum{a}),'b*')    
else
    errorbar(abs(Epochs_lightAmp(a)),Gsum_mean(a),sqrt(Gsum_var(a))/length(Gsum{a}),'ro') 
end
hold on
end
title('G inh sum')
xlabel('Abs contrast')
ylabel('G sum (nS)')

subplot(4,2,6) % mean vs. variance of indiivdual sums
plot(Gsum_mean(decInhEpochs), Gsum_var(decInhEpochs),'b*')
hold on
plot(Gsum_mean(incInhEpochs), Gsum_var(incInhEpochs),'r*')
title('G sum mean vs. var')
xlabel('mean individ G sum')
ylabel('var individ G sum')

subplot(4,2,7) % plot mean variance of the individual inh traces
plot(abs(Epochs_lightAmp(decInhEpochs)), Gstepvar_mean(decInhEpochs),'b*')
hold on
plot(Epochs_lightAmp(incInhEpochs), Gstepvar_mean(incInhEpochs),'r*')
title('Gstep var mean')
xlabel('abs contrast')
ylabel('mean var of G')

subplot(4,2,8) % plot mean variance of the individual inh traces divided by the mean sum
plot(abs(Epochs_lightAmp(decInhEpochs)), Gstepvar_mean(decInhEpochs)./mean(Gmean(decInhEpochs,StartPnt:EndPnt),2)','b*')
hold on
plot(Epochs_lightAmp(incInhEpochs), Gstepvar_mean(incInhEpochs)./mean(Gmean(incInhEpochs,StartPnt:EndPnt),2)','r*')
title('Gstep var mean/stepmean')
xlabel('abs contrast')
ylabel('mean var of G/ mean G')

% Plot comparison of rise time between OFF and ON inhibition
figure
subplot(3,2,1) % plot risetime of Gmean
plot(abs(Epochs_lightAmp(decInhEpochs)), Gmean_riseTime(decInhEpochs),'b*')
hold on
plot(Epochs_lightAmp(incInhEpochs), Gmean_riseTime(incInhEpochs),'r*')
title('Gmean inh rise time')
xlabel('Abs contrast')
ylabel('Gmean risetime (sample points)')
legend('decs','incs')

subplot(3,2,3) % plot rise time of individual traces
for a = length(Epochs_str_exc)+1:length(Epochs_str) % for each set of inh epochs
if sum(ismember(decInhEpochs,a)) % if the epoch set is an inh decrement  
    plot(repmat(abs(Epochs_lightAmp(a)),size(GriseTime{a})), GriseTime{a},'b*')
    hold on
    plot(abs(Epochs_lightAmp(a)),GriseTime_mean(a),'b+')
else
    plot(repmat(abs(Epochs_lightAmp(a)),size(GriseTime{a})), GriseTime{a},'ro')
    hold on
    plot(abs(Epochs_lightAmp(a)),GriseTime_mean(a),'r+')
end
end
title('G inh rise time')
xlabel('Abs contrast')
ylabel('G risetime (points)')

subplot(3,2,5) % plot rise time of mean of individs with SEM
for a = length(Epochs_str_exc)+1:length(Epochs_str) % for each set of inh epochs
if sum(ismember(decInhEpochs,a)) % if the epoch set is an inh decrement  
    errorbar(abs(Epochs_lightAmp(a)),GriseTime_mean(a),sqrt(GriseTime_var(a))/length(GriseTime{a}),'b*')    
else
    errorbar(abs(Epochs_lightAmp(a)),GriseTime_mean(a),sqrt(GriseTime_var(a))/length(GriseTime{a}),'ro') 
end
hold on
end
title('G inh rise time')
xlabel('Abs contrast')
ylabel('G rise time (points)')


% plot to compare measures of excitiation vs. inhibition 
figure
subplot(3,2,1) % plot Gmean peaks
plot(Epochs_lightAmp_exc,Gmean_peak(1:length(Epochs_lightAmp_exc)),'b*')
hold on
plot(Epochs_lightAmp_inh,Gmean_peak(length(Epochs_lightAmp_exc)+1:end),'r*')
title('peak Gmean')
xlabel('step contrast')
ylabel('peak of Gmean (nS)')
legend('exc','inh')

subplot(3,2,2) % plot ratio of exc/ inh peak
plot(Epochs_lightAmp_exc,Gmean_peak(1:length(Epochs_lightAmp_exc))./Gmean_peak(length(Epochs_lightAmp_exc)+1:end),'b*')
title('ratio peak Gmean exc/inh ')
xlabel('step contrast')
ylabel('peak of Gmean exc/inh (nS)')

subplot(3,2,3) % plot individual G peaks
for a = 1:length(Epochs_str)
    if a<=length(Epochs_str_exc)
    plot(repmat(Epochs_lightAmp(a),size(Gpeak{a})),Gpeak{a},'b*')
    hold on
    plot(Epochs_lightAmp(a),Gpeak_mean(a),'+')    
    else 
    plot(repmat(Epochs_lightAmp(a),size(Gpeak{a})),Gpeak{a},'ro')
    hold on
    plot(Epochs_lightAmp(a),Gpeak_mean(a),'+r')        
    end
end
title('individual G peaks')
xlabel('step contrast')
ylabel('peak of G (nS)')

subplot(3,2,4) % plot ratio of mean exc/inh peak
plot(Epochs_lightAmp_exc,Gpeak_mean(1:length(Epochs_lightAmp_exc))./Gpeak_mean(length(Epochs_lightAmp_exc)+1:end),'b*')
title('ratio peak Gpeak mean exc/inh ')
xlabel('step contrast')
ylabel('peak of mean indivdl G exc/inh (nS)')

subplot(3,2,5) % plot mean G peak w/ SEM error bars
for a = 1:length(Epochs_str) % for each set of epochs
if a<=length(Epochs_str_exc)
    errorbar(Epochs_lightAmp(a), Gpeak_mean(a), sqrt(Gpeak_var(a))/length(Gpeak{a}),'b*') 
else
    errorbar(Epochs_lightAmp(a), Gpeak_mean(a), sqrt(Gpeak_var(a))/length(Gpeak{a}),'ro') 
end
hold on
end
title('individual G peaks mean w/ SEM')
xlabel('step contrast')
ylabel('peak of G (nS)')

subplot(3,2,6) % plot mean of G peak vs var of G peak
for a = 1:length(Epochs_str) % for each set of epochs
if a<=length(Epochs_str_exc)
    plot(Gpeak_mean(a), Gpeak_var(a),'b*') 
else
    plot(Gpeak_mean(a), Gpeak_var(a),'r*') 
end
hold on
end
title('Gpeak mean vs variance')
xlabel('mean individ G peak (nS)')
ylabel('var individ G peak (nS^2)')


% as figure above but with sum G instead of peak G 
figure
subplot(4,2,1) % plot Gmean sum
plot(Epochs_lightAmp_exc,Gmean_sum(1:length(Epochs_lightAmp_exc)),'b*')
hold on
plot(Epochs_lightAmp_inh,Gmean_sum(length(Epochs_lightAmp_exc)+1:end),'r*')
title('sum Gmean')
xlabel('step contrast')
ylabel('sum of Gmean (nS)')
legend('exc','inh')

subplot(4,2,2) % plot ratio of exc/ inh sum
plot(Epochs_lightAmp_exc,Gmean_sum(1:length(Epochs_lightAmp_exc))./Gmean_sum(length(Epochs_lightAmp_exc)+1:end),'b*')
title('ratio sum Gmean exc/inh ')
xlabel('step contrast')
ylabel('sum of Gmean exc/inh (nS)')

subplot(4,2,3) % plot individual G sum
for a = 1:length(Epochs_str)
    if a<=length(Epochs_str_exc)
    plot(repmat(Epochs_lightAmp(a),size(Gsum{a})),Gsum{a},'b*')
    hold on
    plot(Epochs_lightAmp(a),Gsum_mean(a),'+')    
    else 
    plot(repmat(Epochs_lightAmp(a),size(Gsum{a})),Gsum{a},'ro')
    hold on
    plot(Epochs_lightAmp(a),Gsum_mean(a),'+r')        
    end
end
title('individual G sum')
xlabel('step contrast')
ylabel('sum of G (nS)')

subplot(4,2,4) % plot ratio of mean exc/inh sum
plot(Epochs_lightAmp_exc,Gsum_mean(1:length(Epochs_lightAmp_exc))./Gsum_mean(length(Epochs_lightAmp_exc)+1:end),'b*')
title('ratio peak Gsum mean exc/inh ')
xlabel('step contrast')
ylabel('sum of mean indivdl G exc/inh (nS)')

subplot(4,2,5) % plot mean G sum w/ SEM error bars
for a = 1:length(Epochs_str) % for each set of epochs
if a<=length(Epochs_str_exc)
    errorbar(Epochs_lightAmp(a), Gsum_mean(a), sqrt(Gsum_var(a))/length(Gsum{a}),'b*') 
else
    errorbar(Epochs_lightAmp(a), Gsum_mean(a), sqrt(Gsum_var(a))/length(Gsum{a}),'ro') 
end
hold on
end
title('individual G sum mean w/ SEM')
xlabel('step contrast')
ylabel('peak of G (nS)')

subplot(4,2,6) % plot mean of G sum vs var of G sum
for a = 1:length(Epochs_str) % for each set of epochs
if a<=length(Epochs_str_exc)
    plot(Gsum_mean(a), Gsum_var(a),'b*') 
else
    plot(Gsum_mean(a), Gsum_var(a),'r*') 
end
hold on
end
title('Gsum mean vs variance')
xlabel('mean individ G sum (nS)')
ylabel('var individ G sum (nS^2)')

subplot(4,2,7) % plot mean variance of G step as function of contrast
plot(Epochs_lightAmp_exc,Gstepvar_mean(1:length(Epochs_lightAmp_exc)),'b*')
hold on
plot(Epochs_lightAmp_inh,Gstepvar_mean(length(Epochs_lightAmp_exc)+1:end),'r*')
title('Gstep var mean')
xlabel('contrast')
ylabel('mean var of G')

subplot(4,2,8) % plot mean variance of the individual traces divided by the mean sum
plot(Epochs_lightAmp_exc, Gstepvar_mean(1:length(Epochs_lightAmp_exc))./mean(Gmean(1:length(Epochs_lightAmp_exc),StartPnt:EndPnt),2)','b*')
hold on
plot(Epochs_lightAmp_inh, Gstepvar_mean(length(Epochs_lightAmp_exc)+1:end)./mean(Gmean(length(Epochs_lightAmp_exc)+1:end,StartPnt:EndPnt),2)','r*')
title('Gstep var mean/stepmean')
xlabel('abs contrast')
ylabel('mean var of G/ mean G')

% plot individual conductances with peak, upper and lower cross points marked
for a = 1:length(Epochs_str) 
    figure
    plot([1:length(G{a})],G{a})
    hold on
    plot(Gpeakpnt{a}, Gpeak{a}, 'k*')
    for b = 1:length(GlowercrossPnt{a})
    plot(GlowercrossPnt{a}(b), G{a}(b,GlowercrossPnt{a}(b)),'g*')
    plot(GuppercrossPnt{a}(b), G{a}(b,GuppercrossPnt{a}(b)),'r*')
    plot(GreturnPnt{a}(b),G{a}(b,GreturnPnt{a}(b)),'y*')
    end
    clear b
end 
clear a

% subplot(2,2,2)
% plot(Epochs_lightAmp_exc,Gmean_peakpnt(1:length(Epochs_lightAmp_exc)),'b*')
% hold on
% plot(Epochs_lightAmp_inh,Gmean_peakpnt(length(Epochs_lightAmp_exc)+1:end),'r*')
% xlabel('step contrast')
% ylabel('sample point time of mean G peak (nS)')
% legend('exc','inh')
% 
% subplot(2,2,3)
% plot(Epochs_lightAmp_exc,Gmean_lowercrossPnt(1:length(Epochs_lightAmp_exc)),'b*')
% hold on
% plot(Epochs_lightAmp_inh,Gmean_lowercrossPnt(length(Epochs_lightAmp_exc)+1:end),'r*')
% xlabel('step contrast')
% ylabel('lower cross pnt mean G (nS)')
% legend('exc','inh')
% 
% subplot(2,2,4)
% plot(Epochs_lightAmp_exc,Gmean_riseTime(1:length(Epochs_lightAmp_exc)),'b*')
% hold on
% plot(Epochs_lightAmp_inh,Gmean_riseTime(length(Epochs_lightAmp_exc)+1:end),'r*')
% xlabel('step contrast')
% ylabel('rise time of mean G (nS)')
% legend('exc','inh')


% figure 
% for a = 1:length(Epochs_str) 
% if a<=length(Epochs_str_exc) & Epochs_lightAmp(a)<0  ;      % if the epoch set is an excitatory decrement epoch set... 
%     subplot(2,2,1)
%     title('Gvar exc decs')
% elseif a<=length(Epochs_str_exc) & Epochs_lightAmp(a)>0 ;   % exc inc 
%     subplot(2,2,2)
%     title('Gvar exc incs')
% elseif a>length(Epochs_str_exc) & Epochs_lightAmp(a)<0 ;    % inh dec
%     subplot(2,2,3)
%     title('Gvar inh decs')
% else                                                        % inh inc
%     subplot(2,2,4)
%     title('Gvar inh incs')
% end
%     plot(Gvar(a,:))
%     hold on
% end 
% clear a
% 
% figure
% for a = 1:length(Epochs_str)
% subplot(2,2,1)
%     if a<=length(Epochs_str_exc)
%     plot(repmat(Epochs_lightAmp(a),size(Gpeak{a})),Gpeak{a},'.')
%     hold on
%     plot(Epochs_lightAmp(a),Gpeak_mean(a),'+')    
%     else 
%     plot(repmat(Epochs_lightAmp(a),size(Gpeak{a})),Gpeak{a},'.r')
%     hold on
%     plot(Epochs_lightAmp(a),Gpeak_mean(a),'+r')        
%     end
%     xlabel('step contrast')
%     ylabel('peak of G (nS)')
%     legend('exc','inh')
% subplot(2,2,2)
%     if a<=length(Epochs_str_exc)
%     plot(repmat(Epochs_lightAmp(a),size(Gpeakpnt{a})),Gpeakpnt{a},'.')
%     hold on
%     plot(Epochs_lightAmp(a),Gpeakpnt_mean(a),'+')    
%     else 
%     plot(repmat(Epochs_lightAmp(a),size(Gpeakpnt{a})),Gpeakpnt{a},'.r')
%     hold on
%     plot(Epochs_lightAmp(a),Gpeakpnt_mean(a),'+r')        
%     end
%     xlabel('step contrast')
%     ylabel('peakpnt G (nS)')
%     legend('exc','inh')
% subplot(2,2,3)
%     if a<=length(Epochs_str_exc)
%     plot(repmat(Epochs_lightAmp(a),size(GlowercrossPnt{a})),GlowercrossPnt{a},'.')
%     hold on
%     plot(Epochs_lightAmp(a),GlowercrossPnt_mean(a),'+')    
%     else 
%     plot(repmat(Epochs_lightAmp(a),size(GlowercrossPnt{a})),GlowercrossPnt{a},'.r')
%     hold on
%     plot(Epochs_lightAmp(a),GlowercrossPnt_mean(a),'+r')        
%     end
%     xlabel('step contrast')
%     ylabel('lower cross pnt G (nS)')
%     legend('exc','inh')
% subplot(2,2,4)
%     if a<=length(Epochs_str_exc)
%     plot(repmat(Epochs_lightAmp(a),size(GriseTime{a})),GriseTime{a},'.')
%     hold on
%     plot(Epochs_lightAmp(a),GriseTime_mean(a),'+')    
%     else 
%     plot(repmat(Epochs_lightAmp(a),size(GriseTime{a})),GriseTime{a},'.r')
%     hold on
%     plot(Epochs_lightAmp(a),GriseTime_mean(a),'+r')        
%     end   
%     xlabel('step contrast')
%     ylabel('rise time G (nS)')
%     legend('exc','inh')
% end
% clear a






