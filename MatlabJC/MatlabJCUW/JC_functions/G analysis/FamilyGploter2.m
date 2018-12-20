function [Imeans, ImeansOffset, ImeansOffsetMinusInh, ImeanOffsetMinusInh_sumExc, ImeanOffsetMinusInh_peakExc, Gexc_sum, Gexc_peak] = FamilyGploter2(CellInfo_str,epochs_str,epochCond_num,epochHold) ;

% This function takes a set of epochs recorded at various holding
% potentials, offsets them according to their standing exitation and
% inhibition (using a nonlinearity factor to offset the exc), and pulls out
% conductance - voltage ralationships
% JC 1/3/08 (modified from previous version)

% % WHEN RUN AS A SCRIPT
% CellInfo_str = '102307Ac3_a' ;
% epochs_str = {'[109:113]','[114:118]','[119:123]','[124:128]','[129:133]','[134:138]','[139:143]'} ;
% epochCond_num = [1,1,1,1,1,1,1,1] ;
% epochHold = [-80,-60,-40,-20,0,20,40] ; % holding potential at which the epoch set was collected 

% PARAMETERS
Eexc = 20; % exc reversal potential
Einh = -60; % inh reversal potential

excstart = 5000 ; % sample points where exc current is best isolated
excend = 7000 ; 

% FROM PARAMETERS
ExcIndex = find(epochHold == Einh) ; 
InhIndex = find(epochHold == Eexc) ; % index of epochHold which isolated inhibition  
EpochsIndexReversed = find(epochHold>Eexc) ; %find the index of the epochsets that reverse the exc current to an outward current

% Get appropriate cell structure format
cd ~/data_analysis/Index ;
load (CellInfo_str) ;

CellInfo = LoadSCIData(CellInfo,1) ;  % for two amps used

EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ; 

% get mean data from listed epochs
for a=1:length(epochs_str) % for each set of epochs...
clear epoch_Idata % clear previous epoch data 
    
epochs = str2num(epochs_str{a}) ; % change string to number

for b=1:length(epochs) ; % for each epoch in that set...
epoch_Idata(b)= find(EpochCondition(epochCond_num(a)).EpochNumbers == epochs(b)) ; % find the index of those epochs
end

Imeans(a,:) = mean(EpochCondition(epochCond_num(a)).EpochData.Data(epoch_Idata,:)) ; 
end
clear a b

% OFFSET CURRENTS BY THIER STANDING CONDUCTANCES

% calculate the Ginh and subtract it from the mean current to isolate Iexc
GinhMean_preoffset = Imeans(InhIndex,:)/(Eexc - Einh) ; % calculate the conductance of the mean Inhibitory current

for a = 1:size(Imeans,1) ; %for each recorded current mean...
ImeansMinusInh_preoffset(a,:) = Imeans(a,:) - (GinhMean_preoffset*(epochHold(a)-Einh)) ;   % subtract the inhibtory current (scaled by driving force)
ImeansMinusInh_peakExc_preoffset(a) = min(ImeansMinusInh_preoffset(a,excstart:excend)) ;      %find the peak exc current between excstart and excend 
end
clear a

for a = 1:length(EpochsIndexReversed) ; % for each epoch that the current is reversed..
    ImeansMinusInh_peakExc_preoffset(EpochsIndexReversed(a)) = max(ImeansMinusInh_preoffset(EpochsIndexReversed(a),excstart:excend)) ; % we want the max not the min of the peak
end
clear a

% find offset for standing excitation using the nonlinearity found in "ImeansMinusInh_peakExc" vs epochHold (IV curve)
HighPntExc = max(Imeans(ExcIndex,:)) ; % find the highest point on the exc trace 

Gpeaks_preoffset = ImeansMinusInh_peakExc_preoffset./(epochHold - Eexc) ; % find the coductance for the peak in order to describe the conductance nonlineraity 

OffsetFactors = (Gpeaks_preoffset./Gpeaks_preoffset(ExcIndex)) ; % this factor accounts for the NMDA nonlinearity

Ioffset_exc = (OffsetFactors*(HighPntExc/(Einh-Eexc))).*(epochHold - Eexc) ; % these offsets account for both expected linearity (driving force) and the NMDA nonlinearity  
Ioffset_exc(InhIndex) = 0 ; % this is necessary because the coductance of the nonlinearity at this point is unkowable and results in (-Inf) but this doesn't matter bc the driving force = 0

% find offset of standing inh (assuming linearity)
LowPntInh = min(Imeans(InhIndex,:)) ; % find the lowest point in the inh trace

Ioffset_Inh = (LowPntInh/(Eexc-Einh))*(epochHold - Einh) ; % these points assume linearity 

% offset currents by standing excitiation and inhibition
for a = 1:size(Imeans,1) ; % for each recorded current mean ...
ImeansOffset(a,:) = Imeans(a,:) - Ioffset_Inh(a) - Ioffset_exc(a) ; 
end


% NOW THAT THE CURRENTS ARE OFFSET BY THEIR ASSUMED STANDING CONDUCTANCES ...

% calculate the Ginh and subtract it from the mean current to isolate Iexc
Ginh_mean = ImeansOffset(InhIndex,:)/(Eexc-Einh) ; % calculate the conductance of the mean Inhibitory current

for a = 1:size(Imeans,1) ; %for each recorded current mean...
ImeansOffsetMinusInh(a,:) = ImeansOffset(a,:) - Ginh_mean*(epochHold(a)-Einh) ;   % subtract the inhibtory current (scaled by driving force)
ImeanOffsetMinusInh_peakExc(a) = min(ImeansOffsetMinusInh(a,excstart:excend)) ;      %find the peak exc current between excstart and excend 
end
clear a

for a = 1:length(EpochsIndexReversed) ; % for each epoch that the current is reversed..
    ImeanOffsetMinusInh_peakExc(EpochsIndexReversed(a)) = max(ImeansOffsetMinusInh(EpochsIndexReversed(a),excstart:excend)) ; % we want the max not the min of the peak
end
clear a

ImeanOffsetMinusInh_sumExc = sum(ImeansOffsetMinusInh(:,excstart:excend),2)' ; % sum of the Iexc

% change IV relations into GV relations
Gexc_sum = ImeanOffsetMinusInh_sumExc./(epochHold-Eexc) ;
Gexc_sum(InhIndex) = NaN ; % cannot assess G at Eexc

Gexc_peak = ImeanOffsetMinusInh_peakExc./(epochHold-Eexc) ;
Gexc_peak(InhIndex) = NaN ; % cannot assess G at Eexc


% FIGURES TO PLOT 
% the mean currents
figure
subplot(2,1,1)
plot([1:length(Imeans)],Imeans)
title({CellInfo_str, 'mean currents'})
legend(epochs_str)
ylabel('current (pA)')
xlabel('sample points')

% the mean current offset by their stading g  
subplot(2,1,2)
plot([1:length(Imeans)],ImeansOffset)
title({CellInfo_str, 'mean currents offset'})
legend(epochs_str)
ylabel('current (pA)')
xlabel('sample points')

% the mean currents with Iinh subtracted (Iexc isolated)  
figure
subplot(2,2,1:2)
plot([1:length(Imeans)],ImeansOffsetMinusInh)
title({CellInfo_str, 'mean currents offset minus Inhibition'})
legend(epochs_str)
ylabel('current (pA)')
xlabel('sample points')

% the sum of the assumed exc current during the period evaluated of the recroding 
subplot(2,2,3)
plot(epochHold, ImeanOffsetMinusInh_sumExc) 
title({CellInfo_str, 'sum'})
ylabel('current Integrals (pA)')
xlabel('holding potential of each epoch set')

% the peak of the currents
subplot(2,2,4)
plot(epochHold,ImeanOffsetMinusInh_peakExc,'r') % the peak exc between the time period of exc analyzed
title({CellInfo_str, 'peak mean currents offset minus Inhibition vs holding potential'})
ylabel('current Integrals (pA)')
xlabel('holding potential of each epoch set')
text('units','normalized') %this sets the text in units where (1,1) is upper right
text(.1, .9,['exc start = ' num2str(excstart)],'units','normalized')
text(.1, .85,['exc end = ' num2str(excend)],'units','normalized')

% the sum of the assumed exc conductance during the period evaluated of the recroding   
figure
subplot(2,2,1)
plot(epochHold,Gexc_sum,'b*') ;
title({CellInfo_str,'sum Gexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('sum Gexc during evaluation time')

% the peak Gexc 
subplot(2,2,2)
plot(epochHold,Gexc_peak,'b*')
title({CellInfo_str,'peak Gexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('peak Gexc during evaluation time')
hold on
plot(epochHold,Gpeaks_preoffset,'r*')
legend('postoffset','preoffset')

end % end function

















