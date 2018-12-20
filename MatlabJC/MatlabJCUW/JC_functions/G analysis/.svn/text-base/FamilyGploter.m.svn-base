function [Imeans, ImeansOffset, ImeansOffsetMinusInh, ImeanOffsetMinusInh_sumExc, ImeanOffsetMinusInh_peakExc, Gexc_sum, Gexc_peak] = FamilyGploter(CellInfo_str,epochs_str,epochCond_num,epochHold) ;

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

% Get appropriate cell structure format
cd ~/data_analysis/Index ;
load (CellInfo_str) ;

CellInfo = LoadSCIData(CellInfo,1) ;  % for two amps used

EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ; 

for a=1:length(epochs_str) % for each set of epochs...
clear epoch_Idata % clear previous epoch data 
    
epochs = str2num(epochs_str{a}) ; % change string to number

for b=1:length(epochs) ; % for each epoch in that set...
epoch_Idata(b)= find(EpochCondition(epochCond_num(a)).EpochNumbers == epochs(b)) ; % find the index of those epochs
end

Imeans(a,:) = mean(EpochCondition(epochCond_num(a)).EpochData.Data(epoch_Idata,:)) ; 
end
clear a b

% find lowest point on trace which isolated inhibition
InhIndex = find(epochHold == Eexc) ; % index of epochHold which isolated inhibition  
LowPntInh = min(Imeans(InhIndex,:)) ; % find the lowest point on the trace
Goffset_Inh = LowPntInh/(Eexc - Einh) ; % the conductance offset will = the current / the driving force for that current

% use this GoffsetInh value to offset all traces appropriately for standing Inh 
for c = 1:size(Imeans,1) ; %for each recorded current mean...
    Ioffset_Inh(c) = Goffset_Inh*(epochHold(c) - Einh) ; % calculate the appropriate current offset for each trace 
    ImeansOffset_Inh(c,:) = Imeans(c,:) - Ioffset_Inh(c) ; % offset each trace appropriatley
end
clear c

% find the highest point on the trace which isolated excitation
ExcIndex = find(epochHold ==Einh) ; % index of epochHold which isolated excitation
HighPntExc = max(Imeans(ExcIndex,:)) ; % find the highest point on the trace
Goffset_Exc = HighPntExc/(Einh - Eexc) ; % the conductance offset (as above)

% use this GoffsetExc value to offset all trace appropriately for standing
% exc, so that after this loop everthing above zero is inh and all below
% exc
for c = 1:size(Imeans,1) ; %for each recorded current mean...
    Ioffset_Exc(c) = Goffset_Exc*(epochHold(c) - Eexc) ; % calculate the appropriate current offset for each trace 
    ImeansOffset(c,:) = ImeansOffset_Inh(c,:) - Ioffset_Exc(c) ; % offset each trace appropriatley
end
clear c

% calculate the Ginh and subtract it from the mean current to isolate Iexc
Ginh_mean = ImeansOffset(InhIndex,:)/(Einh - Eexc) ; % calculate the conductance of the mean Inhibitory current

for c = 1:size(Imeans,1) ; %for each recorded current mean...
ImeansOffsetMinusInh(c,:) = ImeansOffset(c,:) - Ginh_mean*(Einh - epochHold(c)) ;   % subtract the inhibtory current (scaled by driving force)
ImeanOffsetMinusInh_peakExc(c) = min(ImeansOffsetMinusInh(c,excstart:excend)) ;      %find the peak exc current between excstart and excend 
end
clear c

EpochsIndexReversed = find(epochHold>Eexc) ; %find the index of the epochsets that reverse the current to an outward current
for c = 1:length(EpochsIndexReversed) ; % for each epoch that the current is reversed..
    ImeanOffsetMinusInh_peakExc(EpochsIndexReversed(c)) = max(ImeansOffsetMinusInh(EpochsIndexReversed(c),excstart:excend)) ; % we want the max not the min of the peak
end
clear c

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

% the mean current offset by their stading inh  
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
plot(epochHold,Gexc_peak,'r*')
title({CellInfo_str,'peak Gexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('peak Gexc during evaluation time')

end % end function