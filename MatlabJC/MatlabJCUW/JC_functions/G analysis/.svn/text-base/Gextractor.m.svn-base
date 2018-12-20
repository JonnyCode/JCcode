Eexc = 20; % exc reversal potential

for R = 1:2

if R ==1    
% control
CellInfo_str = '102307Ac3_a' ;
epochs_str = {'[57:61]','[52:56]','[62:66]','[67:71]','[72:76]','[77:83]','[84:88]'} ;
epochCond_num = [1,1,1,1,1,1,1,1] ;
epochHold = [-80,-60,-40,-20,0,20,40] ; % holding potential at which the epoch set was collected 
end

if R ==2
% drug
CellInfo_str = '102307Ac3_a' ;
epochs_str = {'[109:113]','[114:118]','[119:123]','[124:128]','[129:133]','[134:138]','[139:143]'} ;;
epochCond_num = [1,1,1,1,1,1,1,1] ;
epochHold = [-80,-60,-40,-20,0,20,40] ; % holding potential at which the epoch set was collected 
end

if R ==3
% wash
CellInfo_str = '110707Ac1_a' ;
epochs_str = {'[220:224]','[225:229]','[230:234]','[235:239]','[240:244]','[245:249]','[250:264]'} ;
epochCond_num = [1,1,1,1,1,1,1,1] ;
epochHold = [-80,-60,-40,-20,0,20,40] ; % holding potential at which the epoch set was collected 
end

[Imeans{R}, ImeansOffset{R}, ImeansOffsetMinusInh{R}, ImeanOffsetMinusInh_sumExc{R}, ImeanOffsetMinusInh_peakExc{R}, Gexc_sum{R}, Gexc_peak{R}]...
    = FamilyGploter2(CellInfo_str,epochs_str,epochCond_num,epochHold) ;

end % end R loop

% calculate difference in current between control and drug treatment
sumIExc_diff = ImeanOffsetMinusInh_sumExc{1} - ImeanOffsetMinusInh_sumExc{2} ;
peakIExc_diff = ImeanOffsetMinusInh_peakExc{1} - ImeanOffsetMinusInh_peakExc{2} ;

% make the above from IV into GV relations
sumGexc_diff = sumIExc_diff./(epochHold - Eexc) ;
peakGexc_diff = peakIExc_diff./(epochHold - Eexc) ;


% FIGURES
figure
% sum Iexc
subplot(2,2,1)
plot(epochHold,ImeanOffsetMinusInh_sumExc{1})
hold on
plot(epochHold,ImeanOffsetMinusInh_sumExc{2},'r')
plot(epochHold,ImeanOffsetMinusInh_sumExc{3},'g')
legend ('control','drug','wash')
title({CellInfo_str, 'sum'})
ylabel('current Integrals (pA)')
xlabel('holding potential of each epoch set')

% peak Iexc
subplot(2,2,2)
plot(epochHold,ImeanOffsetMinusInh_peakExc{1}) % the peak exc between the time period of exc analyzed
hold on
plot(epochHold,ImeanOffsetMinusInh_peakExc{2},'r')
plot(epochHold,ImeanOffsetMinusInh_peakExc{3},'g')
legend('control', 'drug','wash')
title({CellInfo_str, 'peak mean currents offset minus Inhibition vs holding potential'})
ylabel('current Integrals (pA)')
xlabel('holding potential of each epoch set')

% the sum Gexc
subplot(2,2,3)
plot(epochHold,Gexc_sum{1},'b*') ;
hold on
plot(epochHold,Gexc_sum{2},'r*') ;
plot(epochHold,Gexc_sum{3},'g*') ;
legend('control','drug','wash')
title({CellInfo_str,'sum Gexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('sum Gexc during evaluation time')

% the peak Gexc 
subplot(2,2,4)
plot(epochHold,Gexc_peak{1},'b*')
hold on
plot(epochHold,Gexc_peak{2},'r*')
plot(epochHold,Gexc_peak{3},'g*')
legend('control','drug','wash')
title({CellInfo_str,'peak Gexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('peak Gexc during evaluation time')

figure
% I difference between drug and control
subplot(2,2,1)
plot(epochHold,sumIExc_diff)
title({CellInfo_str,'diff sumIexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('Icontrol-Idrug')

subplot(2,2,2)
plot(epochHold,peakIExc_diff)
title({CellInfo_str,'diff peakIexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('Icontrol-Idrug')

subplot(2,2,3)
plot(epochHold,sumGexc_diff,'b*')
title({CellInfo_str,'diff sumGexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('Icontrol-Idrug/df')

subplot(2,2,4)
plot(epochHold,peakGexc_diff,'b*')
title({CellInfo_str,'diff peakGexc vs holding pot'})
xlabel('holding potential of each epoch set')
ylabel('Icontrol-Idrug/df')

