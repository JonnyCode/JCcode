%function savedata(src,event)
function savedataJC(tempdata, initialtime, newTimeStamps, newdat)

% parameters
HSgain = 10^4 ; %(pA/V) headstage gain with amp set at LOW probe gain 
extFixedVoltageGain = 5 ; % the gain set on the external filter 
extFixedCurrentGain = 10 ; % 

% load data recorded earlier in the same trial
try
    load(tempdata.voltagefile) ;
    load(tempdata.currentfile) ;
catch
    subplot(2,1,1)
    hold off
    subplot(2,1,2)
    hold off
end

% define voltage and current gains based on AM systems amplifier
gain = mean(newdat(:,6));
mode = mean(newdat(:,7));

% gain values - amp output in V but data saved as mV and pA
variableGains = [1,2,5,10,20,50,100] ;

if mode<=2 ; % voltage clamp modes (Vtest, Vcomp, Vclamp)
     tempdata.tempcurrent = newdat(:,2) ;
     tempdata.tempvoltage = newdat(:,4) ;
     tempdata.voltagegain = 1000/(extFixedVoltageGain*10) ; % bessell filter gain at extFixedVoltageGain with fixed ouputs x10Vm
    
     [minV,minI] = min(abs(([1.5:.5:4.5]-gain))) ; 
     tempdata.currentgain = HSgain/variableGains(minI) ;
     
elseif mode>2 ; % current clamp (I=0, Iclamp, Iresist, Ifollow)
     tempdata.tempvoltage = newdat(:,2) ;
     tempdata.tempcurrent = newdat(:,3) ;
     tempdata.currentgain = HSgain/(extFixedCurrentGain*10)  ; % bessell filter gain at extFixedCurrentGain (*10 may be issue on JC rig only)
     
     [minV,minI] = min(abs(([3:.5:6]-gain))) ;
     tempdata.voltagegain = 1000/variableGains(minI) ;     
     
end

% % correct gain for voltage and current
tempdata.tempvoltage = tempdata.tempvoltage*tempdata.voltagegain;
tempdata.tempcurrent = tempdata.tempcurrent*tempdata.currentgain;

if exist('voltage','var') ;
    voltage = [voltage;tempdata.tempvoltage];
    current = [current;tempdata.tempcurrent];
else
    voltage = tempdata.tempvoltage ;
    current = tempdata.tempcurrent ;
end

% analog outputs
Ao0 = newdat(:,1) ;
Ao1 = newdat(:,5) ;

% figure
subplot(5,1,1:2)
plot(newTimeStamps,tempdata.tempvoltage)
hold on
title(['trial = ',num2str(tempdata.n)])
xlabel('time (sec)')
ylabel('mV')

subplot(5,1,3:4)
plot(newTimeStamps,tempdata.tempcurrent)
hold on
xlabel('time (sec)')
ylabel('pA')

subplot(5,1,5)
plot(newTimeStamps,Ao0,'k-')
hold on
plot(newTimeStamps,Ao1,'r-')
xlabel('time (sec)')
ylabel('Ao0 and Ao1 (volts)')


save([tempdata.voltagefile ],'voltage');
save([tempdata.currentfile ],'current');
save([tempdata.Ao0file ],'Ao0');
save([tempdata.Ao1file ],'Ao1');
save([tempdata.timefile ], 'initialtime') ;

end
 

