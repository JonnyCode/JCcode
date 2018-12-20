%function savedata(src,event)
function savedata(tempdata, initialtime, newTimeStamps, newdat)
%global data
size(newdat);
tempdata.tempvoltage = newdat(:,2); %connect the voltage output from the amplifier to channel 1 and so on...
tempdata.tempcurrent = newdat(:,3);
tempdata.tempodorpulse =newdat(:,4);
tempdata.gain=newdat(:,6);
tempdata.mode=newdat(:,7);
% %%define voltage and current gains based on AM systems amplifier
 gain = mean(tempdata.gain);
 mode = mean(tempdata.mode);
 if ((mode-1.8)<0.2)  % voltage clamp
 elseif ((mode-3.8)<0.2)% current clamp
     if ((gain-4.0)<0.1)
		tempdata.voltagegain=10;       % if telegraphed gain is 4.0, then voltage gain in current clamp mode is 1000mV/mV and voltage clips at ~100 mV
     end
     
    if ((gain-5.3)<0.5)				 % if telegraphed gain is 3.5, then voltage gain in current clamp mode is 500mV/mV and voltage clips at ~200 mV
 		tempdata.voltagegain=20;			
     end
     
 	if ((gain-5)<0.1)				 % if telegraphed gain is 3.0, then voltage gain in current clamp mode is 200mV/mV
 		tempdata.voltagegain=50;
     end
        
 	if ((gain-2.5)<0.1)				 % if telegraphed gain is 2.5, then voltage gain in current clamp mode is 100mV/mV
 		tempdata.voltagegain=100;			%  warning: at this gain, bit noise may become a noticeable
    end
 else
     tempdata.voltagegain=1; %%%% COMMENT THIS LINE OUT DURING ACTUAL EXPERIMENT. THIS IS JUST TO assign some value to voltage gain when 
%                         %%%% not using the amplifier
 sprintf('amplifier is off')
 end
 tempdata.voltagegain=1; %%%% COMMENT THIS LINE OUT DURING ACTUAL EXPERIMENT. THIS IS JUST TO assign some value to voltage gain when 
                         %%%% not using the amplifier
 tempdata.currentgain=100;%this is fixed for current clamp.
% % correct gain for voltage and current
 tempdata.tempvoltage = tempdata.tempvoltage*tempdata.voltagegain;
 tempdata.tempcurrent = tempdata.tempcurrent*tempdata.currentgain;
% 
 tempdata.voltage = [tempdata.voltage;tempdata.tempvoltage];
 tempdata.current = [tempdata.current;tempdata.tempcurrent];
tempdata.odorpulse = [tempdata.odorpulse;tempdata.tempodorpulse];
 voltage = tempdata.voltage; 
 current = tempdata.current;
 odorpulse = tempdata.odorpulse;
 save([tempdata.voltagefile ],'voltage');
 save([tempdata.currentfile ],'current');
end
 