function [Voltage_unique,fV,fV_sem] = fVgClampFinder(voltage,gexc,ginh,Eexc,Einh,Iadd,Capacitance,VbinSize,SI) ;

% this function will find fV (a function of voltage representing intricsic
% properties of the cell), given dynamic clamp data: the voltage, the
% conductances presented (gexc, ginh), their reversal potentials
% (Eexc,Einh), any other constant current injected (pA), the cells capacitance(nF), the voltage bin size (how many volts per bin)
% and the sample interval (sec).  voltage(mV) and gexc,inh(nS) should be in rows for each trial. 
% The units of fV will be mV/sec, see Badel et. al. 2008.

% JC 5/9/11

Iapp = gexc.*(voltage-Eexc) + ginh.*(voltage-Einh) + Iadd ; % current injected into the cell

dVdt = diff(voltage,1,2)/SI ; % voltage derivative

Iion = Iapp - [zeros(size(dVdt,1),1),dVdt]*Capacitance ; % current attributed intrinsic properties

Voltage_ints = round(voltage(:)/VbinSize)*VbinSize ; % values of voltage, binning by VbinSize
Voltage_unique = unique(Voltage_ints) ; % unique voltage values

for a=1:length(Voltage_unique) ; % for every unique voltage integer
    Vbin = find(Voltage_ints==Voltage_unique(a)) ; % find those values that are that unique voltage 
    Idyn(a) = mean(Iion(Vbin)) ; % Iion averaged over VbinSize as a function of voltage
    Idyn_std(a) = std(Iion(Vbin)) ;
    VbinSize(a) = length(Vbin) ;
end

fV = -Idyn/Capacitance ; % intrinsic properties impact on voltage as function of voltage 
fV_sem = Idyn_std./(Capacitance*sqrt(VbinSize)) ; % sem of fV 

