function [HHspikes]=HHspikefinder(V,I)

spiketime=zeros(1,length(V))

dV=diff(V)
for a=1:length(V)
    if dV(a)=0 & dV(a-1)>0 
        spiketime(1,a)=1
    end
    if dV(a)>0 & dV(a+1)<0
        spiketime(1,a)=1
    end
end
    