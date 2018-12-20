% understanding the global orientation selective index

% jc 2/9/2016

xx=[0:359] ;
pref = [45,45+180]

rr=zeros(1,length(xx)) ;
for a=1:length(pref) ;
    rr = rr+ exp(-((xx-pref(a)).^2)/(2*10^2)) ;
end

gosi = sqrt((sum(rr.*sind(2*xx)).^2) +(sum(rr.*cosd(2*xx)).^2))/sum(rr) 
