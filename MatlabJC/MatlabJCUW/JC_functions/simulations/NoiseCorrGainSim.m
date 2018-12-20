ge = [0:100] ;
gi = [0:100] ;

Nc = 10 ;
Ne = 10 ;
Ni = 10 ;

for a=1:length(ge) ;
    for b=1:length(gi) ;
        
        c(a,b) = (ge(a)*gi(b)*Nc)/sqrt((ge(a)^2*gi(b)^2*Nc^2)+(ge(a)^2*Nc*Ni)+(Ne*gi(b)^2*Nc)+(Ne*Ni)) ;
    end
end

plot(c(2,:))
hold on
plot(c(3,:),'r')


