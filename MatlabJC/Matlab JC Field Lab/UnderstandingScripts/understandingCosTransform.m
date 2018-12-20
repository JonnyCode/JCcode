xx=[0:.0001:1] ;
shift= 5*pi/180 ;
fr=2*pi*1 ; 
yy=cos((fr*xx)+shift) ;
zz=SigFun2([3,.5,.1, 0],yy) ;
[px,ps,p] = PowerPhaseFinder(zz,10000) ;
cc = cos((2*pi*xx)+p(2)/(180/pi)) ;