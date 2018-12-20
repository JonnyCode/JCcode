% greg S. sent this to me to help implement a butterworth and avoid ringing

Wn = F*SampleInterval; %normalized frequency cutoff
[z, p, k] = butter(1,Wn,'high'); %high can be changed to low
[sos,g]=zp2sos(z,p,k); % convesion?
myfilt=dfilt.df2sos(sos,g);
Xfilt = filter(myfilt,X'); %implementation
Xfilt = Xfilt';   