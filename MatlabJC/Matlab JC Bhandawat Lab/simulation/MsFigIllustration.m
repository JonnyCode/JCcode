% figure illustration

X = 0:800 ;

% bg odor
tau = 0 ;
ofs = 10 ;
Y(1,:) = exp(-X*tau)+ofs ;

% lfp
tau = .01 ;
ofs = 10 ;
Y(2,:) = exp(-X*tau)+ofs ;

% ORN spikes - PN synaptic
tau = .01 ;
ofs = .5 ;
Y(3,:) = exp(-X*tau)+ofs ;

% PN spikes 
tau = .01 ;
ofs = .2 ;
Y(4,:) = exp(-X*tau)+ofs ;

for a=1:4 ;
    Y(a,:) = Y(a,:)/max(Y(a,:)) ;
end

Y = [zeros(4,1000),Y] ; 

plot(Y')

ForIgor.bgRs = Y' ;

% pulses 

X = 0:800 ;

tau = .001 ;
R(1,:) = ones(1,length(X))+2 ;
R(2,:) = ones(1,length(X)) ;
R(3,:) = 1*exp(-X*tau) ; % control
R(4,:) = .7*exp(-X*tau) ; % ORN LFP - ORN spikes
R(5,:) = .4*exp(-X*tau) ; % PN syn - PN spikes

Rplus = [2*ones(1,1000);zeros(4,1000)] ;
R = [Rplus,R,Rplus] ; 
plot(R')

ForIgor.pulseRs = R' ;
