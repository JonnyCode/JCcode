% two cell model testing how SNR impacts pop error for correlated
% population when independence is assumed (or not)

% JC 2018-10-30

nt = 10000 ; % time points

% stim values
s1 = 0 ;
s2 = 1 ; 
stim = randi([s1,s2],1,nt) ; % stim

% cell params
g1 = 1 ; % gain cell 1
g2set = [0:.2:3] ; % gain cell 2

Cvar = .5 ; % std

for loop = 1:length(g2set) ;
    
    g2 = g2set(loop) ; 

    % mean response
    C1mean = stim*g1 ; % high snr
    C2mean = stim*g2 ; % low snr

    % noise
    C1n = normrnd(0,Cvar,1,nt) ; % noise cell 1

    C2n = C1n ; % nc=1
    %C2n = normrnd(0,Cvar,1,nt) ; % nc=0 

    % response
    C1rsp = C1mean + C1n ;
    C2rsp = C2mean + C2n ;

    % bayes decoding
    sigma = [Cvar,0;0,Cvar] ; % assume mc=0
    %sigma = [Cvar,Cvar;Cvar,Cvar] ; % assume nc=1
    prgs1 = mvnpdf([C1rsp',C2rsp'],[s1*g1,s1*g2],sigma) ; % prob rsp given stim 1
    prgs2 = mvnpdf([C1rsp',C2rsp'],[s2*g1,s2*g2],sigma) ; % prob rsp given stim 2
    pr = prgs1+prgs2 ; % probability of the response

    ps1gr = prgs1./pr ;
    ps2gr = prgs2./pr ;

    Sest = ps1gr*s1 + ps2gr*s2 ; % stim est

    % error
    mse(loop) = mean((stim'-Sest).^2) ; 

end

% figure

figure
plot(stim)
hold on
plot(Sest,'r--')

figure % 
plot(g2set,mse)
xlabel('snr cell 2')
ylabel('mse')

