% Simulates a general nonlinear IF neuron driven by OU conductance

dt=0.01 ; %ms
Tmax = 500; 

tlist=0:dt:Tmax;
nsteps=length(tlist);
Vlist=zeros(1,nsteps);  %list of voltages that come from simulation
gelist=zeros(1,nsteps);
gilist=zeros(1,nsteps);
Iapplist=zeros(1,nsteps);
tslist=[] ; %list of spike times

%define current-voltage relationship here
C=1;  %capacitance
f = @(V) 1/C*(-V +  V^2) ;  %QIF

%eventually seek to fit:
%C dV/dt = I_ion + I_app
%rewrite
%dV / dt = f(V) + I_app /C
%so we must find both f(V) and C!


%define after-spike current here
tau_AHP=5; %ms
h = @(t) 0.5*1*-exp(-t/tau_AHP)/C ;

%define threshold and reset
Vt=1; Vr=0;
Vlist(1) = Vr;

%define time consts of conductance, mean levels, standard deviations, and reversal potentials
tau_ge=1; Ve=2.5; me=1; stde=.3;  %make sure Ve is over thresh, or get no spikes
tau_gi=2; Vi=0; mi=1; stdi=.3;

%formula for OU:
%dx =  tau(m - x) dt + sigma dW
%gives mean m, time const tau, standard deviation=sigma/sqrt(2 tau)
%so, set sigma = std*sqrt(2 tau)

ts=-Inf;  %set up time of last spike to be back in past, initially

for n=1:(nsteps-1)
    
    Iapplist(n) = gelist(n)*(Ve - Vlist(n)) + gilist(n)*(Vi - Vlist(n)) ;
    
    Vlist(n+1) = Vlist(n) + dt*( f(Vlist(n)) +  h(tlist(n)-ts) + 1/C*(Iapplist(n))) ;

    gelist(n+1) = gelist(n) + dt*tau_ge*(me - gelist(n)) + stde*sqrt(2* tau_ge)*sqrt(dt)*randn ;
    gilist(n+1) = gilist(n) + dt*tau_gi*(me - gilist(n)) + stdi*sqrt(2* tau_gi)*sqrt(dt)*randn ;    

    %spike detect / reset / AHP
    if Vlist(n+1) > Vt ;
        Vlist(n+1) = Vr ;
        ts = tlist(n) ; %time of last spike;
        tslist=[tslist ts];
    end
    
end

figure
subplot(311)
set(gca,'FontSize',16)
plot(tlist,Vlist)
xlabel('t'); ylabel('V')
subplot(312)
set(gca,'FontSize',16)
plot(tlist,gelist)
xlabel('t'); ylabel('ge')
subplot(313)
set(gca,'FontSize',16)
plot(tlist,gilist)
xlabel('t'); ylabel('gi')


%save data if want to use in later model construction
save cell_data tlist gelist gilist Ve Vi Vlist tslist C Iapplist