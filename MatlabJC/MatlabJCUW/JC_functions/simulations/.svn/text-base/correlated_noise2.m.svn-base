% this script was modified from "correlated_noise". 

V = 10 ; % intermediate holding potential

b=1 ;
for  sigma_common = 0:10:100 ; % std of common noise
    
a=1 ;
for hp = [-60,V,20] ; % for each holding potential
    
num_trials = 10 ; % the number of trials    
    
for trial = 1:num_trials ;

% construct gaussian conductance noise
mu = 0 ;                  % mean
sigma_exc = 50 ;             % std
sigma_inh = 50 ;             % std
length = 20000 ;         % length of conductance     
filter = 1 ;               % filter (1 negates filter)
G_noise_exc = conv(normrnd(mu,sigma_exc,1,length),filter) ; % convolved w/ filter
G_noise_inh = conv(normrnd(mu,sigma_inh,1,length),filter) ;
G_common_noise = conv(normrnd(mu,sigma_common,1,length),filter) ;


% add common noise
G_noise2_exc = G_noise_exc + G_common_noise ; 
G_noise2_inh = G_noise_inh + G_common_noise ;

% construct stable conductance
G_stable_exc = sin(.001*[1:length])*800+800 ; % stable exc g
G_stable_inh = sin(.001*[1:length]+(45*pi))*800+800 ; % stable inh g

% add noise to stable conductances
G_exc = G_stable_exc + G_noise2_exc ;
G_inh = G_stable_inh + G_noise2_inh ;

% calculate recorded currents
rp_exc = 20 ; % exc rev potential
rp_inh = -60 ; %inh rev pot
I{a}(trial,1:length) = G_exc*(hp-rp_exc) + G_inh*(hp-rp_inh) ;

end
a = a+1 ;
end

% calculate residuals 
mean_exc = repmat(mean(I{1}),num_trials,1) ; % repmat replicates vector of mean
mean_int = repmat(mean(I{2}),num_trials,1) ;
mean_inh = repmat(mean(I{3}),num_trials,1) ;

resid_exc = I{1} - mean_exc ;
resid_int = I{2} - mean_int ;
resid_inh = I{3} - mean_inh ;

%variance of exc, inh, and exc+inh currents
var_exc(b) = mean(var(resid_exc,0,2)) ;
var_int(b) = mean(var(resid_int,0,2)) ;
var_inh(b) = mean(var(resid_inh,0,2)) ;


% calculate the predicted ratio based on pseudo analytical equation

R(b) = (((sigma_exc^2)*((rp_exc-rp_inh)^2))+((sigma_inh^2)*((rp_exc-rp_inh)^2))+(2*(sigma_common^2)*((rp_exc-rp_inh)^2)))...
    /(((sigma_exc^2)*((rp_exc-V)^2))+((sigma_inh^2)*((V-rp_inh)^2))+abs(((sigma_common*(V-rp_inh))-(sigma_common*(rp_exc-V)))^2)) ;

b = b+1 ;
end

plot([0:10:100],var_int) ;
xlabel('std common noise')
ylabel('variance of intermediate current')

figure, plot([0:10:100],var_exc) ;
xlabel('std common noise')
ylabel('variance of exc current')

figure, plot([0:10:100],var_exc+var_inh)
hold on, plot([0:10:100],var_int,'r')
xlabel('std common noise')
ylabel('variance of currents')
legend('variance of summed current','variance of intermediate current')

figure, plot([0:10:100],var_exc+var_inh)
hold on, plot([0:10:100],4*var_int,'r')
xlabel('std common noise')
ylabel('variance of currents')
legend('variance of summed current','4 x variance of intermediate current')

figure, plot([0:10:100],(var_exc+var_inh)./(var_int))
xlabel('std common noise')
ylabel('ratio of variances of currents')
hold on
plot([0:10:100],R,'r')
legend('simulated', 'analytical')
