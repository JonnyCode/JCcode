b=1 ;
for  amp_common_noise = 0:.1:1 ; % fraction of noise that is combo of correlated noise
    
a=1 ;
for hp = [-60,-20,20] ; % for each holding potential
    
num_trials = 10 ; % the number of trials    
    
for trial = 1:num_trials ;

% construct gaussian conductance noise
mu = 0 ;                  % mean
sigma = 50 ;             % std
length = 20000 ;         % length of conductance     
filter = 1 ;               % filter (1 negates filter)
G_noise_exc = conv(normrnd(mu,sigma,1,length),filter) ; % convolved w/ filter
G_noise_inh = conv(normrnd(mu,sigma,1,length),filter) ;
G_common_noise = conv(normrnd(mu,sigma,1,length),filter) ;


% add common noise
G_noise2_exc = G_noise_exc*(1-amp_common_noise) + G_common_noise*(amp_common_noise) ; 
G_noise2_inh = G_noise_inh*(1-amp_common_noise) + G_common_noise*(amp_common_noise) ;

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

b = b+1 ;
end

figure, plot([0:.1:1],var_int) ;
xlabel('fraction correlated noise')
ylabel('variance of intermediate current')

figure, plot([0:.1:1],var_exc) ;
xlabel('fraction correlated noise')
ylabel('variance of exc current')

figure, plot([0:.1:1],(var_exc+var_inh)) 
hold on, plot([0:.1:1], var_int,'r')
xlabel('fraction of correlated noise')
ylabel('variance of summed currents and intermediate current')
legend('summed currrents','intermediate currents')

figure, plot([0:.1:1],(var_exc+var_inh)) 
hold on, plot([0:.1:1], var_int*4,'r')
xlabel('fraction of correlated noise')
ylabel('variance of summed currents and 4 x intermediate current')
legend('summed currrents','intermediate currents')

figure, plot([0:.1:1],(var_exc+var_inh)./(var_int*4))
xlabel('fraction of correlated noise')
ylabel('ratio of variances of currents')
axis([0 .9 0 5])

