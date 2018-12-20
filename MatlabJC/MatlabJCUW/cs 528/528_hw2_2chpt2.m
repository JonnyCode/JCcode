% question 2 chpter 2 Dayan and Abbott
% find kernel that provides best estimate of linear firing rate

kernel = 45*one_spike_sta/(var(stim));      % see eq 2.6 of D & A

spike_rate=conv(kernel,stim);   %convolve stimulus with kernel

rho_rate=sum(rho)/length(rho);  %find average firing rate of response
rho_rate=rho_rate*500 ;


% find firing rate of rho for comparison

window = ones(1,20) ;
real_firing_rate=conv(window,rho)/20/.002 ;
figure, plot(real_firing_rate)

% normalize r(est)(spike_rate)
norm_spike_rate = spike_rate.*.002 ;

% create spikes based on normalized r(est)
spike_est = poissrnd(norm_spike_rate);
index=find(spike_est>0);
spike_est=zeros(1,length(spike_est));
spike_est(index)=1 ;
figure, plot(spike_est(75:end))

% create r(o) based on the firing rate of r(est)

est_rate=sum(spike_est)/length(spike_est);  %find average firing rate of response
est_rate=est_rate*500 ;



% ploting comparison of estimated spikes and real spikes, first 4000 points 
figure, subplot(4,1,1)
plot(spike_est(100:2100),'r')
subplot(4,1,2)
plot(rho(1:2001))
subplot(4,1,3)
plot(spike_est(2100:4100),'r')
subplot(4,1,4)
plot(rho(2000:4000))
xlabel('time')
ylabel('spike')


% auttocorrelations

auttocorr_est= xcorr(spike_est,'coeff') ;
x_values= 1:length(auttocorr_est) ;
x_values = x_values - (length(x_values)+1)/2;
figure, plot(x_values,auttocorr_est); 

auttocorr_rho=xcorr(rho, 'coeff');
x_values= 1:length(auttocorr_rho) ;
x_values = x_values - (length(x_values)+1)/2;
figure, plot(x_values,auttocorr_rho); 

% isi histograms

i = find(spike_est==1)
isi_synthetic=diff(i)*2  
figure, hist(isi_synthetic,1000)

j = find(rho==1)
isi_rho=diff(j)*2
figure, hist(isi_rho,1000)

% calc Cv

cv_synthetic=std(isi_synthetic)/mean(isi_synthetic)
cv_rho=std(isi_rho)/mean(isi_rho)


% last question 528 (not in Dayan and Abbott)
%pick random waves from the stimulus
rand_point = zeros(length(stm),151);
for a=1:length(stm)         % for each spike preeding wave...
    rand_point = floor(rand(1)*(length(stim)-150));  % choose a random number along the length of the stimulus 
    rand_wave(a,:)=stim(rand_point:rand_point+150,1);  % select the wave 151 points before the random point
    a
end

stm_cov = cov(stm - repmat(one_spike_sta,length(stm),1));  % covariance of residuals (waves precedings spikes minus sta)
stim_cov = cov(rand_wave) ;                             % covariance of random waves from stimulus

[EigVec, EigVal] = eig(stm_cov - stim_cov);        % find eigen values, vectors 

plot(sum(EigVal))
eig_outliers = find(sum(EigVal)<-1850)              % find the two outliers
plot(EigVec(:,1))                                   % ploting top two              
hold on, plot(EigVec(:,2),'r')                      

%project prior waves onto top 2 eigenmodes
project_top_prior = rand_wave*EigVec(:,1);      
project_top2_prior = rand_wave*EigVec(:,2);
figure, plot(project_top_prior,project_top2_prior,'r*');

% project wave preceding spike onto top 2 eigenmodes
project_top = stm*EigVec(:,1);      
project_top2 = stm*EigVec(:,2);
hold on, plot(project_top,project_top2,'*');


% make pdf of projections from the first eigenmode
[priorproj1n,bins] = hist(project_top_prior,50);
priorproj1pdf = priorproj1n/sum(priorproj1n);
[proj1n,bins] = hist(project_top,bins);
proj1pdf = proj1n/sum(proj1n);

% I have several distributions now: P[proj|spike],P[proj],P[spike]what 
% I want is P[spike|proj], which I can get from 
% Bayes rule:P[spike|proj] = P[spike]*P[proj|spike]/P[proj]

figure
bar(bins,priorproj1pdf,'r')
hold on
bar(bins,proj1pdf)

figure
decision1 = proj1pdf./priorproj1pdf*mean(rho)/.002;
bar(bins,decision1)

% make pdf of projections from second eigenmode

[priorproj2n,bins] = hist(project_top2_prior,50);
priorproj2pdf = priorproj1n/sum(priorproj2n);
[proj2n,bins] = hist(project_top2,bins);
proj2pdf = proj2n/sum(proj2n);

% Bays law again

figure
bar(bins,priorproj2pdf,'r')
hold on
bar(bins,proj2pdf)

figure
decision2 = proj2pdf./priorproj2pdf*mean(rho)/.002;
bar(bins,decision2)

%multiply projections from the spike triggered stimuli


for a=1:length(project_top);
    for b=1:length(project_top2);
    projection_product(a,b)=project_top(a)*project_top2(b);
    end
end



% find combo of density difference between prior projection and projection
% smooth trick from T Azevedo

proj=[project_top,project_top2];
cntmat2 = zeros(50,50);

for i = 1:size(proj,1)          % for each set of coordinates...
e = proj(i,:);              

[y,j] = min((bins-e(1)).^2);        % find which bin the first one is closest too
[y,k] = min((bins-e(2)).^2);        % find which bin the second one is closest too

cntmat2(j,k) = cntmat2(j,k)+1;         %add a point to the appropriate bin 

end
bar3(cntmat2)


projprior=[project_top_prior,project_top2_prior];
cntmatprior = zeros(50,50);

for i = 1:size(projprior,1)          % for each set of coordinates...
f = proj(i,:);              

[z,l] = min((bins-f(1)).^2);        % find which bin the first one is closest too
[z,m] = min((bins-f(2)).^2);        % find which bin the second one is closest too

cntmatprior(l,m) = cntmatprior(l,m)+1;         %add a point to the appropriate bin 

end
figure
bar3(cntmatprior)



