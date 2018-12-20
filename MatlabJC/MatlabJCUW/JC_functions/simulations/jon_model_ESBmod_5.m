% Can we get correlated signals and anticorrelated noise in simple model
% for neural integration?  Consider split field experiment - half of inputs
% see +stm, other half -stm.  Exc and Inh are rectified.  Noise in inputs
% is spatially uncorrelated.  

% FR and JC 
% JC modified 1/11/11

clear all

NumTrials = 5000;  %500
NumPts = 1000; %1000
Thresh = -.2;          % rectification

% FIX TIME-DEP SIGNAL PRESENTED OVER MANY TRIALS
% signal - will be + for group1 and - for group2
% time along rows (horiz), trials down cols (vert)
signal=normrnd(0, 1, 1, NumPts);
%this is matrix with same sig repeated over multiple rows (trials)
sig = ones(NumTrials,1)*signal;

full_field_stim=0;

if full_field_stim
    sig1=sig;
    sig2=sig;
else %default to jon's original half-field
    sig1=sig;
    sig2=-sig;
end


%weights 
w_e_pre=1;
w_i_pre=1;
w_e_post=1;
w_i_post=1;


clear response_exc response_inh resid_exc resid_inh

    noise1 = normrnd(0, 1, NumTrials , NumPts);      % group 1 noise
    noise2 = normrnd(0, 1, NumTrials , NumPts);      % group 2 noise

    % create group 1 exc and inh
    exc1 = w_e_pre*(sig1 + noise1);
    exc1(find(exc1 < Thresh)) = 0;
    inh1 = w_i_pre*(-sig1 - noise1);
    inh1(find(inh1 < Thresh)) = 0;

    % create group 2 exc and inh
    exc2 = w_e_pre*(sig2 + noise2);
    exc2(find(exc2 < Thresh)) = 0;
    inh2 = w_i_pre*(-sig2 - noise2);
    inh2(find(inh2 < Thresh)) = 0;

    % store integrated responses
    response_exc = exc1 + exc2;  
    response_inh = inh1 + inh2;
    
    %shuffle across trials.  After random shuffling of rows, each column 
    %will still represent a single stimulus value, but noises will be
    %shuffled across trials, so that have no corr between exc and inh paths

    %append a first random col
    response_exc_temp=[rand(size(response_exc,1),1) response_exc];
    response_exc_shuffled = sortrows(response_exc_temp,1);
    %ditch that first col
    response_exc_shuffled=response_exc_shuffled(:,2:end);
    

    %append a first random col
    response_inh_temp=[rand(size(response_inh,1),1) response_inh];
    response_inh_shuffled = sortrows(response_inh_temp,1);
    %ditch that first col
    response_inh_shuffled=response_inh_shuffled(:,2:end);

%need covariance b/w e and i responses for given signal values
for i=1:NumPts
    covar_temp=cov(response_exc(:,i),response_inh(:,i)) ;
    var_exc(i)=covar_temp(1,1);
    var_inh(i)=covar_temp(2,2);
    covar(i)=covar_temp(1,2);   
end
    

response = w_e_post*response_exc + w_i_post*response_inh; % WHY ADD?
mean_response=mean(response);
mean_exc=mean(response_exc);
mean_inh=mean(response_inh);

var_response=var(response);
std_response=sqrt(var_response);

response_shuffled = w_e_post*response_exc_shuffled - w_i_post*response_inh_shuffled;
var_response_shuffled=var(response_shuffled);
std_response_shuffled=sqrt(var_response_shuffled);


% figure
% set(gca,'FontSize',12)
% errorbar(signal,mean_response,std_response,'.')
% xlabel('signal value')
% ylabel('response and std dev')


%compute linear fisher information that resp carries about stim.  this is
%the inverse of the variance of the optimal unbiased linear readout of
%signal sig, and equals  [d mean(sig)/d sig]^2/ var(sig):  square of deriv of mean
%response w.r.t signal sig divide by variance of responses that occur at
%that signal.

%sort responses in order of signal values
data_matrix_sorted=sortrows([signal' mean_response' var_response' var_response_shuffled' var_exc' var_inh' covar' mean_exc' mean_inh'],1); 
signal_sorted = data_matrix_sorted(:,1)';
mean_response_sorted = data_matrix_sorted(:,2)';
var_response_sorted = data_matrix_sorted(:,3)';
var_response_shuffled_sorted = data_matrix_sorted(:,4)';
var_exc_sorted=data_matrix_sorted(:,5)';
var_inh_sorted=data_matrix_sorted(:,6)';
covar_sorted=data_matrix_sorted(:,7)';
mean_exc_sorted=data_matrix_sorted(:,8)';
mean_inh_sorted=data_matrix_sorted(:,9)';


%smooth variance and mean (in particular so can sensibly take derivative)
poly_order=6;
pmean=polyfit(signal_sorted,mean_response_sorted,poly_order);
mean_response_sorted_smoothed=polyval(pmean,signal_sorted);
pvar=polyfit(signal_sorted,var_response_sorted,poly_order);
var_response_sorted_smoothed=polyval(pvar,signal_sorted);
pvar_shuffled=polyfit(signal_sorted,var_response_shuffled_sorted,poly_order);
var_response_shuffled_sorted_smoothed=polyval(pvar_shuffled,signal_sorted);
pexc=polyfit(signal_sorted,var_exc_sorted,poly_order);
var_exc_sorted_smoothed=polyval(pexc,signal_sorted);
pinh=polyfit(signal_sorted,var_inh_sorted,poly_order);
var_inh_sorted_smoothed=polyval(pinh,signal_sorted);
pcovar=polyfit(signal_sorted,covar_sorted,poly_order);
covar_sorted_smoothed=polyval(pcovar,signal_sorted);
pexc=polyfit(signal_sorted,mean_exc_sorted,poly_order);
mean_exc_sorted_smoothed=polyval(pexc,signal_sorted);
pinh=polyfit(signal_sorted,mean_inh_sorted,poly_order);
mean_inh_sorted_smoothed=polyval(pinh,signal_sorted);

mean_respose_deriv = diff(mean_response_sorted_smoothed)./diff(signal_sorted); %fwd diff
mean_exc_deriv = diff(mean_exc_sorted_smoothed)./diff(signal_sorted); %fwd diff
mean_inh_deriv = diff(mean_inh_sorted_smoothed)./diff(signal_sorted); %fwd diff

lin_fisher_info_response = mean_respose_deriv.^2 ./ var_response_sorted_smoothed(1:end-1);

for k=1:(length(signal_sorted)-1)
    fprime=[mean_exc_deriv(k) ; mean_inh_deriv(k)] ;
    C=[var_exc_sorted_smoothed(k) covar_sorted_smoothed(k) ; covar_sorted_smoothed(k) var_inh_sorted_smoothed(k)];
    lin_fisher_info_response_full(k) = fprime'*inv(C)*fprime;
    Cshuffled=[var_exc_sorted_smoothed(k) 0 ; 0 var_inh_sorted_smoothed(k)];
    lin_fisher_info_response_shuffled_full(k) = fprime'*inv(Cshuffled)*fprime;
end


avg_lin_fisher_info_response=sum(lin_fisher_info_response)/length(lin_fisher_info_response)
avg_lin_fisher_info_response_full=sum(lin_fisher_info_response_full)/length(lin_fisher_info_response_full)
avg_lin_fisher_info_response_shuffled_full=sum(lin_fisher_info_response_shuffled_full)/length(lin_fisher_info_response_shuffled_full)

lin_fisher_info_shuffled_response = mean_respose_deriv.^2 ./ var_response_shuffled_sorted_smoothed(1:end-1);
avg_lin_fisher_info_shuffled_response=sum(lin_fisher_info_shuffled_response)/length(lin_fisher_info_shuffled_response);


figure
subplot(121)
set(gca,'FontSize',12)
errorbar(signal_sorted,mean_response_sorted,sqrt(var_response_sorted))
xlabel('signal value')
ylabel('response and std dev')
subplot(122)
set(gca,'FontSize',12)
errorbar(signal_sorted,mean_response_sorted,sqrt(var_response_shuffled_sorted))
xlabel('signal value')
ylabel('response and std dev')
title('shuffled')

figure
subplot(311)
plot(signal_sorted,mean_inh_sorted_smoothed); hold on
plot(signal_sorted,mean_inh_sorted,'.'); hold on
plot(signal_sorted,mean_exc_sorted_smoothed); hold on
plot(signal_sorted,mean_exc_sorted,'.'); hold on
title('mean exc and inh')
xlabel('signal')
subplot(312)
plot(signal_sorted,var_exc_sorted_smoothed); hold on
plot(signal_sorted,var_exc_sorted,'.'); hold on
plot(signal_sorted,var_inh_sorted_smoothed); hold on
plot(signal_sorted,var_inh_sorted,'.'); hold on
title('var exc and inh')
xlabel('signal')
subplot(313)
plot(signal_sorted,covar_sorted_smoothed); hold on
plot(signal_sorted,covar_sorted,'.'); hold on
title('covar')
xlabel('signal')

figure
set(gca,'FontSize',16)
%drop firstand last ndrop points
ndrop=round(length(signal)/50);
plot(signal_sorted(ndrop:end-ndrop-1),lin_fisher_info_response_full(ndrop:end-ndrop),'r-','LineWidth',3)
xlabel('signal value')
ylabel('linear fisher inf')
hold on
plot(signal_sorted(ndrop:end-ndrop-1),lin_fisher_info_response_shuffled_full(ndrop:end-ndrop),'k--','LineWidth',3)
xlabel('signal value')
ylabel('linear fisher inf')
legend('full','shuffled')
title(sprintf('fisher=%g, fisher shuffled=%g',avg_lin_fisher_info_response_full,avg_lin_fisher_info_response_shuffled_full))
maxval1=max([ lin_fisher_info_response_full(ndrop:end-ndrop) lin_fisher_info_response_shuffled_full(ndrop:end-ndrop)])*1.05;
ylim([0 maxval1])
xlim([-Inf Inf])
%----------------------------

% how correlated are signals?  
% signals are averages across many trials for a given sig.
% Average down cols (in direction 1)
cc_temp=corrcoef(mean(response_exc), mean(response_inh)) ;
cc_signal = cc_temp(1,2) 


resid_exc = response_exc - ones(NumTrials,1)*mean(response_exc);
resid_inh = response_inh - ones(NumTrials,1)*mean(response_inh);

% how correlated are residuals?  Note, this is avg across many signals AND
% trials
cc_temp = corrcoef(resid_exc, resid_inh) ;
cc_noise = cc_temp(1,2)

figure
subplot(121)
plot((mean(response_exc)),(mean(response_inh)),'k.')
xlabel('exc mean'),ylabel('inh mean')
axis square
subplot(122)
plot(resid_exc,resid_inh,'k.')
xlabel('exc resid'),ylabel('inh resid')
axis square
