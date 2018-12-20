% This script runs a function precision model with different conductance
% parameters.
clear all
count=1                             % set counter

% set conductance parameters
for f_inh = .03:.02:.05 ;             % frequency of inhibitory conductance 
for f_exc = .03:.02:.05 ;             % frequency of excitatory conductance

for amp_inh = 0:5:10 ;        % amplitude of inhibitory synaptic coductance
for amp_exc = 5:5:10 ;        % amplitude of excitatory synaptic conductance

for phase = pi/4:pi/4:(3*pi)/4 ;          % relative phase of sinusiods of conductances
   
% run function    
stats(:,count) = precision_model(f_inh,f_exc,amp_inh,amp_exc,phase) ; 

count= count+1
end
end
end
end
end

% stats=[f_inh;f_exc;amp_inh;amp_exc;phase;average_spkd_norm;stdev_spkd_norm];
% plot one variable data (ie. paramter vs spk d norm)
figure, plot(stats(6,:),stats(1,:),'*')
title('frequency inhibitory')

figure, plot(stats(6,:),stats(2,:),'*')
title('frequency excitatory')

figure, plot(stats(6,:),stats(3,:),'*')
title('amp inhibitory')

figure, plot(stats(6,:),stats(4,:),'*')
title('amp excitatory')

figure, plot(stats(6,:),stats(5,:),'*')
title('phase')

% more plotting

figure, plot(stats(6,:),(stats(4,:)./stats(3,:)),'*')
title('amp exc/amp inh')

figure, plot3(stats(4,:),stats(3,:),stats(6,:),'*')

figure, plot3(stats(5,:),(stats(4,:)./stats(3,:)),stats(6,:),'*')




