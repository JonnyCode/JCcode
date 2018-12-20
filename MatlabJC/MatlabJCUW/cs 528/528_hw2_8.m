% finding one spike sta for cs 528 hw#2 q8

temp_spike=find(rho==1)'; %find index of each spike
c=find(temp_spike>150);
spike=temp_spike(c);

for a=1:length(spike)   %for each spike
        stm(a,:)=stim(spike(a)-150:spike(a),1) ;  % go into stim and pull out previous 150 points (300ms) 
end 

one_spike_sta=mean(stm)

% finding 2 spike sta for 528 hw#2 q9
for e=1:50;  % for each distance were interested in

dist=[1,zeros(1,e),1];   % number of points between spikes

temp_spike2 = find(conv(rho,dist)==2); % indices of all events with two spikes separated by 2 x dist ms
c2=find(temp_spike2>150);
spike2=temp_spike2(c2);

for d=1:length(spike2)   %for each spike
        stm2(d,:)=stim(spike2(d)-150:spike2(d),1) ; 
end
two_spike_sta{e}=mean(stm2);

% calculating sum of one spike sta
sum_one_spike_sta{e}=[one_spike_sta(1,e+1:end),zeros(1,e)]+one_spike_sta;

% subtract sum from two spike sta
residual_sta{e}=two_spike_sta{e}-sum_one_spike_sta{e};

sum_residual_sta{e}=sum(abs(residual_sta{e}));

clear stm2



% dot product for correlation
for e=1:50
dot_sta(e)=(two_spike_sta{e}*sum_one_spike_sta{e}')/(sqrt(two_spike_sta{e}*two_spike_sta{e}')*sqrt(sum_one_spike_sta{e}*sum_one_spike_sta{e}'));
end



