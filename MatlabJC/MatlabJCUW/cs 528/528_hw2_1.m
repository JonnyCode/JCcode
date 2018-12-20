% dayan and Abbott q1 chpter 1

spike(1)=0 ;
for a=2:2000 ;                % for each spike I want to make
spike(a)=spike(a-1)-log(rand)/100 ;  % take the time of the previous spike and...
end

spike2=zeros(1,max(spike)*1000+1) ;%make a zeros matrix to put our spike times in
spike_place=ceil(spike*1000) ;      % put these times in zeros matrix
spike_place2=spike_place(2:end) ;
spike2(1,spike_place2)=1 ;


figure
plot(spike2)

isis=diff(spike);
cv=std(isis)/mean(isis);

%calculating fano

for b=1:100 ;            % for each bin size
    bin_num = floor((length(spike2))/b) ;    % each bin in spike 2 of that size
        for c = 1:bin_num-1 ;               % for each bin 
        num_spikes{b}(c) = sum(spike2(1,c*b:c*b+b));   % number of spikes in each bin 
        end
fano(b)=var(num_spikes{b})/mean(num_spikes{b});
end


figure,
hist(isis*1000,100)

