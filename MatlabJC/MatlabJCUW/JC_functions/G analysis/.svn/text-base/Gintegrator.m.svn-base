% 

CellInfo_str = '031208Bc4_a' ;
epochs_str = {'[42:51]','[52:61]','[62:71]','[72:81]','[82:91]','[92:101]','[102:112]','[113:122]','[135:139]','[140:144]','[146:149]','[152:154]','[158:164]','[168:169,171:174]','[178:184]','[185:186]'} ;
InhEpochsBeyond = 8 ; % beyond this epochnumber there are inh epochs only (this is the last exc epoch)
epochCond_num = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] ;

Color =  {'b','r','g','k','y','c','m','b--','b:','r:','g:','k:','y:','c:','m:','b-.'} ;
fig_num = 1 ;

clear epoch_Idata epochs trace
% Get appropriate cell structure format
load (CellInfo_str) ;

CellInfo = LoadSCIData(CellInfo,1) ;  % for two amps used

EpochCondition = LoadAndSmoothEpochCondition(CellInfo,5000,1) ; 

for a=1:length(epochs_str) % for each set of epochs
    
epochs = str2num(epochs_str{a}) ; % change string to number

for b=1:length(epochs) ; % for each epoch in that set
epoch_Idata{a}(b)= find(EpochCondition(epochCond_num(a)).EpochNumbers == epochs(b)) ; % find the index of those epochs
end

trace{a} = mean(EpochCondition(epochCond_num(a)).EpochData.Data(epoch_Idata{a},:),1) ; % the mean of the epochs

% offset traces
if a>InhEpochsBeyond ; % if its an inhibitory epoch 
    trace{a} = trace{a} - min(trace{a}) ;
    abstrace{a} = trace{a} ;
else
    trace{a} = trace{a} - max(trace{a}) ;
    abstrace{a} = -1*trace{a} ;
end

% take the cumulative sum of the trace during the light flash
CumSumResponse{a} = cumsum(abstrace{a}(5000:10000)) ;

NormCumSumResponse{a} = CumSumResponse{a}/max(CumSumResponse{a}) ;


figure(fig_num)
plot(trace{a},Color{a}) ; 
hold on

figure(fig_num+1)
plot(CumSumResponse{a},Color{a})
hold on

figure(fig_num+2)
plot(NormCumSumResponse{a},Color{a})
hold on

end

title(CellInfo_str)
legend(epochs_str)
ylabel('current (pA)')
xlabel('sample points')




