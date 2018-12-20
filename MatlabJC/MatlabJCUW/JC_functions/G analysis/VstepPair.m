function ForIgor = VstepPair(Input,Parameters,id,A) ; 

% 12/1/10 this function will get gj conductance between a pair of v-clamped
% cells, one being voltage step while the other records current change

% get data
try
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);
catch
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);
end

epochs = str2num(Input(A).(id)) ;

for a = 1:length(epochs) ; % for each epoch

    [voltageCommand1(a,:), error] = ITCReadEpochStm(epochs(a), 0,fp) ; % get voltage command
    [voltageCommand2(a,:), error] = ITCReadEpochStm(epochs(a), 1,fp) ; 
    
    [data1(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
    [data2(a,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;
    
    [prePnts(a), error] = ITCGetStmPrePts(epochs(a), 0, 0, fp) ;
    [postPnts(a), error] = ITCGetStmTailPts(epochs(a), 0, 0, fp) ;
    
    data1(a,:) = data1(a,:)-mean(data1(a,1:prePnts(a))) ;
    data2(a,:) = data2(a,:)-mean(data2(a,1:prePnts(a))) ;
    
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

gainFactor = 20 ; % from multiclamp settings
voltageCommand1 = voltageCommand1(:,1:length(data1))*gainFactor ;
voltageCommand2 = voltageCommand2(:,1:length(data1))*gainFactor ;

if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end

for a = 1:length(epochs) ; % for each epoch
    Vstep1(a) = round(voltageCommand1(a,prePnts(a)+1)-voltageCommand1(a,1)) ; % if there is a voltage step
    Vstep2(a) = round(voltageCommand2(a,prePnts(a)+1)-voltageCommand2(a,1)) ;
end

UniqueVsteps1 = unique(Vstep1) ;
UniqueVsteps2 = unique(Vstep2) ;

Ievoked1 = cell(1,length(UniqueVsteps2)) ; 
Ievoked2 = cell(1,length(UniqueVsteps1)) ;
Iintrinsic1 = cell(1,length(UniqueVsteps1)) ;
Iintrinsic2 = cell(1,length(UniqueVsteps2)) ; 
Inull1 = [] ;
Inull2 = [] ;

    
for a = 1:length(epochs) ; % for each epoch
    if Vstep2(a)==0 ; % if cell 2 did not have a voltage step 
        Vi = find(UniqueVsteps1==Vstep1(a)) ;
        Ievoked2{Vi}=[Ievoked2{Vi}; data2(a,:)] ; % get current in cell 2
        Iintrinsic1{Vi}=[Iintrinsic1{Vi}; data1(a,:)] ; % get current in cell 1 
    end
    if Vstep1(a)==0 ;
        Vi = find(UniqueVsteps2==Vstep2(a)) ;
        Ievoked1{Vi}=[Ievoked1{Vi}; data1(a,:)] ;
        Iintrinsic2{Vi}=[Iintrinsic2{Vi}; data2(a,:)] ;  
    end
end


for a=1:length(Ievoked1) ;
    if ~isempty(Ievoked1{a}) ;
        Ievoked_mean1(a,:) = mean(Ievoked1{a}) ;
        Iintrinsic_mean2(a,:) = mean(Iintrinsic2{a}) ;
    end
    
end
        
for a=1:length(Ievoked2) ;
    if ~isempty(Ievoked2{a}) ;
        Ievoked_mean2(a,:) = mean(Ievoked2{a}) ;
        Iintrinsic_mean1(a,:) = mean(Iintrinsic1{a}) ;
    end
end

Ievoked_lastPnt1 = Ievoked_mean1(:,end-postPnts(1)) ;
Ievoked_lastPnt2 = Ievoked_mean2(:,end-postPnts(1)) ;

coefs1 = polyfit(UniqueVsteps1',Ievoked_lastPnt2,1) ;
coefs2 = polyfit(UniqueVsteps2',Ievoked_lastPnt1,1) ;

g1to2 = coefs1(1) ;
g2to1 = coefs2(1) ;

figure
subplot(2,2,1)
plot([1:length(Ievoked_mean2)],Ievoked_mean2)

subplot(2,2,2)
plot(UniqueVsteps1,Ievoked_lastPnt2,'*')
hold on
plot(UniqueVsteps1,UniqueVsteps1*coefs1(1)+coefs1(2))

subplot(2,2,3)
plot([1:length(Ievoked_mean1)],Ievoked_mean1)

subplot(2,2,4)
plot(UniqueVsteps2,Ievoked_lastPnt1,'r*')
hold on
plot(UniqueVsteps2,UniqueVsteps2*coefs2(1)+coefs2(2),'r')    
    
    
    
    
    
    
    
    
    
    
    
    