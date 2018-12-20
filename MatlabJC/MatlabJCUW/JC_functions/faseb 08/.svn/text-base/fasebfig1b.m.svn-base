% figure 1b for faseb poster '08

% JC 6/26/08

% data = {Cell,condition1 epochs, condition2 epochs, etc}
conditions = {'control dec','apb dec','wash dec','control inc','apb inc','wash inc'}

data{1} = {'052908Bc1','[433:8:505]','[603:2:621]','[843:2:861]','[432:8:505]','[602:2:621]','[842:2:861]'} ;
data{2} = {'052908Bc3','[629:8:733]','[963:2:981]','[1103:2:1121]','[628:8:733]','[962:2:981]','[1102:2:1121]'} ;
data{3} = {'052208Bc3','[1158:8:1230]','[1252:8:1310]','[]','[1157:8:1230]','[1251:8:1320]','[]'} ;
data{4} = {'050808Bc1','[594:8:666]','[758:8:830]','[]','[593:8:666]','[757:8:830]','[]'} ;
data{5} = {'060908Bc4','[158:162]','[242:244]','[]','[187:191]','[207:211]','[]'} ;
data{6} = {'060908Bc1','[243:247]','[303:307]','[]','[203:207]','[308:312]','[]'} ;
data{7} = {'061008Bc4','[172:176]','[202:206]','[]','[167:171]','[207:211]','[]'} ;
data{8} = {'061008Bc3','[198:202]','[265:269]','[]','[230:234]','[260:264]','[]'} ;

Color = {'b','r','g','k','y','c','m','b-.','r-.','g-.','k-.','y-.','c-.','m-.','b:','r:','g:','k:'} ;

for cell = 1:length(data) ;             % for each cell
    [fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',data{cell}{1}]) ; % get the file
    
    for condition = 1:length(conditions) ;              % for each condition in this cell
        epochNums = str2num(data{cell}{condition+1}) ;        % change epoch string into numbers 
        if isempty(epochNums) ;                                 % if this cell does not have this condition
            Gmean(cell,condition) = nan ;                           % then make a nan 
        else                                                    % otherwise...
            for epoch = 1:length(epochNums) ;                       % for each epoch of that condition
            [Data{condition}{cell}(epoch,:), error] = ITCReadEpoch(epochNums(epoch), 0, fp) ;  % 
            end
            Data_mean{condition}{cell} = mean(Data{condition}{cell},1) ;
            Data_mean{condition}{cell} = Data_mean{condition}{cell} - mean(Data_mean{condition}{cell}(1)) ;
            Imean(cell,condition) = mean(Data_mean{condition}{cell}(5000:10000)) ;
            Gmean(cell,condition) =  Imean(cell, condition)/61.6 ;
        end
    
figure(cell)
if ~isnan(Gmean(cell,condition))
plot(Data_mean{condition}{cell},Color{condition})
hold on
end


    end
end

for a= 1:8
GmeanNorm(a,1:3) = Gmean(a,1:3)/Gmean(a,1) ;
GmeanNorm(a,4:6) = Gmean(a,4:6)/Gmean(a,4) ;
end

figure
plot(1,GmeanNorm(:,2),'b.')
hold on
plot(2,GmeanNorm(:,5),'r.')
plot(1,mean(GmeanNorm(:,2)),'b+')
plot(2,mean(GmeanNorm(:,5)),'r+')

figure
plot(1,Gmean(:,1),'b.')
hold on
plot(2,Gmean(:,4),'r.')
plot(1,mean(Gmean(:,1)),'b+')
plot(2,mean(Gmean(:,4)),'r+')

% put in a struct and igor export

ForIgor.NormDec = GmeanNorm(:,2)' ;
ForIgor.NormInc = GmeanNorm(:,5)' ;
ForIgor.NormDecMean = mean(GmeanNorm(:,2)) ;
ForIgor.NormIncMean = mean(GmeanNorm(:,5)) ;

ForIgor.Dec = Gmean(:,1)' ;
ForIgor.Inc = Gmean(:,4)' ;
ForIgor.DecMean = mean(Gmean(:,1)) ;
ForIgor.IncMean = mean(Gmean(:,4)) ;

% export this structure as a HDF5 which will be read into igor
exportStructToHDF5(ForIgor,'StepResponse','~/Desktop') ;


