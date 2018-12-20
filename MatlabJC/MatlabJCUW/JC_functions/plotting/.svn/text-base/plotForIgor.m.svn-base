ploted = plotForIgor(ForIgor) ;

% script for dealing with ForIgor plots

ploted = 1 ; % output (mean nothing, all I want is the graphs)


%% RASTOR PLOTS
% plot rastor of spiketimes
for cell = 1:20 ; % for each possible cell
cell = num2str(cell)
    
for b=1:5 ; % for each possible light stim
figure(b)
hold on
r=1 ;

for a = b:10:100 % the same light stim
id = ['rastor',num2str(a),'CA',cell] ;
if isfield(ForIgor,id)
plot(ForIgor.(id),r,'k.')
r=r+1 ;
end
end

for a = b:10:100
id = ['rastor',num2str(a),'IC',cell] ;
if isfield(ForIgor,id)
plot(ForIgor.(id),r,'y.')
r=r+1 ;
end
end

for a = b:5:100
id = ['rastor',num2str(a),'DC',cell] ;
if isfield(ForIgor,id)
plot(ForIgor.(id),r,'b.')
r=r+1 ;
end
end

for a = b:5:100
id = ['rastor',num2str(a),'DCnoInh',cell] ;
if isfield(ForIgor,id)
plot(ForIgor.(id),r,'r.')
r=r+1 ;
end
end

for a = b:5:100
id = ['rastor',num2str(a),'DCApb',cell] ;
if isfield(ForIgor,id)
plot(ForIgor.(id),r,'g.')
r=r+1 ;
end
end

for a = b:5:100
id = ['rastor',num2str(a),'DCnoExc',cell] ;
if isfield(ForIgor,id)
plot(ForIgor.(id),r,'c.')
r=r+1 ;
end
end

axis([0 8 -10 r+10])
end
end

clearvars -except ForIgor

%% looking at one type of data over many cells
c=['b','r','g','y','c','k','b','r','k','b','r','g','y','c','k','b','r','k'] ;
r=1
for a=1:10
    for b=1:10:100
        id = ['rastor',num2str(b),'CA',num2str(a)] ;
        if isfield(ForIgor,id)
            plot(ForIgor.(id),r,[c(a),'.'])
            r=r+1 ;
            hold on
        end    
        end
end



%% ISI DIS PLOTS
% looking at many types of data across one cell
for cell = [1,4,6:9,13] ;
figure
hold on

idy = ['ISIdis','CA1'] ;
idx = ['ISIbins','CA1'] ;
if isfield(ForIgor,idx)
plot(ForIgor.(idx),ForIgor.(idy),'k')
end

idy = ['ISIdis','CA4'] ;
idx = ['ISIbins','CA4'] ;
if isfield(ForIgor,idx)
plot(ForIgor.(idx),ForIgor.(idy),'k--')
end

idy = ['ISIdis','DC',num2str(cell)] ;
idx = ['ISIbins','DC',num2str(cell)] ;
if isfield(ForIgor,idx)
plot(ForIgor.(idx),ForIgor.(idy),'b')
end


idy = ['ISIdis','DCnoInh',num2str(cell)] ;
idx = ['ISIbins','DCnoInh',num2str(cell)] ;
if isfield(ForIgor,idx)
plot(ForIgor.(idx),ForIgor.(idy),'g')
end


idy = ['ISIdis','DCnoExc',num2str(cell)] ;
idx = ['ISIbins','DCnoExc',num2str(cell)] ;
if isfield(ForIgor,idx)
plot(ForIgor.(idx),ForIgor.(idy),'r')
end


idy = ['ISIdis','DCApb',num2str(cell)] ;
idx = ['ISIbins','DCApb',num2str(cell)] ;
if isfield(ForIgor,idx)
plot(ForIgor.(idx),ForIgor.(idy),'c')
end

end

%% looking at isi across many cells
figure

for cell = [1:10] ;
idy = ['ISIdis','CA',num2str(cell)] ;
idx = ['ISIbins','CA',num2str(cell)] ;
if isfield(ForIgor,idx)
plot(ForIgor.(idx),ForIgor.(idy),'k')
hold on
end
end

%% NN plots

for cell = [7:9,13] ;
figure
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on

idy = ['NNcdN','CA',num2str(cell)] ;
idy2 = ['NNcd','CA',num2str(cell)] ; 
idx = ['NNbins','CA',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'k')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'k')
end

idy = ['NNcdN','DC',num2str(cell)] ;
idy2 = ['NNcd','DC',num2str(cell)] ; 
idx = ['NNbins','DC',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'b')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'b')
end

idy = ['NNcdN','DCnoInh',num2str(cell)] ;
idy2 = ['NNcd','DCnoInh',num2str(cell)] ; 
idx = ['NNbins','DCnoInh',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'g')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'g')
end

idy = ['NNcdN','DCnoExc',num2str(cell)] ;
idy2 = ['NNcd','DCnoExc',num2str(cell)] ; 
idx = ['NNbins','DCnoExc',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'r')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'r')
end

idy = ['NNcdN','DCApb',num2str(cell)] ;
idy2 = ['NNcd','DCApb',num2str(cell)] ; 
idx = ['NNbins','DCApb',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'c')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'c')
end

end

% reverse nn plots

for cell = [1,4,6:9,13] ;
figure
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on

idy = ['NNcdRN','CA',num2str(cell)] ;
idy2 = ['NNcdR','CA',num2str(cell)] ; 
idx = ['NNbinsR','CA',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'k')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'k')
end

idy = ['NNcdRN','DC',num2str(cell)] ;
idy2 = ['NNcdR','DC',num2str(cell)] ; 
idx = ['NNbinsR','DC',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'b')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'b')
end

idy = ['NNcdRN','DCnoInh',num2str(cell)] ;
idy2 = ['NNcdR','DCnoInh',num2str(cell)] ; 
idx = ['NNbinsR','DCnoInh',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'g')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'g')
end

idy = ['NNcdRN','DCnoExc',num2str(cell)] ;
idy2 = ['NNcdR','DCnoExc',num2str(cell)] ; 
idx = ['NNbinsR','DCnoExc',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'r')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'r')
end

idy = ['NNcdRN','DCApb',num2str(cell)] ;
idy2 = ['NNcdR','DCApb',num2str(cell)] ; 
idx = ['NNbinsR','DCApb',num2str(cell)] ;
if isfield(ForIgor,idx)
    subplot(2,1,1)
    plot(ForIgor.(idx),ForIgor.(idy),'c')
    subplot(2,1,2)
    plot(ForIgor.(idx),ForIgor.(idy2),'c')
end

end


%% IFRTrain plots
% across one cell
cell='1' ;

for a=1:5
figure(a)
hold on

idx = ['time','CA1'] ;
idy = ['ifrTrain',num2str(a),'CA1'] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'k')
end

idx = ['time','IC',cell] ;
idy = ['ifrTrain',num2str(a),'IC',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'y')
end

idx = ['time','DC',cell] ;
idy = ['ifrTrain',num2str(a),'DC',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'b')
end

idx = ['time','DCnoInh',cell] ;
idy = ['ifrTrain',num2str(a),'DCnoInh',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'g')
end

idx = ['time','DCnoExc',cell] ;
idy = ['ifrTrain',num2str(a),'DCnoExc',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'r')
end

idx = ['time','DCApb',cell] ;
idy = ['ifrTrain',num2str(a),'DCApb',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'c')
end

end

%% plot light to G normalized transfer functions

for cell=[1] ;

figure    
hold on
cell = num2str(cell) ;
   
idx = ['xTFun','CA',cell] ;
idy = ['TFunN','CA',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'k*-')
end 
    
idx = ['xTFun','Exc',cell] ;
idy = ['TFunN','Exc',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'g*-')
end

idx = ['xTFun','Inh',cell] ;
idy = ['TFunN','Inh',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'r*-')
end

idx = ['xTFun','InhApb',cell] ;
idy = ['TFunN','InhApb',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idx),ForIgor.(idy),'c*-')
end

a=gca ;
set(a,'xScale','log')
set(a,'yScale','log')

end


