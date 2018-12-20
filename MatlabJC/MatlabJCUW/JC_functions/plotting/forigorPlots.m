ploted = plotForIgor(ForIgor) ;

% script for dealing with ForIgor plots

% RASTOR PLOTS
% looking at many types of data in one cell
cell='15' ;

for b=1:5
figure(b)
hold on
r=1 ;

for a = b:10:100
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

% looking at one type of data over many cells
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

% ISI DIS PLOTS
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

% IFRTrain plots
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

% LIF plots
% across one cell
cell='15' ;

for a=1:10
figure(a)
hold on

idx = ['LIFtrace',num2str(a),'DC',cell] ;
idy = ['LIFtrace',num2str(a),'DCnoInh',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idy),'r')
hold on
plot(ForIgor.(idx),'b')
end

end


% LIF vs DC
% across one cell
cell = '15' ;

for a=1:10
figure(a)


idx = ['DCtrace',num2str(a),'DCnoInh',cell] ;
idy = ['LIFtrace',num2str(a),'DCnoInh',cell] ;
if isfield(ForIgor,idy)
plot(ForIgor.(idy),'r')
hold on
plot(ForIgor.(idx),'b')
end

end





