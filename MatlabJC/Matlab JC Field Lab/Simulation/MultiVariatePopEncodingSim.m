% simulation to explore signal and noise correlation integration

% JC 3/19/2018 

% time and 'trials'
ns = 6 ; % number of discreet stim
nt = 1000 ; % number of trials

% stimulus tuning
TcS_type = {'linear','flat'} ; % options: linear, pdf, cdf, flat
TcS_params = {3,3} ; % parameters for tuning curves

% distractor tuning
TcD_type = {'randPlus','randPlus'} ; % options: linear, pdf, cdf, flat
TcD_params = {.01,.01} ; % parameters for tuning curves

% variables
S = repmat([1:ns], nt, 1); % signal (each column is a repeat stim value X) 
D1 = normrnd(ns/2,ns/2,nt,ns) ; % distractor 1

NumCells = length(TcS_type) ;

% Responses (cell,trial,stim)
corrD1 = D1(randperm(nt),randperm(ns)) ; %for use with 'randPlus'

for N=1:NumCells ; % for each cell
    % response to S
    if strcmp(TcS_type{N},'linear')
        Rs(N,:,:) = S*TcS_params{N} ;
    elseif strcmp(TcS_type{N},'flat')
        Rs(N,:,:) = ones(nt,ns)*TcS_params{N} ;
    elseif strcmp(TcS_type{N},'pdf')
        Rs(N,:,:) = normpdf(S,TcS_params{N},TcS_params{N}) ;
        Rs(N,:,:) = ns*Rs(N,:,:)/normpdf(TcS_params{N},TcS_params{N},TcS_params{N}) ; % max same as linear model
    elseif strcmp(TcS_type{N},'pdfV1')
        Rs(N,:,:) = normpdf(S,TcS_params{N},2) ;
        Rs(N,:,:) = ns*Rs(N,:,:)/normpdf(TcS_params{N},TcS_params{N},2) ; % max same as linear model
    elseif strcmp(TcS_type{N},'cdf')    
        Rs(N,:,:) = ns*normcdf(S,TcS_params{N},TcS_params{N}) ;
    elseif strcmp(TcS_type{N},'cdfV1')    
        Rs(N,:,:) = ns*normcdf(S,TcS_params{N},2) ;
    elseif strcmp(TcS_type{N},'randPlus')
        % not finished
    end
    Rs(N,:,:) = round(Rs(N,:,:)) ; % spike number needs to be int
    
    % response to D1
    if strcmp(TcD_type{N},'linear')
        Rd(N,:,:) = D1*TcD_params{N} ;
    elseif strcmp(TcD_type{N},'flat')
        Rd(N,:,:) = ones(nt,ns)*TcD_params{N} ;
    elseif strcmp(TcD_type{N},'pdf')
        Rd(N,:,:) = normpdf(D1, TcD_params{N}, TcD_params{N}) ;
        Rd(N,:,:) = ns*Rd(N,:)/normpdf(TcD_params{N},TcD_params{N},TcD_params{N}) ;
    elseif strcmp(TcD_type{N},'cdf')    
        Rd(N,:,:) = ns*normcdf(D1, TcD_params{N}, TcD_params{N}) ;
    elseif strcmp(TcD_type{N},'rand')
        Rd(N,:,:) = D1(randperm(nt),randperm(ns))*TcD_params{N} ;
    elseif strcmp(TcD_type{N},'randPlus')
        Rd(N,:,:) = corrD1*(1-TcD_params{N})+ D1(randperm(nt),randperm(ns))*TcD_params{N} ;
    end
    Rd(N,:,:) = round(Rd(N,:,:)) ; % spike number needs to be int
    
    % response
    R(N,:,:) = Rs(N,:,:).*Rd(N,:,:) ;
end
    
% information between individual cell and S
Hs = log2(ns) ;
for N=1:NumCells ; % for each cell 
    Is(N)=get_MI_from_resp(squeeze(Rs(N,:,:)),[1:ns],1,[1:nt]) ;
    Id(N)=get_MI_from_resp(squeeze(Rd(N,:,:)),[1:ns],1,[1:nt]) ;
    Ir(N)=get_MI_from_resp(squeeze(R(N,:,:)),[1:ns],1,[1:nt]) ;
end

% information with two cells
for N=1:NumCells ; % for each cell
    Rpair = {squeeze(R(1,:,:)),squeeze(R(N,:,:))} ;
    Ipair(N) = get_MI_from_pair_resp(Rpair,[1:ns],1,[1:nt]) ;
    Redundancy(N) = (Ipair(N) - (Ir(1) + Ir(N)))/Ipair(N) ;
    Atemp = squeeze(Rd(1,:,:)) ;
    Btemp = squeeze(Rd(N,:,:)) ;
    CovRd1with(N) = corr(Atemp(:),Btemp(:)) ;
    Atemp = squeeze(Rs(1,:,:)) ;
    Btemp = squeeze(Rs(N,:,:)) ;
    CovRs1with(N) = corr(Atemp(:),Btemp(:)) ;
end

% Ishuff
Rmean = squeeze(mean(R,2)) ; % average over trials
RmeanDelta = diff(Rmean,1,2) ;

for st = 1:ns-1 ; % for each stim change
    RCov = cov(squeeze(R(:,:,st))') ;
    Fi(st) = RmeanDelta(:,st)'*inv(RCov)*RmeanDelta(:,st) ;
    RCovInd = RCov.*eye(NumCells) ;
    FiInd(st) = RmeanDelta(:,st)'*inv(RCovInd)*RmeanDelta(:,st) ;
end

Ishuff = sum(Fi) - sum(FiInd) ;
IshuffRel = Ishuff/sum(Fi) ;
  
% Idiag - Averbeck et al 2006 ;
for st = 1:ns-1 ; % for each stim change
    RCov = cov(squeeze(R(:,:,st))') ;
    RCovInd = RCov.*eye(NumCells) ;
    Denom = RmeanDelta(:,st)'*inv(RCovInd)*RCov*inv(RCovInd)*RmeanDelta(:,st) ;
    FiDiag(st) = FiInd(st)^2/Denom ; 
end

% change Fi, FiInd, FiDiag into percent correct classification
PcFi = 1 - erfc(sqrt(Fi)/2) ;
PcFiInd = 1 - erfc(sqrt(FiInd)/2) ;
PcFiDiag = 1 - erfc(sqrt(FiDiag)/2) ;

% Idiag - OLE/OQE decoding all stim;
% shuffled trials
Rshuff = nan(size(R)) ;
for N=1:NumCells ; % for each cell
    for st = 1:ns ; % for each stim change
        Rshuff(N,:,st) = R(N,randperm(nt),st) ;
    end
end

% fit/test OLE 
Sreshape = reshape(S,[ns*nt,1]) ;
Rreshape = reshape(R, [NumCells,ns*nt])' ;
W = Rreshape\Sreshape ; % R(cell,trial,stim)

OleSest = Rreshape*W ;
OleMse = mean((Sreshape-OleSest).^2) ;

% fit/test OLE with Rshuff
Sreshape = reshape(S,[ns*nt,1]) ;
RshuffReshape = reshape(Rshuff, [NumCells,ns*nt])' ;
Wshuff = RshuffReshape\Sreshape ; % R(cell,trial,stim)

OleShuffSest = Rreshape*Wshuff ;
OleShuffMse = mean((Sreshape-OleShuffSest).^2) ;

% fit/test OQE
J = triu(ones(NumCells),1); % binary array of which elements of product matrix are "not same cell twice"
RcrossProds = nans(nt*ns,sum(J(:))) ; %prep mat
for tr = 1:nt*ns ; % for each trial/stim
    temp = Rreshape(tr,:)'*Rreshape(tr,:) ; % cross products
    RcrossProds(tr,:) = temp(J==1) ;
end
Rall = [Rreshape,RcrossProds] ;   
Woqe = Rall\Sreshape ;

OqeSest = Rall*Woqe ;
mean((Sreshape-OqeSest).^2)

% fit/test OQE with Rshuff
J = triu(ones(NumCells),1); % binary array of which elements of product matrix are "not same cell twice"
RcrossProds = nans(nt*ns,sum(J(:))) ; %prep mat
for tr = 1:nt*ns ; % for each trial/stim
    temp = RshuffReshape(tr,:)'*RshuffReshape(tr,:) ; % cross products
    RcrossProds(tr,:) = temp(J==1) ;
end
RallShuff = [RshuffReshape,RcrossProds] ;   
WoqeShuff = RallShuff\Sreshape ; % train on shuffle

OqeShuffSest = Rall*WoqeShuff ; % test on actual
mean((Sreshape-OqeShuffSest).^2)

% fit/test mod OLE
% Rcov mean
RCovMean = zeros(NumCells) ;
for st = 1:ns ; % for each stim change
    RCovMean = RCovMean + cov(squeeze(R(:,:,st))') ;
end
RCovMean = RCovMean/ns ;

modA = -alpha*(RCovMean*R)'*W ;
OleSest + Rreshape*W



% figure
figure % Idiag - OLE
subplot(2,1,1) % S and Sest
plot(OleShuffSest,'b')
hold on
plot(OleSest,'r:')
plot(Sreshape,'k-')

subplot(2,1,2) % W
plot(Wshuff,'b')
hold on
plot(W,'r')

figure % Idiag - OQE
subplot(2,1,1) % S and Sest
plot(OqeShuffSest,'b')
hold on
plot(OqeSest,'r:')
plot(Sreshape,'k-')

subplot(2,1,2) % W
plot(WoqeShuff,'b')
hold on
plot(Woqe,'r')

figure
subplot(3,1,1) % info as function of distractor covariance
plot(CovRd1with,Redundancy,'*')
hold on
xlabel('distractor covariance')
ylabel('Redundancy')

subplot(3,1,2) % info as function of signal covariance
plot(CovRs1with,Redundancy,'*')
hold on
xlabel('signal covariance')
ylabel('Redundancy')

subplot(3,1,3) % info as function of signal covariance
plot(CovRs1with,Ipair/Hs,'*')
hold on
xlabel('signal covariance')
ylabel('fraction signal entropy')
