function [midAngle, sets, I] = angleDiscriminationTask_VariableN(T,N,noiseCorrProp,noiseAmp,corrMatrix,angleDelta,Npositions,poissNoise,Nperpopulation,Ntrials)
%T, tuning curves (from makeTuningCurves)
%N number of cells to consider at one in information calculation
%noiseCorrProp, proportion of noise that is shared
%corrMatrix, sets which cells are correlated (0 if uncorrelated)
%noiseAmp, noise amplitude
%angleDelta, for discrimination task, must be integer divisible by 2
%Npositions, number of positions around the circle to test, must be divisor
%of 360
%poissNoisem, make noise variance proportional to mean response
%Nperpopulation, number of cells with independent noise in each population
%(to reduce noise by sqrt(N)
%Ntrials, number of trials in simulation

phaseOffset = 360./Npositions;
Ncells = size(T,1);
sets = combnk(1:Ncells,N);
Nsets = size(sets,1);

midAngle = zeros(1,Npositions);

%set up correlated noise streams
randStreams = setdiff(unique(corrMatrix),0);
L = length(randStreams);
corrNoise = zeros(L,Ntrials*2);
for i=1:L
    randn('seed',randStreams(i));
    corrNoise(i,:) = randn(1,Ntrials*2);
end

%information
I = zeros(Nsets,Npositions);

for p=1:Npositions
    midAngle(p) = (p-1)*phaseOffset+1;
    A1 = midAngle(p)-angleDelta/2;
    if A1 < 1, A1 = A1+360; end
    A2 = midAngle(p)+angleDelta/2;
    if A2 > 360, A2 = A2-360; end
    
    for i=1:Nsets        
        %make responses 
        c = sets(i,:); %cell identities for this set        
        R = zeros(N,Ntrials); %responses
        Nc = zeros(N,Ntrials*2); %correlated noise
        Nuc = zeros(N,Ntrials*2); %uncorrelated noise
        Nfull1 = zeros(N,Ntrials); %total noise, angle 1
        Nfull2 = zeros(N,Ntrials); %total noise, angle 2
        for n=1:N
            %responses
            R1(n,:) = ones(1,Ntrials).*T(c(n),A1);
            R2(n,:) = ones(1,Ntrials).*T(c(n),A2);
            
            %correlated noise
            r = corrMatrix(c(n)); %which noise seed
            if r>0
                Nc(n,:) = noiseCorrProp.*corrNoise(r,:);
            end
            
            %random noise
            a = clock+n; randn('seed',a(end));
            randNoise = randn(1,Ntrials*2);
            Nuc(n,:) = (1-noiseCorrProp).*randNoise;
            
            if poissNoise
                %make the noise variance equal to the mean
                Nfull1(n,:) = (Nc(n,1:Ntrials) + Nuc(n,1:Ntrials)).* sqrt(R1(n,1)) ./ sqrt(Nperpopulation);
                Nfull2(n,:) = (Nc(n,Ntrials+1:end) + Nuc(n,Ntrials+1:end)).* sqrt(R2(n,1)) ./ sqrt(Nperpopulation);
            else
                %scale by noiseAmp
                Nfull1(n,:) = (Nc(n,1:Ntrials) + Nuc(n,1:Ntrials)).* noiseAmp ./ sqrt(Nperpopulation);
                Nfull2(n,:) = (Nc(n,Ntrials+1:end) + Nuc(n,Ntrials+1:end)).* noiseAmp ./ sqrt(Nperpopulation);
            end
            
            %add noise
            R1(n,:) = R1(n,:) + Nfull1(n,:);
            R2(n,:) = R2(n,:) + Nfull2(n,:);
        end %for each member of set
        
        %calculate mutual information
                
        P = zeros(2^(N-1),2); %P(r,s), two stimuli are angle 1 and angle 2
        
        %make predicates to get probabilities
        bitResp = cell(2^(N-1),2); %bits represention the reponse we are looking for (like '101')        
        for z=1:2^(N-1)
            bitResp{z,1} = dec2bin(z-1,N);
        end
        ind = 2^N-1;
        for z=1:2^(N-1)
            bitResp{z,2} = dec2bin(ind,N);
            ind = ind-1;
        end
        
        for j=1:size(bitResp,1)
            for k=1:2                
                predicate = '';
                b = bitResp{j,k};
                for n=1:N
                    if str2num(b(n))
                        predicate = [predicate 'R1(' num2str(n) ',:)>R2(' num2str(n) ',:)'];
                    else
                        predicate = [predicate 'R1(' num2str(n) ',:)<=R2(' num2str(n) ',:)'];
                    end                    
                    if n==N
                        predicate = [predicate ';'];
                    else
                        predicate = [predicate '&'];
                    end
                end
                %predicate
                P(j,k) = sum(eval(predicate));
            end            
        end
        
        %normalize and do the MI calculation with P
        P = P./Ntrials;        
        Ps = [0.5 0.5]; %equal prob of 2 angles
        I(i,p) = mutualInformation(P, Ps);
        
    end %for each set 
end %for each position





