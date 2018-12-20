% function [covar,EigVal,EigVec] = contPCA(response,stimulus,lag,minCovPnts) ;

% % analog covar matrix
% uniqueResponse = unique(response(:)) ; % unique responses
% 
% round = 1 ;
% num_responseValues = 0 ;
% covar = zeros(lag,lag) ;
% for trial = 2:length(uniqueResponse) ; % for each unique response
% 
% [responseRow,responseColumn] = find(response>uniqueResponse(trial-round) & response<=uniqueResponse(trial)); % find those indicies of the unique respones value
% 
% responseColumn = responseColumn(responseColumn>lag/2 & responseColumn<length(response)-lag/2) ; % use only those indicies whithin the range you can still lag the stim
% responseRow = responseRow(responseColumn>lag/2 & responseColumn<length(response)-lag/2) ;
% 
% if length(responseRow)<minCovPnts ; % if you don't have enough points to estimate covariance
%     round = round+1 ; % run another round
% else % otherwise calculate covariance and wieght by response
%     num_responseValues = num_responseValues + 1 ; % number of different covariance matricies calculated
%     
%     responseValue =mean(response(sub2ind(size(response),responseRow,responseColumn))) ; % mean of all values
%     
%     stimSet = 0 ;
%     for a = 1:length(responseColumn) ; % for each response pnt
%         for Tau = -lag/2:lag/2-1 ;
%             stimSet = stimSet+stimulus(responseRow(a),responseColumn(a)+Tau) ;
%         end
%     end
%     stimSetMean = stimSet/(length(responseColumn)*length(-lag/2:lag/2-1)-1) ; % mean of selected stim set
%     
%     for a = 1:length(responseColumn) ; % for each response pnt
%         for Tau = -lag/2:lag/2-1 ;
%             for Tau2 = -lag/2:lag/2-1 ;   
%                 covar(Tau+lag/2+1,Tau2+lag/2+1) = covar(Tau+lag/2+1,Tau2+lag/2+1) + responseValue*(stimulus(responseRow(a),responseColumn(a)+Tau)-stimSetMean)*(stimulus(responseRow(a),responseColumn(a)+Tau2)-stimSetMean) ;      
%             end
%         end
%     end
% 
%     covar = covar/(length(responseColumn)*length(-lag/2:lag/2-1)^2) ; % covariance for this response value
%     round = 1 ;
% end
% end
% 
%     
% figure(1);
% mesh(covar);
% 
% [EigVec, EigVal] = eig(covar);
% figure(2);
% semilogy(abs(EigVal), 'o');



function [covar,EigVal,EigVec] = contPCA(response,stimulus,lag) ;
% reproduced from FR? 

covar = zeros(lag,lag) ;

for trial = 1:size(response,1) ; % for each response

for t=lag/2+1:length(response)-lag/2 ;
    for Tau = -lag/2:lag/2-1 ;
        for Tau2 = -lag/2:lag/2-1 ;   
            covar(Tau+lag/2+1,Tau2+lag/2+1) = covar(Tau+lag/2+1,Tau2+lag/2+1) + response(trial,t)*stimulus(trial,t-Tau)*stimulus(trial,t-Tau2) ;
        end
    end
end

end
covar = covar/size(response,1) ;

figure(1);
mesh(covar);

[EigVec, EigVal] = eig(covar);
figure(2);
semilogy(abs(EigVal), 'o');

% [EigVec, EigVal] = eigs(covar, 2);
% 
% figure(3);
% plot(EigVec);