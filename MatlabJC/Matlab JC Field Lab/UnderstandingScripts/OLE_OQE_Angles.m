%%%%%%%%%%%%%%%%%% Explanation of script and data requirements
%this script requires that you load a vector called "angles" (one angle for each trial):
%Ntrialsx1 -- these values should be in radians

%you also need to load a matrix "responses" of firing rates on each trial: Ntrials x Ncells

%this script will make the OLE for direction vectors
%and then test that OLE to see how well it predicts directions from firing
%rates

%it willl also train and test the optimal quadratic estimator
%(OQE)

%key outcomes are the OLE and OQE (weightings in the estimators)
%and the angular deviations of each...

%%%%%%%%%%%%%%%%%%%%
%first split the data into a random group for training vs testing
%(that way we don't over-estimate performance by testing on the same data
%we used to train the OLE)
[Ntrials] = length(angles);
permutor = randperm(Ntrials);
train_set = permutor(1:ceil(0.8*Ntrials)); %use 80% for training (these percentages can obv. be changed)
test_set = permutor(ceil(0.8*Ntrials)+1:end); %use the other 20% for testing 
%%%%%%%%%%%%%%%%%%% data parsing



%%%%%%%%%%%%%%%%%%%%% this part trains the OLE
%make the response angles into vectors
xvalues = cos(angles(train_set)); 
yvalues = sin(angles(train_set));
%assemble the vectors into a Ntrials x 2 matrix
stim_train(:,1) = xvalues;
stim_train(:,2) = yvalues;

%now assemble the matrix of neural responses...
resps_train = responses(train_set,:);

%append to this a vector of all ones to use in the regression
%this allows it to learn the constant offset
predictors_train = [resps_train ones(ceil(0.8*Ntrials),1)];

%now do regression to get the estimator coefficients
%this will give a Ncells+1 x 2 matrix of prediction coefficients -- the OLE
OLE = (predictors_train\stim_train) ;
%%%%%%%%%%%%%%%%%%%%% train the OLE
disp('OLE trained')

%%%%%%%%%%%%%%%%this part tests the OLE
%make the response angles into vectors
xvalues = cos(angles(test_set)); 
yvalues = sin(angles(test_set));
%assemble the vectors into a Ntrials x 2 matrix
stim_test(:,1) = xvalues;
stim_test(:,2) = yvalues;

%now assemble the matrix of neural responses...
resps_test = responses(test_set,:);

%append to this a vector of all ones to use in the regression
%this allows it to use the learned constant offset
predictors_test = [resps_test ones(floor(0.2*Ntrials),1)];

%predict the direction vectors -- and normalize them to be unit length
%(just estimating *direction* here, not speed)
predicted_vectors = predictors_test*OLE;
predicted_vectors = diag(1./sqrt(sum(predicted_vectors.^2,2)))*predicted_vectors;

%now compute angle deviations between predicted and true vectors
dot_products = diag(predicted_vectors*stim_test');
angle_deviations_OLE = acos(dot_products); %for each trial, tells you the angle the OLE is off by
%%%%%%%%%%%%%%test the OLE
disp('OLE tested, mean error is (Deg.)')
mean(angle_deviations_OLE)*180/pi


%%%%%%%%%%%%%for the OQE, need an arrays that's Ntrials x (Ncells +
%%%%%%%%%%%%%Ncells^2), containing (for each trial) the concatenation of
%%%%%%%%%%%%%the firing rates, and the products of them.
[NCells] = length(responses(1,:));
disp('making augmented data for OQE; this may take a few minutes')
disp('percentage complete: 0')
 J = triu(ones(NCells),1); % binary array of which elements of product matrix are "not same cell twice"
for ii = 1:Ntrials
    prodmat= responses(ii,:)'*responses(ii,:);
    prods_notsamecells = prodmat(J==1);
    responses_augmented(ii,:) = [responses(ii,:) prods_notsamecells'];
    
    if(mod(ii,floor(Ntrials/10)) == 0)
        fprintf('percentage complete: %s\n',num2str(100*ii/Ntrials));
    end
end
disp('augmented data ready for OQE')

%%%%%%%%%%%%%%%%%%%%% this part trains the OQE
%make the response angles into vectors
xvalues = cos(angles(train_set)); 
yvalues = sin(angles(train_set));
%assemble the vectors into a Ntrials x 2 matrix
stim_train(:,1) = xvalues;
stim_train(:,2) = yvalues;

%now assemble the matrix of neural responses...
resps_train = responses_augmented(train_set,:);

%append to this a vector of all ones to use in the regression
%this allows it to learn the constant offset
predictors_train = [resps_train ones(ceil(0.8*Ntrials),1)  ];

%now do regression to get the estimator coefficients
%this will give a matrix of prediction coefficients -- the OQE
OQE = (predictors_train\stim_train) ;
%%%%%%%%%%%%%%%%%%%%% train the OQE
disp('OQE trained')

%%%%%%%%%%%%%%%%this part tests the OQE
%make the response angles into vectors
xvalues = cos(angles(test_set)); 
yvalues = sin(angles(test_set));
%assemble the vectors into a Ntrials x 2 matrix
stim_test(:,1) = xvalues;
stim_test(:,2) = yvalues;

%now assemble the matrix of neural responses...
resps_test = responses_augmented(test_set,:);

%append to this a vector of all ones to use in the regression
%this allows it to use the learned constant offset
predictors_test = [resps_test ones(floor(0.2*Ntrials),1)];

%predict the direction vectors -- and normalize them to be unit length
%(just estimating *direction* here, not speed)
predicted_vectors = predictors_test*OQE;
predicted_vectors = diag(1./sqrt(sum(predicted_vectors.^2,2)))*predicted_vectors;

%now compute angle deviations between predicted and true vectors
dot_products = diag(predicted_vectors*stim_test');
angle_deviations_OQE = acos(dot_products); %for each trial, tells you the angle the OLE is off by
%%%%%%%%%%%%%%test the OQE
disp('OQE tested, mean angle error is (Deg.)')
mean(angle_deviations_OQE)*180/pi



