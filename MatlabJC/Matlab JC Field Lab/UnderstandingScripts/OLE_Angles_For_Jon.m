%%%%%%%%%%%%%%%%%% Explanation of script and data requirements
%this script requires that you load a vector called "angles" (one angle for each trial):
%Ntrialsx1 -- these values should be in radians

%you also need to load a matrix "responses" of firing rates on each trial: Ntrials x Ncells

%this script will make the OLE for direction vectors
%and then test that OLE to see how well it predicts directions from firing
%rates

%outcome of the analysis will be the fraction of explained variance (R^2)
%%%%%%%%%%%%%%%%%%%%


%first split the data into a random group for training vs testing
%(that way we don't over-estimate performance by testing on the same data
%we used to train the OLE)
[Ntrials] = length(angles);
permutor = randperm(Ntrials);
train_set = permutor(1:ceil(0.8*Ntrials)); %use 80% for training (these percentages can obv. be changed)
test_set = permutor(ceil(0.8*Ntrials)+1:end); %use the other 20% for testing 

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
angle_deviations = acos(dot_products); %for each trial, tells you the angle the OLE is off by

%and use this to get the fraction of explained variance in the angles (FEV)
%which is 1 - MSE/Var, where MSE is mean squared error of estimator
%and var is the variance in the angles
FEV = 1 - mean(angle_deviations.^2)/var(angles(test_set));
%%%%%%%%%%%%%%test the OLE