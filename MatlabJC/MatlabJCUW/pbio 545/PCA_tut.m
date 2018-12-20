
    %-------------------------------------------------------------
% Principle Components Analysis Tutorial for NeuBeh/PBIO 545
%
% This tutorial introduces Principle Components Analysis.  The aims are:
% (1) introduce dimensional reduction and its importance
% (2) introduce the covariance and covariance matrix
% (3) show how the eigensystem of the covariance matrix helps rank dimensions of
%	  a probability distribution by the amount of variance they capture
% (4) illustrate these concepts with a real practical application
%
% This builds on two key concepts introduced earlier in the course: probability 
% distributions and eigensystems.  If you are uncomfortable with either of those
% take another pass through the appropriate sections of the stochastic processes
% and/or linear algebra tutorials.
%
% Created 2/03 Fred Rieke
% Revised 2/07 FMR
%       cleanup up a bit, elaborated a few sections
%-------------------------------------------------------------

% define axis fonts
set(0, 'DefaultAxesFontName','Helvetica')
set(0, 'DefaultAxesFontSize', 16)

%-------------------------------------------------------------
% Part 1: Motivation for dimensional reduction
%
%	What do we mean by dimensional reduction and why should you care
% about it?  With finite experimental data we always face a tradeoff between
% trying to identify all the important structure in the data while not being
% misled by randomness introduced because we don't have an infinite number of 
% samples.  A common example is construction of a post-stimulus time histogram 
% for spike responses to a repeated stimulus: in choosing the time bin for the 
% histogram we face a tradeoff of smoothing over temporal structure in the spike 
% train and avoiding lots of spurious structure due to finite data.  

% Here's an example of a set of spike responses from a retinal ganglion cell to a 
% dim flash of light:  
load rgc-spike-response

% plot several examples
figure(1);
clf;
plot(0);
hold on;
for resp=1:10
	plot([-1:.001:2.999], RGCSpikes(resp, :)/1.2 + resp)
end
hold off
xlabel('time (sec)')
ylabel('trial')

% This data is discretized into 1 msec time bins - so we have already thrown some
% information about spike times away.  Still a psth without any extra smoothing looks
% pretty ugly:

figure(1)
subplot(1, 2, 1)
plot([-1:.001:2.999], mean(RGCSpikes))
xlabel('time (sec)')
ylabel('spike probability')

% You can at least see there is an extra density of spikes after the flash
% (at time 0), but you also have some sense that most of the structure in
% the histogram is just noise.  In effect we are retaining too many
% variables - in this case each of the 4000 time bins that specify the
% spike times.  This is equivalent to representing each spike train in a
% 4000 dimensional space, where each dimension corresponds to one time bin.
% Let's resample the data to produce larger time bins and repeat

ResampleFactor = 100;		% 100 msec bins
clear ResampledRGCSpikes;
for resp=1:size(RGCSpikes, 1)
	ResampledRGCSpikes(resp, :) = decimate(RGCSpikes(resp, :), ResampleFactor, 1);
end
figure(1)
subplot(1, 2, 2)
plot([-1:.001*ResampleFactor:2.999], mean(ResampledRGCSpikes));
xlabel('time (sec)')
ylabel('spike probability')

% This looks better.  Now we are in effect using a 40 dimensional space,
% with each dimension a 100 msec time bin.  But we might start to worry
% that we are smoothing too much and losing some real structure in the
% data.  Part of the problem here is unavoidable: with finite data there is
% only so much structure we can identify and separate from sampling errors.
% But we would like to make the best use possible of our finite data - i.e.
% we would like to choose the low dimensional space used to represent the
% data carefully so that the structure of the data is captured as
% effectively as possible.

% You face this issue every time you do an experiment.  By choosing to
% filter and sample your data at a specific frequency you are choosing what
% parts of the data you consider important and what parts unimportant.  

% PCA is about finding a relatively low dimensional representation of a set
% of data - and doing so systematically rather than using some ad hoc
% procedure.  In the example above, the strategy of smoothing over bins to
% produce a better looking psth is somewhat arbitrary and is not guided by
% the structure of the cell's responses itself.  What PCA will do is
% provide a set of components (think of this as a set of axes in the space
% the data is represented in - our 4000 dimensional space above).
% Importantly, PCA tells us how much structure in the data is captured by
% each dimension.  This will, hopefully, become more concrete below.

%-------------------------------------------------------------
% Part 2: Covariance and covariance matrix
%
%	Consider a signal characterized by two variables, x and y.  We have
% multiple samples of x,y pairs of the signal.  The covariance matrix for x
% and y is defined as:
%	Covar = [<xx>-<x><x>  <xy>-<x><y>
%			 <xy>-<x><y>  <yy>-<y><y>]
% where <x> is the average of x across samples, etc.  The diagonal elements
% contain the variances of x (upper left) and y (lower right).  The off
% diagonal terms contain the covariance between x and y - i.e. how much do
% x and y change together.  Note that each element is corrected for the
% expectation based on the mean values alone - i.e. we take <xy> and
% subtract what we expect from the averages of x and y alone - <x><y> -
% assuming x and y are uncorrelated.
%
% This is closely related to r values and the correlation matrix.  To
% compute the correlation between x and y, you would take 
%   (<xy>-<x><y>) / sqrt(<xx><yy>)
% in other words you take the off diagonal term of the covariance matrix
% above and normalize it by the square root of the product of the diagonal
% terms.  In general the i,j element of the correlation matrix is
%		CorM(i,j) = Cov(i,j)/sqrt[Cov(i,i)Cov(j,j)]
% or alternatively the covariance matrix is
%		Cov(i,j) = Cor(i,j) * sqrt[var_i var_j]

% Let's look at some examples.
%
% Example 1: x and y independent

% each row of Dist contains a separate x-y pair, x is the first column, y the second
clear Dist;
Dist(:, 1) = normrnd(0,0.3,10000, 1);			% normal distribution, sd 0.3		
Dist(:, 2) = normrnd(0,1,10000, 1);				% normal distribution, sd 1

% plot sample points in x-y plane
figure(1)
clf
subplot(1, 2, 1)
plot(Dist(:, 1), Dist(:, 2), '.');
axis([-3 3 -3 3])
axis square
xlabel('x')
ylabel('y')

% construct the covariance matrix explicitly
varx = mean(Dist(:, 1).^2) - mean(Dist(:, 1))^2;		% <xx> - <x><x>
vary = mean(Dist(:, 2).^2) - mean(Dist(:, 2))^2;		% <yy> - <y><y>
covarxy = mean(Dist(:, 1) .* Dist(:, 2)) - mean(Dist(:, 1)) * mean(Dist(:, 2));		% <xy> - <x><y>

Covar = [varx 		covarxy
		 covarxy 	vary]

% alternatively Matlab has a covariance command

cov(Dist)

% The diagonal of the covariance matrix is the variance of x and y,
% and the off diagonal elements are small (and in fact are not 0 only 
% due to the finiteness of our sample set - rerun the code above several 
% times to convince yourself of this). 

% Compare this to the correlation matrix

corrcoef(Dist)

% Again we see structure along the diagonal (now nicely normalized) and not
% much off diagonal.

% The lack of covariance (and correlation) is what we might expect from the
% shape of the plotted distribution.  The spread of points along the x axis
% appears uncorrelated with the spread along the y axis - or in other
% words, if I told you the x axis value, you could not predict the y axis
% value.  

% Example 2: x and y correlated

clear Dist;
Dist(:, 1) = normrnd(0,0.3,10000, 1);			
Dist(:, 2) = normrnd(0,1,10000, 1) + 1.5 * Dist(:, 1);

figure(1)
subplot(1, 2, 2)
plot(Dist(:, 1), Dist(:, 2), '.');
axis([-3 3 -3 3])
axis square
xlabel('x')
ylabel('y')

% construct the covariance matrix 
varx = mean(Dist(:, 1).^2) - mean(Dist(:, 1))^2;
vary = mean(Dist(:, 2).^2) - mean(Dist(:, 2))^2;
covarxy = mean(Dist(:, 1) .* Dist(:, 2)) - mean(Dist(:, 1)) * mean(Dist(:, 2));

Covar = [varx 		covarxy
		 covarxy 	vary]
cov(Dist)

% Now note that the distribution shows a correlation between x and y - i.e.
% the distribution is angled in the x-y plane.  Correspondingly the
% covariance has off diagonal elements.  

% Again compare this to the correlation matrix

corrcoef(Dist)

% You can think of lots of cases where x and y might be correlated.  For
% example, x might be fraction of lectures Mike gives, and y might be
% fraction of lectures containing superfluous PowerPoint `features' like
% lightening bolts, clapping, etc (he's toned way down too).  Another
% example is two cells with common input and hence correlated spike trains:
% x could be the firing rate in one cell and y the firing rate in the
% other.  If the cells are correlated, when cell x generates a high firing
% rate, so does cell y.  In this case the value of x gives you some ability
% to predict y.

% Here is another use of the covariance matrix.  Let's say we want to take
% our distribution Dist above and "spherize" it (this is what happens when
% you let the statisticians take over language - you get verbs like "to
% sphere").  What that means is we want to create a distribution with new
% axes x' and y' where the two are uncorrelated, and the standard deviation
% along each axis is 1.  We can do that using the inverse of the SQUARE
% ROOT of the covariance matrix.  Why the square root?  Remember the
% entries of the covariance are variances or covariances - e.g. the upper
% right element is the variance of x.  So to get the units right we need
% something that is more like a standard deviation - hence the square root.

% We considered just the scaling part in the linear algebra tutorial, so
% let's go back over that.  Consider distribution of data points
% characterized by two parameters:

clear Dist;
Dist(1, :) = normrnd(0,0.5,1, 5000);
Dist(2, :) = normrnd(0,1.5,1, 5000);

% plot points
figure(1); clf
subplot(1, 2, 1)
plot(Dist(1, :), Dist(2, :), '.');
axis([-5 5 -5 5])
axis square
xlabel('x')
ylabel('y')

% To measure distance of a given data point from the origin in standard
% deviations, we rescale the axes by the standard deviation:

M = [0.5 0 				% matrix with sd along axis
	 0 1.5]

NewDist = inv(M) * Dist;	% rescale distribution using inverse of M
subplot(1, 2, 2)			
plot(NewDist(1, :), NewDist(2, :), '.');
axis([-5 5 -5 5])
axis square
xlabel('x / \sigma_x')
ylabel('y / \sigma_y')

% now the sd along each axis is 1 and the distribution is symmetrical, and
% moving a distance 1 along either the x or the y axis moves us 1 sd

% OK, now back to our problem.  We have the additional issue that our
% distribution was not lined up with the x and y axes:

clear Dist;
Dist(:, 1) = normrnd(0,0.3,10000, 1);			
Dist(:, 2) = normrnd(0,1,10000, 1) + 1.5 * Dist(:, 1);
Covar = cov(Dist)

figure(1); clf
subplot(1, 2, 1)
plot(Dist(:, 1), Dist(:, 2), '.');
axis([-3 3 -3 3])
axis square
xlabel('x')
ylabel('y')

% Inspired by our previous scaling matrix, let's try the inverse of the
% square root of the covariance matrix:

Dist2 = Dist * inv(Covar^0.5);

figure(1)
subplot(1, 2, 1)
plot(Dist(:, 1), Dist(:, 2), '.');
axis([-3 3 -3 3])
axis square
xlabel('x')
ylabel('y')
subplot(1, 2, 2)
plot(Dist2(:, 1), Dist2(:, 2), '.');
axis([-3 3 -3 3])
axis square

% What the inverse of the square root of the covariance has done is rotated
% and scaled the axes.  Now in the x', y' coordinate system the
% distribution is circular (i.e. it has been spherized).

% The covariance generalizes to higher dimensions.  Let's consider a few
% more examples:

% Example 3: 6-d, no correlation between any components

clear Dist;
NumDim = 6
for dim = 1:NumDim
	Dist(:, dim) = normrnd(0,dim,10000, 1);			
end

cov(Dist)

% Again in this case the covariance matrix is diagonal, with elements equal to the 
% variance, and the off diagonal elements are not systematically different
% from 0 (run again to convince yourself of that).

% Example 4: 6-d, some elements correlated

clear Dist;
NumDim = 6
for dim = 1:NumDim
	Dist(:, dim) = normrnd(0,dim,10000, 1);			
	if (dim > 1)
		Dist(:, dim) = Dist(:, dim) + 0.5 * Dist(:, dim-1);			
	end
end

cov(Dist)

% Like the simple 2-d example, when we create a correlation between the variables
% corresponding to the different axes or dimensions of the data, we introduce
% non-zero off diagonal components.

% Keep this picture of the covariance matrix in mind when we turn to applications
% to problems in neuroscience below.

% Question set #1:
%	(1) Explain the value of the components of the 2 dimensional covariance matrix in 
%		Example 2 above.  
%	(2) Which of the off diagonal elements in the covariance matrix for Example 4 above
%		are nonzero?  Why?
%	(3) Can you separate inv(Covar^0.5) as above into the product of a rotation and 
%		a scaling matrix?


%-------------------------------------------------------------
% Part 3: Eigensystem of the covariance matrix
%
%	The eigenvectors of the covariance matrix turn out to provide a useful
% coordinate system. In particular, the eigenvalues determine how much of
% the variance of the distribution falls along the associated
% eigenvector.  Let's look at this for our 2d cases above:

% Example 1: x,y independent

clear Dist;
Dist(:, 1) = normrnd(0,0.3,10000, 1);			
Dist(:, 2) = normrnd(0,1,10000, 1);

C = cov(Dist);
[EigVec, EigVal] = eig(C);

% remember EigVec is a square matrix whose columns are the eigenvectors, and EigVal
% is a diagonal matrix with the eigenvalues along the diagonal.  Also remember the 
% eigenvectors v satisfy M v = l v, where l is the eigenvalue (a scalar).  Review this
% stuff in the linear algebra tutorial if it seems all new.

% superimpose eigenvectors on distribution
figure(1)
clf
subplot(1, 2, 1)
hold on
plot(Dist(:, 1), Dist(:, 2), '.');
axis([-3 3 -3 3])
axis square
xlabel('x')
ylabel('y')
PlotVector(EigVec(:, 1), 'r');
PlotVector(EigVec(:, 2), 'g');
hold off

% So the first eigenvector (red) points along the x-axis and the second (green)
% along the y-axis.  This makes sense since the covariance matrix is near diagonal.
% The eigenvectors should be close to (1 0) and (0 1).  Look at the associated eigenvalues:

EigVal

% They are equal to the variance associated with each eigenvector - i.e. the variance along 
% the x and y axes in this example.

% Example 2: x,y correlated

clear Dist;
Dist(:, 1) = normrnd(0,0.3,10000, 1);			
Dist(:, 2) = normrnd(0,1,10000, 1) + 1.5 * Dist(:, 1);

C = cov(Dist);
[EigVec, EigVal] = eig(C);

figure(1)
subplot(1, 2, 2)
hold on
plot(Dist(:, 1), Dist(:, 2), '.');
axis([-3 3 -3 3])
axis square
xlabel('x')
ylabel('y')
PlotVector(EigVec(:, 1), 'r');
PlotVector(EigVec(:, 2), 'g');
hold off

% now the eigenvectors no longer point along the x and y axis.  Instead,
% the eigenvector with the largest eigenvalue (green in this case) points
% in the direction along which the cloud of points is the most extended.
% The second eigenvalue points in an orthogonal direction.  Why do the
% eigenvectors have this property?  Remember we could think about the
% covariance matrix as a rotation and scaling matrix.  The rotation part
% would generate a new coordinate system x', y' in which the covariance was
% diagonal - like the independent case above.  If we view the eigenvectors
% in this coordinate system they should line up with the axes. In our
% original coordinate system these eigenvectors are at angles determined by
% the angle of the distribution.

% Again the amount of variance associated with each axis is given by the
% eigenvalue. Let's check this last statement.  We will project each point
% along each of the two eigenvectors and measure the associated variance:

% remember each row of Dist is an x-y pair (i.e. specifies the x and y
% components of the vector pointing to that data point).  Thus we can
% compute the projections along the eigenvectors by multiplying the Dist
% matrix by each eigenvector:

Proj1 = Dist * EigVec(:, 1);
Proj2 = Dist * EigVec(:, 2);

% now check the variance of each set of projections
var(Proj1)
var(Proj2)

% this should be equal to the associated eigenvalues

EigVal

% this again extends to higher dimensions:
clear Dist;
NumDim = 6
for dim = 1:NumDim
	Dist(:, dim) = normrnd(0,dim,10000, 1);			
	if (dim > 1)
		Dist(:, dim) = Dist(:, dim) + 0.5 * Dist(:, dim-1);			
	end
end
C = cov(Dist);
[EigVec, EigVal] = eig(C);
EigVal

% so in this case we see that the last eigenvector is associated with the
% most variance.  Again we can check that the associated eigenvalues are
% indeed equal to the variance:

for dim = 1:NumDim
	fprintf(1, 'Eigenvector %d: variance = %d\n', dim, var(Dist*EigVec(:, dim)));
end

% Why is all this useful?  Remember that the eigenvectors provide a
% coordinate system or a set of axes to represent the space spanned by the
% columns of a matrix.  But this is not an arbitrary coordinate system.
% The eigensystem is providing a coordinate system for our data where each
% axis is ranked in terms of how much variance is associated with it.  It
% is this ranking of the axes that provides a systematic way to reduce the
% number of dimensions describing some piece of data.  That is we can make
% a logical choice of which axes to keep and which to discard, and know how
% much of the variance we are capturing and failing to capture when we
% limit the number of dimensions.  

%-------------------------------------------------------------
% Part 4: Applications
%
% Principle components of single photon responses
%
%	This section illustrates how these tools can be used to provide a low
% dimensional representation of single photon responses from rod
% photoreceptors.  We know that the rod single photon responses vary from
% one to the next, and would like to have a compact way to characterize
% these response variations.  Let's see how we can do this using PCA.

% load in rod data - this should create two matrices: Singles and Failures.
% Each row of these matrices contains a single response (sampled every 10
% msec).  Singles contains single photon responses with the flash at time
% point 20, failures contains sections of dark record. Both singles and
% failures contain noise produced by the rod in darkness.  In addition, the
% single photon responses themselves vary.  Thus we want to identify and
% characterize the extra variability present in Singles compared to
% Failures.

load RodData

% plot mean and time-dependent variance of singles and failures

tme = 1:size(Singles, 2);
tme = (tme - 20) * 0.01;

figure(1);
clf;
subplot(1, 2, 1)
plot(tme, mean(Singles),tme, mean(Failures))
xlabel('time')
ylabel('current (pA)')
axis tight
subplot(1, 2, 2)
plot(tme, var(Singles),tme, var(Failures))
xlabel('time')
ylabel('variance (pA^2)')
axis tight

% Blue is the mean and variance of the singles, green of the failures.
% Look at the variance plot.  There is additional variance of the singles
% above that explained by the variance of the dark records.  This is the
% variance attributable to fluctuations in the single photon response.  The
% time-dependent variance, however, does not provide a complete description
% of how the responses vary because all it tells us is how much variance is
% present at one time point, but it fails to tell us how correlated the
% response is at two different times.  Another way of saying this is that
% specifying the mean and variance as in figure 1 does not allow us to
% generate the responses - we are missing information about how to go from
% the variance to a single response.  This is like telling you the square
% of a signal and asking you to estimate the signal itself; you can't
% because there are multiple signals that have the same square.  Likewise
% multiple distributions of single photon responses can have the same
% time-dependent variance.  PCA can help us here.  

% compute covariance matrix of singles corrected for covariance of failures
SinglesCovar = cov(Singles) - cov(Failures);

% look at eigenvalues
SinglesEigVal = eig(SinglesCovar);
figure(1)
clf
plot(SinglesEigVal, 'o')
xlabel('component')
ylabel('variance')

% So there are several large eigenvalues, and a bunch of small ones.  This looks promising: we
% probably want to keep the large ones, and hopefully can discard the small ones without losing too
% much.  We'll test this by comparing the variance of the singles and failures along each eigenvector
% and computing the variance of each

[EigVec, EigVal] = eig(SinglesCovar);
VarFailures = var(Failures * EigVec);
VarSingles = var(Singles * EigVec);

% plot variance of singles and failures projected along each eigenvector

figure(1)
plot([1:length(SinglesEigVal)], VarSingles, 'or', [1:length(SinglesEigVal)], VarFailures, '*')
xlabel('component')
ylabel('variance')

% Note that the variance of the projections along the last few eigenvectors
% is larger for the singles than failures.  This tells us that there is
% significant structure in the singles in this direction (i.e. along this
% axis) that is not present in the failures.  If the variance of the
% singles and failures is similar we reject that axis because it is not
% picking up structure beyond that expected from the dark fluctuations in
% current. Lets focus in on the region of interest

figure(1)
xlim([110, 120]);

% It looks like most of the action is in the last 2 or 3 components.  Let's
% look at those. We will reduce each response to a single point in a 3-d
% space with axes given by the last 3 eigenvectors.  Each response is
% specified by its projections along these axes. First lets use 'eigs' to
% return the 3 eigenvectors associated with the largest eigenvalues - these
% are the 3 axes we are interested in:

% get a few selected eigenvectors (by default eigs gives us those with the
% largest eigenvalues):

[EigVec, EigVal] = eigs(SinglesCovar, 3);

% look at these
figure(1)
clf
plot(EigVec);
xlabel('time')

% they should look reasonable.  In particular they should not have much 
% structure in the first 20 points (which are prior to the flash). Note,
% however, that none of these looks like the mean single photon response. 

% Now we will describe each individual response by its projection along
% these responses - in other words by the weight of each component.  First
% lets look at a couple real responses and their descriptions in terms of
% the eigenvectors:

figure(2)
clf
start = 10;				% first response to look at
for n = 1:3
	resp(1:size(Singles, 2)) = 0;
	fprintf(1, '\nResponse %d ', n+start);
	% add up weighted versions of each eigenvector
	for comp = 1:3
		resp = resp + Singles(n+start, :) * EigVec(:, comp) * EigVec(:, comp)';
		fprintf(1, 'coefficient %d = %d ', comp, Singles(n+start, :) * EigVec(:, comp));
	end
	subplot(1, 3, n)
	plot([1:size(Singles, 2)], Singles(n+start, :), [1:size(Singles, 2)], resp)
end

% So we are expanding each response in a series with the basis set given by
% the eigenvectors and the expansion coefficients given by the projection
% of the response along the eigenvector.  The blue lines are the original
% responses, the green traces the fits using the projections along each of
% the three eigenvectors.  Remember there is also dark noise in the
% recorded responses, so we don't expect the agreement to be perfect.

% Now lets look at the distributions of singles and failures along the
% three axes:

figure(2)
clf
plot3(Singles * EigVec(:, 1), Singles * EigVec(:, 2), Singles * EigVec(:, 3), 'or')
hold on
plot3(Failures * EigVec(:, 1), Failures * EigVec(:, 2), Failures * EigVec(:, 3), '*')
xlabel('component 1')
ylabel('component 2')
zlabel('component 3')
grid on
rotate3d on 

% Rotate the plot around (using the mouse) to get a sense of how the
% structures of the two distributions differ.  You should see two things.
% First, there is a difference between the means of the two distributions.
% That's good - because there is a lot of structure in the mean single
% photon response that should get captured.  Second, the spread of the
% singles distribution looks larger than that of the failures, clearly
% along the component 1 and 2 axes, and probably along the component 3
% axis.  

% What does this do for us?  We have described the set of single photon
% responses in a relatively low dimensional space.  This helps us in a
% couple ways.  First, we can describe the nature of the variations in the
% responses about the mean, and compare these to various models for how
% these variations are generated.  Second, we can use this procedure to
% make a generative model that produces responses with characteristics
% close to those of the actual single photon responses.  

% Question set #2:
%	(1) Can you think of another case where dimensional reduction is important?
%	(2) How much of the variance of the singles are we capturing with the first
%	    3 components above?  How much more do we get with 10 components?
%   (3) Why do we subtract the covariance of the failures from that of
%       the singles? 
%   (4) Why don't any of the eigenvectors we recover look like the mean single
%       photon response?