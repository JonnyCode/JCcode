% TLS tutorial
% It's often the case that we'd like to see how a linear combination of
% inputs acounts for some outcome.  A typical situations is one in which
% we have n-1 inputs that we're manipulating, and m different
% observations of those variables.  We also observe m different outcomes,
% one for each set of input values.  Alternatively we can view the input and outcomes
% as n variables.

% This tutorial deals with the problem of finding a set of linear
% coefficients that minimize the error in the output, given the inputs.

% The simplest case is one in which we have a single input (x), or a set of
% x_i observations of that input, and what we hypothesize is that the outcome (y), or
% a set of y_i observations, depends linearly on the input.  Then, we're looking for a, b such that
% ax_i+b=y minimizes the total mean-squared error between the actual
% observation, y_i, and the linear prediction from {x_i}; we're looking for the
% best fit line.

% We can extend this situation to any number of variables (n) and any
% number of observations (m).

% This operation is also related to PCA.

%% 1) 2D case, ie, 1 independent variable, 1 supposed dependent variable
clear

% First we generate some data.
% For clarity, we can start with >3 observations, then increase the number.
% We can do this procedure for systems that are under defined, that is we
% have fewer observations than variables, but then the solution we'll find
% isn't unique.  We'll stick to the case in which the system has a
% solution, so at least two observations for 2 variables.  This case isn't
% all that interesting though: we know the best fit line will be the line
% connecting the points.  So we'll start with at least three
% observations. More points make it easier to visualize though.  

numobservations = 100;
xnoise = 1.*randn(numobservations,1);
ynoise = 1.*randn(numobservations,1);

% We'll now build a correlation between the variables, you can change the
% ammount of correlation, or even whether it's a linear correlation.
x = 1*xnoise;
y = .4*ynoise+1.5*xnoise.^2;

close all
figure(1);
ax = gca;
hold on
plot(x,y,'.','Color',[1,.7 .7])
title('Observations of x and y are correlated')


%% The matrix M (m x n) describes the system
% We can make a matrix of all of the observations: [x y].  We have m
% observations, so we have m rows.  Also, n is 2.

% it's easiest to substract off the mean of each observation.  That way,
% each observation is actually a deviation from the mean, and the origin is
% in the centroid of the cloud of observations.  We can add the origin
% back to our solution later, but this reduces the 2D problem to finding
% the single coeficient (a) that relates x to y with the least total error
origin = [mean(x), mean(y)];
M = [x-mean(x),y-mean(y)];

plot(M(:,1),M(:,2),'.','Color',[1 .5 .5])
pbaspect([1 1 1]);
daspect([1 1 1]);
set(ax,'Xgrid','on');
set(ax,'Ygrid','on');
plot(origin(1),origin(2),'+','Color',[1,.7 .7])
plot(0,0,'+','Color',[1 0 0])
plot([origin(1),0],[origin(2),0],'--','Color',[1 .7 .7]);
% set(ax,'XGrid','on','YGrid','on');

title('Move the origin to the mean of input (x) and output (y) values')

%% SVD
% Now we come to the magic.  We want to orthogonalize the cloud of
% observations.  We want to rotate the cloud such that most of the
% variation lies along one axis, with the orthogonal axis having less
% variation.  We can see this well when there are a lot of observations,
% that there will be what looks like a skinny axis and a fat axis.  Our
% hypothesis says that the long axis is the correct relationship, and the
% skinny access results from error.  
% So how do we perfom this rotation? We're going to go about it in a
% slightly round about way.  The first step is to realize that every matrix M can
% be represented by the multiplication of three matrices U, S and V', USV',
% where these three matrices have interesting properties.

% Perform SVD
[U,S,V] = svd(M);

%% Outputs are:
% U (m x m), column vectors form orthonormal basis set with at least m-n
% vectors describing the null space of M in R^m.
% S (m x n)is a diagonal matrix of singular values for M. That is, the diagonal
% elemenents satisfy:
%   Mv = s u, where v is a column vector of V, s is a diagonal value of S
%   and u is a column vector of U;
% In the case of linear regression, when we have more observations than
% variables, there will be at most n singular values and n colums of S.
% V (n x n) is a set of orthonormal basis vectors, ie a unitary matrix.

% illustration of properties of U, S, and V
sum(diag(U*U')) % Unitary, m x m
sum(diag(V*V')) % unitary, n x n
diag(S)         % singular values


%% We can also think of U as orthogonal vectors in R^m (mD space) that describe all the
% observations and S is the matrix that projects those observation into n
% space.  S also scales the axes.  See what happens when we don't scale
% them: ie S(S>0) = 1;
L = S;
L(L>0) = 1;
UinNspace = U*L;
plot(UinNspace(:,1),UinNspace(:,2),'.','color',[.8 .8 1]);
title('Light blue cloud is the unscaled projection of U into 2D space');

%% U*S is n column vectors in R^m that represent a scaled version of the R^m
% basis vectors, ie stretching or squeezing along the right singular
% vectors (most of which go to zero, ie there are n degrees of freedom in
% R^m, and the rest of the directions are 0).

% U*S can now be thought of as the projections of all the observations onto
% the orthonormal basis vectors of V.  V is then the rotation of the
% coordinates of U*S into the cartesian plane.

obsInVBasis = U*S; 
obsInCBasis = obsInVBasis*V';  % rotation of observations into x,y coordinates

plot(obsInVBasis(:,1),obsInVBasis(:,2),'.','color',[.5 .5 1]);
axis tight
pbaspect([1 1 1]);
daspect([1 1 1]);
title('U*S represents observations in an appropriate basis space, that spanned by the columns of V')


%% display both basis vectors in V, and those vectors represented in
% the cartesian basis of x and y as well

compass(1,0,'b');
compass(0, 1,'b');
compass(V(1,1),V(2,1),'r')
compass(V(1,2),V(2,2),'r')
axis tight
pbaspect([1 1 1]);
daspect([1 1 1]);

title('V is the operation that rotate the x,y unit vectors')

%% try a variety of these:
% change ii, the row index, to see how an observation in U*S moves to an
% observation in M when the iith row of U*S acts on V to produce a new
% vector in the cartesian plane.  Compare the blue coordinates to the red
% coordinates
%
%  NOTE: axes might be both rotated and FLIPPED

ii = 14;
obInVBasis = obsInVBasis(ii,:);

% go along the first coordinate in V
plot([0,obInVBasis(1)],[0,0]);

% go along the second coordinate in V
plot([obInVBasis(1),obInVBasis(1)],[0,obInVBasis(2)]);

% go along the first coordinate in V rotated into cartesian plane
plot([0,obInVBasis(1)*V(1,1)],[0,obInVBasis(1)*V(2,1)],'r');
% go along the second coordinate in V rotated into cartesian plane
plot([obInVBasis(1)*V(1,1),obInVBasis(1)*V(1,1)+obInVBasis(2)*V(1,2)],...
    [obInVBasis(1)*V(2,1),obInVBasis(1)*V(2,1)+obInVBasis(2)*V(2,2)],'r');

title(sprintf('Row %g of U*S (blue) rotates to row %g (red) of M by application of V''',ii,ii))

%% Least squares solution
% We can see that we would like the error to lie along the dimension that
% scales U the least.
% So, we estimate that the diagonal value of S that scales U the least is
% ACTUALLY 0. 

fprintf('Diagonal values of S: \n');
fprintf('%g ',diag(S));

[smallsigma,column] = min(diag(S));
Soptimal = S;
Soptimal(column,column) = 0;

fprintf('\nNew Diagonal values of S:\n')
fprintf(' %g ',diag(Soptimal));
fprintf('\n')

%% Now we can compress U into R^n in the same way as before, multiply U*S.
% Now, any component of U*S that strecheted along the vector associated
% with the smallest singular value is 0.
obsInOptimalBasis = U*Soptimal;
plot(obsInOptimalBasis(:,1),obsInOptimalBasis(:,2),'.','LineStyle','-','Color',[.5 1 .5]);

title('Projections onto basis vector with smallest singular value are 0')
%%  Rotate the observations into the cartesian coordinates as before:
Moptimal = obsInOptimalBasis*V';
plot(Moptimal(:,1),Moptimal(:,2),'.','LineStyle','-','Color',[1 1 0]);

title('Rotate space using V''');
%% Best fit linear approximation
% Now it appears that the only thing left is to find the equation of the
% best fit line
% Lets take a look at the basis vector in V that we effectively threw out
% by making its singular value 0:
column

V(:,column) 

% define V = [v1, v2] the column vectors of V.  ie V(:,2) = v2

% We know that M*V = U*S;  When we zeroed the lowest singular value, we
% also made the second column of U*S entirely 0.  Check it out:
sum(obsInOptimalBasis(:,2).*obsInOptimalBasis(:,2))

% This is equivalent to saying that M * v2 = 0
% This means that for all x_i and y_i,  x*V(1,2) + y*V(2,2) = 0, which is
% the equation of a line mx = y:
x1 = min(x)-mean(x);
x2 = max(x)-mean(x);
y1 = min(y)-mean(y);
y2 = max(y)-mean(y);
plot([x1,x2],-V(1,2)*[x1,x2]/V(2,2));

title(sprintf('Column %g of V gives the equation for the best fit line',column));

compass(V(1,2),V(2,2),'b')

% this seems counterintuitive, that you'd take the vector that needed to be
% zeroed as the best fit line, but that's precisely the point.  What we've
% looked for is a direction in n-space, in this case 2D, along which any
% projection we see, we'll call an error, and say that the projection is
% actually 0;  So: we project any x and y values onto that vector, 
% x*v2(1) + y*v2(2) (dot product)
% and make it 0 which again gives us the equation for a line;

% Now try more observations by increasing numobservations in the first cell

%% 2) 3D case, ie, 2 independent variables, 1 supposed dependent variable
clear

numobservations = 100;
xnoise = 1.*randn(numobservations,1);
ynoise = 1.*randn(numobservations,1);
znoise = 1.*randn(numobservations,1);

% We'll now build a correlation between the variables, you can change the
% ammount of correlation, or even whether it's a linear correlation.
x = 1*xnoise;
y = ynoise+1.5*xnoise;
z = .6*x + .79*y + 3*znoise;

close all
figure(1);
ax = gca;
plot3(x,y,z,'.','Color',[1,.7 .7])
hold on
title('Observations of x and y are correlated')
xlabel('x');
ylabel('y');
zlabel('z');


%% The matrix M (m x n) describes the system
% We can make a matrix of all of the observations: [x y].  We have m
% observations, so we have m rows.  Also, n is 2.

% it's easiest to substract off the mean of each observation.  That way,
% each observation is actually a deviation from the mean, and the origin is
% in the centroid of the cloud of observations.  We can add the origin
% back to our solution later, but this reduces the 2D problem to finding
% the single coeficient (a) that relates x to y with the least total error
origin = [mean(x), mean(y), mean(z)];
M = [x-mean(x),y-mean(y),z-mean(z)];

plot3(M(:,1),M(:,2),M(:,3),'.','Color',[1 .5 .5])
pbaspect([1 1 1]);
daspect([1 1 1]);
set(ax,'Xgrid','on');
set(ax,'Ygrid','on');
set(ax,'Zgrid','on');
plot3(origin(1),origin(2),origin(3),'+','Color',[1,.7 .7])
plot3(0,0,0,'+','Color',[1 0 0])
plot3([origin(1),0],[origin(2),0],[origin(3),0],'--','Color',[1 .7 .7]);
% set(ax,'XGrid','on','YGrid','on');

title('Move the origin to the mean of input (x) and output (y) values')

%% SVD
% Now we come to the magic.  We want to orthogonalize the cloud of
% observations.  We want to rotate the cloud such that most of the
% variation lies along one axis, with the orthogonal axis having less
% variation.  We can see this well when there are a lot of observations,
% that there will be what looks like a skinny axis and a fat axis.  Our
% hypothesis says that the long axis is the correct relationship, and the
% skinny access results from error.  
% So how do we perfom this rotation? We're going to go about it in a
% slightly round about way.  The first step is to realize that every matrix M can
% be represented by the multiplication of three matrices U, S and V', USV',
% where these three matrices have interesting properties.

% Perform SVD
[U,S,V] = svd(M);


%% Outputs are:
% U (m x m), column vectors form orthonormal basis set with at least m-n
% vectors describing the null space of M in R^m.
% S (m x n)is a diagonal matrix of singular values for M. That is, the diagonal
% elemenents satisfy:
%   Mv = s u, where v is a column vector of V, s is a diagonal value of S
%   and u is a column vector of U;
% In the case of linear regression, when we have more observations than
% variables, there will be at most n singular values and n colums of S.
% V (n x n) is a set of orthonormal basis vectors, ie a unitary matrix.

% illustration of properties of U, S, and V
sum(diag(U*U')) % Unitary, m x m
sum(diag(V*V')) % unitary, n x n
diag(S)         % singular values

%% We can also think of U as orthogonal vectors in R^m (mD space) that describe all the
% observations and S is the matrix that projects those observation into n
% space.  S also scales the axes.  See what happens when we don't scale
% them: ie S(S>0) = 1;
L = S;
L(L>0) = 1;
UinNspace = U*L;
plot3(UinNspace(:,1),UinNspace(:,2),UinNspace(:,3),'.','color',[.8 .8 1]);
title('Light blue cloud is the unscaled projection of U into 2D space');

%% U*S is n column vectors in R^m that represent a scaled version of the R^m
% basis vectors, ie stretching or squeezing along the right singular
% vectors (most of which go to zero, ie there are n degrees of freedom in
% R^m, and the rest of the directions are 0).

% U*S can now be thought of as the projections of all the observations onto
% the orthonormal basis vectors of V.  V is then the rotation of the
% coordinates of U*S into the cartesian plane.

obsInVBasis = U*S; 
obsInCBasis = obsInVBasis*V';  % rotation of observations into x,y coordinates

plot3(obsInVBasis(:,1),obsInVBasis(:,2),obsInVBasis(:,3),'.','color',[.5 .5 1]);
axis tight
pbaspect([1 1 1]);
daspect([1 1 1]);
title('U*S represents observations in an appropriate basis space, that spanned by the columns of V')

%% display both basis vectors in V, and those vectors represented in
% the cartesian basis of x and y as well

quiver3(0,0,0,1,0,0,'b');
quiver3(0,0,0,0,1,0,'b');
quiver3(0,0,0,0,0,1,'b');
quiver3(0,0,0,V(1,1),V(2,1),V(3,1),'r');
quiver3(0,0,0,V(1,2),V(2,2),V(3,2),'r');
quiver3(0,0,0,V(1,3),V(2,3),V(3,3),'r');
axis tight
pbaspect([1 1 1]);
daspect([1 1 1]);

title('V is the operation that rotate the x,y unit vectors')

%% try a variety of these:
% change ii, the row index, to see how an observation in U*S moves to an
% observation in M when the iith row of U*S acts on V to produce a new
% vector in the cartesian plane.  Compare the blue coordinates to the red
% coordinates
%
%  NOTE: axes might be both rotated and FLIPPED

ii = 14;
obInVBasis = obsInVBasis(ii,:);

% go along the first coordinate in V
plot3([0,obInVBasis(1)],[0,0],[0,0]);

% go along the second coordinate in V
plot3([obInVBasis(1),obInVBasis(1)],[0,obInVBasis(2)],[0,0]);
% go along the third coordinate in V
plot3([obInVBasis(1),obInVBasis(1)],[obInVBasis(2),obInVBasis(2)],[0,obInVBasis(3)]);

% go along the first coordinate in V rotated into cartesian plane
plot3([0,obInVBasis(1)*V(1,1)],[0,obInVBasis(1)*V(2,1)],[0,obInVBasis(1)*V(3,1)],'r');
% go along the second coordinate in V rotated into cartesian plane
plot3([obInVBasis(1)*V(1,1),obInVBasis(1)*V(1,1)+obInVBasis(2)*V(1,2)],...
    [obInVBasis(1)*V(2,1),obInVBasis(1)*V(2,1)+obInVBasis(2)*V(2,2)],...
    [obInVBasis(1)*V(3,1),obInVBasis(1)*V(3,1)+obInVBasis(2)*V(3,2)],'r');
% go along the third coordinate in V rotated into cartesian plane
plot3([obInVBasis(1)*V(1,1)+obInVBasis(2)*V(1,2),obInVBasis(1)*V(1,1)+obInVBasis(2)*V(1,2)+obInVBasis(3)*V(1,3)],...
    [obInVBasis(1)*V(2,1)+obInVBasis(2)*V(2,2),obInVBasis(1)*V(2,1)+obInVBasis(2)*V(2,2)+obInVBasis(3)*V(2,3)],...
    [obInVBasis(1)*V(3,1)+obInVBasis(2)*V(3,2),obInVBasis(1)*V(3,1)+obInVBasis(2)*V(3,2)+obInVBasis(3)*V(3,3)],'r');

title(sprintf('Row %g of U*S (blue) rotates to row %g (red) of M by application of V''',ii,ii))

%% Least squares solution
% We can see that we would like the error to lie along the dimension that
% scales U the least.
% So, we estimate that the diagonal value of S that scales U the least is
% ACTUALLY 0. 

fprintf('Diagonal values of S: \n');
fprintf('%g ',diag(S));

[smallsigma,column] = min(diag(S));
Soptimal = S;
Soptimal(column,column) = 0;

fprintf('\nNew Diagonal values of S:\n')
fprintf(' %g ',diag(Soptimal));
fprintf('\n')

%% Now we can compress U into R^n in the same way as before, multiply U*S.
% Now, any component of U*S that strecheted along the vector associated
% with the smallest singular value is 0.
obsInOptimalBasis = U*Soptimal;
plot3(obsInOptimalBasis(:,1),obsInOptimalBasis(:,2),obsInOptimalBasis(:,3),'.','LineStyle','none','Color',[.5 1 .5]);

%% Projection of octants
octants = [-1,-1,-1
    -1,-1,1
    -1,1,-1
    -1,1,1
    1,1,1
    1,1,-1
    1,-1,1
    1,-1,-1];
   
Uhat = [1 0 0; 0 1 0; 0 0 1];
UhatOptimal = Uhat;
UhatOptimal(:,column) = 0;

planeInOptimalBasis = octants*UhatOptimal;
planeInOptimalBasis = uniquecoordinates(planeInOptimalBasis);
planeInOptimalBasis = shortestPath(planeInOptimalBasis);

h = patch(planeInOptimalBasis(:,1),planeInOptimalBasis(:,2),planeInOptimalBasis(:,3),'b');
set(h,'EdgeAlpha',0.15,'EdgeColor',[.5 1 .5],'FaceAlpha',0.15,'FaceColor',[.5 1 .5])

title('Projections onto basis vector with smallest singular value are 0')
%%  Rotate the observations into the cartesian coordinates as before:
Moptimal = obsInOptimalBasis*V';
plot3(Moptimal(:,1),Moptimal(:,2),Moptimal(:,3),'.','LineStyle','none','Color',[1 1 0]);

title('Rotate space using V''');
%% Best fit linear approximation
% Now it appears that the only thing left is to find the equation of the
% best fit line
% Lets take a look at the basis vector in V that we effectively threw out
% by making its singular value 0:
column

V(:,column) 

% define V = [v1, v2, v3] the column vectors of V.  ie V(:,3) = v3

% We know that M*V = U*S;  When we zeroed the lowest singular value, we
% also made the third column of U*S entirely 0.  Check it out:
sum(obsInOptimalBasis(:,column).*obsInOptimalBasis(:,column))

% This is equivalent to saying that M * v3 = 0
% This means that for all x_i, y_i, and z_i,  x*V(1,3) + y*V(2,3) + z*V(3,3) = 0, which is
% the equation of a plane:

planeInCBasis = planeInOptimalBasis*V';

h = patch(planeInCBasis(:,1),planeInCBasis(:,2),planeInCBasis(:,3),'b');
set(h,'EdgeAlpha',0.15,'EdgeColor',[1 1 0],'FaceAlpha',0.2,'FaceColor',[1 1 0])

title('Projections onto basis vector with smallest singular value are 0')

quiver3(0,0,0,V(1,3),V(2,3),V(3,3),'k','Linewidth',2)

% Finally, the coefficients that make the most sense are:
-V(1,3)/V(3,3)
-V(2,3)/V(3,3)

% this seems counterintuitive, that you'd take the vector that needed to be
% zeroed as the best fit line, but that's precisely the point.  What we've
% looked for is a direction in n-space, in this case 2D, along which any
% projection we see, we'll call an error, and say that the projection is
% actually 0;  So: we project any x and y values onto that vector, 
% x*v2(1) + y*v2(2) (dot product)
% and make it 0 which again gives us the equation for a line;

% Now try more observations by increasing numobservations in the first cell

%% Extension to PCA
% This problem looks a lot like Principle Component Analysis.  In PCA, we
% try to find a basis set that orthogonalizes the variation and then to
% find a subspace that capture the most variation. SVD will give us the
% directions that orthoganalize the variation, that is, find the set of
% basis vectors that capture the variance.  in PCA, we do this by
% calculating the covariation matrix and then finding the orthoganal
% eigenvectors (orthagonal because the covariance matrix is symetric).

% thus, our basis vectors V must be the eigenvectors of the covariance
% matrix!  We just have to see how that is:

% We have just seen that M = U*S*V'.  Then, by multiplication with matrices.
%
% U'*M = U'*U*S*V';
% S'*U'*M = S'*U'*U*S*V';
% V*S'*U'*M = V'*S'*U'*U*S*V';
% M'*M = V*S'*U'*U*S*V';

% Since U is unitary, U'*U is the identity matrix;
sum(diag(U'*U))

% So, M'*M = V*S'*S*V';
% S'*S is a square matrix with only diagonal elements, which looks like this

varMat = S'*S;

diag(varMat)

% We can also recognize M'*M as the covariance matrix!  We'll call this
% Mcov, such that M'*M = Mcov = V *varMat*V';

% We can then write, Mcov * V = V * varMat *V'*V = V*varMat  (V is unitary)
Mcov = M'*M;

% Finally, let's look at the operation V * varMat.  The resulting columns
% are the columns of V multiplied by the diagonal element in the
% corresponding column of varMat.  This means that for each column of V:
fprintf('Compare subsequent vector pairs\n');
disp(Mcov * V(:,1))
disp(varMat(1,1)*V(:,1))

disp(Mcov * V(:,2))
disp(varMat(2,2)*V(:,2))

disp(Mcov * V(:,3))
disp(varMat(3,3)*V(:,3))

% This is the definition of eigenvectors.  This means that by finding V, we
% have the orthonormal basis of eigen vectors of the covariance matrix,
% with variance along those vectors as the square of the diagonal elements
% of S.

%% Goodness of fit
%  We have thought of the last column in M as the dependent variable.  We
%  now have come up with a model that cuts our original n space down to
%  n-1 space, treating the distances to this n-1 space along the direction
%  perpendicular to the plane, as errors.   In the final cell, the regress
%  function, I'm not sure I understand the operation it's doing.  
% essentially I believe it models all of the error along 1 dimension, the
% z-dimension, and finds the plane that minimizes that error.   


IVs = M(:,1:end-1);
DV = M(:,end);

b = regress(DV,IVs);

planeInRBasis = planeInOptimalBasis;
planeInRBasis(:,3) = b(1)*planeInRBasis(:,1)+b(2)*planeInRBasis(:,2);

h = patch(planeInRBasis(:,1),planeInRBasis(:,2),planeInRBasis(:,3),'b');
set(h,'EdgeAlpha',0.15,'EdgeColor',[1 1 1]*.6,'FaceAlpha',0.2,'FaceColor',[1 1 1]*.6)

dvfit = b(1)*IVs(:,1)+b(2)*IVs(:,2);
plot3(IVs(:,1),IVs(:,2),dvfit,'k.');





