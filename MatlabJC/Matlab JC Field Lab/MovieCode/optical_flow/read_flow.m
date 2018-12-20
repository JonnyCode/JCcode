% Function to read the optical flow field (Barron's format)
%
%
% Usage: [U, V] = read_flow (name_in);

function [U, V] = read_flow (name_in)

		%%%%%%%%%%%%%
		% Load Flow %
		%%%%%%%%%%%%%
fid = fopen(name_in, 'r', 'ieee-be');
totx = fread (fid, 1, 'float');
toty = fread (fid, 1, 'float');
numx = fread (fid, 1, 'float');
numy = fread (fid, 1, 'float');
offx = fread (fid, 1, 'float');
offy = fread (fid, 1, 'float');
[X count] = fread (fid,[2*numx numy], 'float');
fclose (fid);
X(X==100) = NaN;
%disp ([totx toty numx numy offx offy]);

		%%%%%%%%%%%
		% Reshape %
		%%%%%%%%%%%
T = permute (reshape(X,[2,numx,numy]),[3 2 1]);
T = T(numy:-1:1,:,:);

		%%%%%%%%%%
		% Offset %
		%%%%%%%%%%
U = NaN.*ones(toty,totx);
V = NaN.*ones(toty,totx);
INDx = (offx+1):(offx+numx);
INDy = (offy+1):(offy+numy);
U(INDy,INDx) = T(:,:,1);
V(INDy,INDx) = T(:,:,2);
