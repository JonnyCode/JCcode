% Example script for computing the optical flow for an image sequence
% It includes visualisation of subsampled flow fields and
% computation of the error statistics


clear all

		%%%%%%%%%%%%%%
		% Parameters %
		%%%%%%%%%%%%%%
gx = 25;        % Number of evaluation points in X dimension
thres_lin = .01;	% Linearity threshold for phase gradient
nc_min = 7;		% Minimum number of valid component velocities

		%%%%%%%%%%%%%%%%%
		% Read Sequence %
		%%%%%%%%%%%%%%%%%
II = read_sequence ('yosemite/yosemite256.', 9, 5);
[sy sx st] = size(II);


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Compute Optical Flow Field %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O = optical_flow (II, gx, thres_lin, nc_min);

		%%%%%%%%%%%%%%%%%%
		% The Art Corner %
		%%%%%%%%%%%%%%%%%%
mag = 1;
pcolor (II(:,:,floor(st/2)+1));
shading flat
colormap gray
axis image
hold on
offset = ceil(9.31648319*sqrt(log(100)));   % Point where Gaussian envelope drops to 10%
[Vx Vy] = vis_flow (O(:,:,1), O(:,:,2), gx, offset, 1, 'm');
hold off



		%%%%%%%%%%%%
		% Evaluate %
		%%%%%%%%%%%%
% Coverage (a bit messy)
if (gx==0)
	jmp = 1;
else
	jmp = floor(sx/gx);
	jmp = jmp + (jmp==0);
end
nv = length(offset+1:jmp:(sy-offset))*length(offset+1:jmp:(sx-offset));
fprintf ('Coverage: %.2f %%\n', sum(sum(~isnan(Vx)))/nv*100);

% Error statistics
[Cx Cy] = read_flow ('yosemite/correct_yosemite256');
[phif,stf,phic,stc,Ef,Ecfin] = eval_flow (O(:,:,1), O(:,:,2), Cx, Cy);
fprintf ('Angular error (full): %.4f deg (%.4f deg)\n', phif, stf)


