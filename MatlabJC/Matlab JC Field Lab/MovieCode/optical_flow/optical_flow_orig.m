% Phase-based Opic Flow Algorithm, described in
% Gautama, T. and Van Hulle, M.M. (2002).  A Phase-based Approach to the
% Estimation of the Optical Flow Field Using Spatial Filtering.
% IEEE Trans. Neural Networks, 13(5), 1127--1136.
%
% Usage: O = optical_flow (II, gx, thres_lin, nc_min)
%	II [sy sx st]	Image Sequence (Y-X-t)
%	gx		Number of velocity vectors along X-axis (0=all)
%	thres_lin	Linearity threshold [.05]
%	nc_min		Minimal number of valid component velocities for
%				computation of full velocity [5]

function O = optical_flow (II, gx, thres_lin, nc_min)

if (nargin<1)
	error ('Please provide an input sequence');
end
if (nargin<2)
	gx = 0;
end
if (nargin<3)
	thres_lin = .05;
end
if (nargin<4)
	nc_min = 5;
end
[sy sx st] = size(II);


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Goal Programming Network Parameters %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ep = 500;
dt = .001;


if (gx==0)
	jmp = 1;
else
	jmp = floor(sx/gx);
	jmp = jmp + (jmp==0);
end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Load Filterbank Parameters %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = [	0.02156825 -0.08049382; ...
	0.05892557 -0.05892557; ...
	0.08049382 -0.02156825; ...
	0.08049382 0.02156825; ...
	0.05892557 0.05892557; ...
	0.02156825 0.08049382; ...
	0.06315486 -0.10938742; ...
	0.10938742 -0.06315486; ...
	0.12630971 0.00000000; ...
	0.10938742 0.06315486; ...
	0.06315486 0.10938742];
S = [9.31648319 9.31648319 9.31648319 9.31648319 9.31648319 9.31648319 ...
	6.14658664 6.14658664 6.14658664 6.14658664 6.14658664]';
nn = size(W,1);

		%%%%%%%
		% Aux %
		%%%%%%%
xx = (1:st);
xx3 = zeros(1,1,st);
xx3(1:st) = 1:st;
Sxx = sum(xx.^2);
Sx = sum(xx);
den = (st.*Sxx-Sx.^2);
pi2 = 2*pi;


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Compute Filter Outputs & Component Velocities %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tclock1 = clock;
AC = zeros(sy,sx,st);
AS = zeros(sy,sx,st);
FV = zeros(nn,sy,sx);	% Filter Component Velocity
LE = zeros(nn,sy,sx);	% MSE of Regression
Ec = zeros(sy,sx,nn);
offset = ceil(max(S)*sqrt(log(100)));	% Point where Gaussian envelope drops to 10%
%offset=0;
[offx offy] = meshgrid(1:sx,1:sy);
for n=1:nn,fprintf('+');end;
for i=offset+1:jmp:(sy-offset),fprintf('.');end;fprintf('\n');
for n=1:nn
	% Generate 1D kernels
	win_len = floor(6*S(n));
	cx = (0:win_len-1)'-win_len/2+.5;
	G = exp(-cx.^2./(2*S(n)*S(n)))./(sqrt(2*pi)*S(n));
	FCx = G.*cos((2*pi*W(n,1)).*cx);
	FCy = G.*cos((2*pi*W(n,2)).*cx);
	FSx = G.*sin((2*pi*W(n,1)).*cx);
	FSy = G.*sin((2*pi*W(n,2)).*cx);

	% Exceptions for null frequencies
	if (sum(FCx.^2)==0)
		FCx = ones(win_len,1);
	end
	if (sum(FCy.^2)==0)
		FCy = ones(win_len,1);
	end
	if (sum(FSx.^2)==0)
		FSx = ones(win_len,1);
	end
	if (sum(FSy.^2)==0)
		FSy = ones(win_len,1);
	end

	% Perform Convolutions, room for improvement (subsampling)
	for t=1:st
		IIp = II(:,:,t)';

		% Sine Filter
		Tsx = conv2(IIp, FSx, 'same')';
		T2 = conv2(Tsx, FCy, 'same');
		Tcx = conv2(IIp, FCx, 'same')';
		T4 = conv2(Tcx, FSy, 'same');
		AS(:,:,t) = T2 + T4;

		% Cosine Filter
		%Tcx = conv2(IIp, FCx, 'same')';
		T2 = conv2(Tcx, FCy, 'same');
		%Tsx = conv2(IIp, FSx, 'same')';
		T4 = conv2(Tsx, FSy, 'same');
		AC(:,:,t) = T2 - T4;
	end

	% Compute and Unwrap Phase
	Mcos = (AC==0);
	P = atan(AS./(AC+Mcos))+pi.*(AC<0);
	P(Mcos) = NaN;
	k = 2;
	while (k<=st)
		D = P(:,:,k) - P(:,:,k-1);
		A = abs(D)>pi;
		P(:,:,k:st) = P(:,:,k:st) - repmat(pi2.*sign(D).*A,[1 1 st-k+1]);
		k = k + (sum(sum(A))==0);
	end

	% Compute Filter Component Velocity
	Sxy = sum(repmat(xx3,[sy sx 1]).*P,3);
	Sy = sum(P,3);
	a = (Sxx.*Sy-Sx.*Sxy)./den;
	b = (st.*Sxy-Sx.*Sy)./den;
	Reg = repmat(a,[1 1 st])+repmat(b,[1 1 st]).*repmat(xx3,[sy sx 1]);
	LE(n,:,:) = mean((Reg-P).^2,3)./abs(b+(b==0));
	FV(n,:,:) = -b./(pi2*sum(W(n,:).^2)).*(W(n,1)+sqrt(-1)*W(n,2));

	fprintf ('*');

end
tclock2 = clock;
time1 = etime(tclock2,tclock1);

		%%%%%%%%%%%%%%%%%%%%%%%%%
		% Compute Full Velocity %
		%%%%%%%%%%%%%%%%%%%%%%%%%
O = NaN.*ones(sy,sx,2);
E = NaN.*ones(sy,sx);
Egpn = NaN.*ones(sy,sx);
for i=offset+1:jmp:(sy-offset)
  for j=offset+1:jmp:(sx-offset)

	% Linearity Check
	IND1 = find(LE(:,i,j)<thres_lin);
	V = FV(IND1,i,j);
	nc = length(IND1);

	if (nc>=nc_min)
		% GPN
		B = abs(V);
		D = [real(V) imag(V)]./repmat(B+(B==0),1,2);
		[R,F,E,U,LR,conv_flag] = gpn_ada (D, B, ones(nc,1), ep, dt, [0 0], 1.0);
		if (sum(R.^2)~=0)
			O(i,j,:) = R;
			Egpn(i,j) = min([NaN;E(E~=0)]);
		end
	end
  end
  fprintf ('*');
end
fprintf ('\n');
tclock3 = clock;
time2 = etime(tclock3, tclock2);
fprintf ('\tElapsed time: %.2f + %.2f = %.2f [sec]\n', time1, time2, time1+time2);

