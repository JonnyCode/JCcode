function h = plot1ras(st1,trialnum,ticksize, tickcolor,tickwidth)
% h = plot1ras(st1,trialnum,ticksize,tickcolor,tickwidth)
% plots a set of ticks (raster) at the times in st1.  The height (ordinate)
% value is trialnum.  The tick size is ticksize (duh!).  Tickcolor is a
% standard color (e.g., 'r') or an rgb 3-vector.  Tickwidth is the width.
% The function returns a handle to the raster, so these parameters (and
% others) can be set using set(h,'Property',value)

% 4/26/98 mns & mem wrote it

if nargin < 3
  ticksize = .5;
end

if nargin < 2
	trialnum = 1;
end		

x = st1(isfinite(st1));
x = x(:);
a = [x x nan*ones(size(x))]';
b = repmat([trialnum-ticksize/2 trialnum+ticksize/2 nan],size(a,2),1)';
h = plot(a(:),b(:));

  
if nargin>3
  if isfinite(tickcolor)
    set(h,'Color',tickcolor)
  end
end
if nargin>4
  if isfinite(tickwidth)
    set(h,'Linewidth',tickwidth)
  end
end

if nargin == 1
	% this is a fast call, so let's help out
	set(gca,'YLim',[.5 1.5])
end
