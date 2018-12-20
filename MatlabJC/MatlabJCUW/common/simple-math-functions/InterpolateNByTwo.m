function Interpolated = InterpolateNByTwo(original, StartWaveLen, EndWaveLen, StepWaveLen)
%
%	function Interpolated = InterpolateNByTwo(original, StartWaveLen, EndWaveLen, StepWaveLen)
%
% This interpolate function works on a nx2 matrix.  The first column is a counter, keeping track
% of time or some other index.  The second column is the data.  This will interpolate between the 
% points in the data and set the appropriate indexes.
% 
%   Revised 5/05 FMR: reordered search to speed up

numpnts = size(original, 1);

OutputPts = (EndWaveLen - StartWaveLen) / StepWaveLen + 1;
for pnt = 1:OutputPts
	CurWaveLen = StartWaveLen + (pnt - 1)*StepWaveLen;
    next = 2;
	while ((original(next, 1) < CurWaveLen) & (next < numpnts))
		next = next + 1;
	end
	prev = next - 1;
	weight = (CurWaveLen - original(prev, 1)) / (original(next, 1) - original(prev, 1));
	Interpolated(pnt, 1) = CurWaveLen;
	Interpolated(pnt, 2) = original(prev, 2) * (1 - weight) + original(next, 2) * weight;
end
