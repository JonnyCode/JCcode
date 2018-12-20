function Interpolated = interpolateNbyFour(original, StartWaveLen, EndWaveLen, StepWaveLen)

numpnts = size(original, 1);

for pnt = 1:EndWaveLen - StartWaveLen + 1
	CurWaveLen = StartWaveLen + (pnt - 1)*StepWaveLen;
	next = 1;
	while ((original(next, 1) < CurWaveLen) & (next < numpnts))
		next = next + 1;
	end
	prev = next - 1;
	weight = (CurWaveLen - original(prev, 1)) / (original(next, 1) - original(prev, 1));
	Interpolated(pnt, 1) = CurWaveLen;
	Interpolated(pnt, 2) = original(prev, 2) * (1 - weight) + original(next, 2) * weight;
	Interpolated(pnt, 3) = original(prev, 3) * (1 - weight) + original(next, 3) * weight;
	Interpolated(pnt, 4) = original(prev, 4) * (1 - weight) + original(next, 4) * weight;
end
