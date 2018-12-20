DataFolder;
[fp, error] = ITCInitializeAnalysis(50000, '012609Ec1');

[NumberEpochs, error] = ITCGetNumberEpochs(fp);

[StartTime, error] = ITCGetEpochTime(0, fp);

for epoch = 0:NumberEpochs-1
    [Comment, error] = ITCGetEpochComment(epoch, fp);
    [Segments, error] = ITCGetEpochSegments(epoch, fp);
    [Time, error] = ITCGetEpochTime(epoch, fp);
    fprintf(1, '*Epoch %d  Time %d sec\n', epoch, Time-StartTime);
    fprintf(1, '     Comment: %s\n', Comment);
    for segment = 1:Segments
        [prepts, error] = ITCGetStmPrePts(epoch, segment-1, 0, fp);
        [stmpts, error] = ITCGetStmPts(epoch, segment-1, 0,  fp);
        [tailpts, error] = ITCGetStmTailPts(epoch, segment-1, 0, fp);
        [amp, error] = ITCGetStmAmp(epoch, segment-1, 0, fp);
        [meanamp, error] = ITCGetStmMean(epoch, segment-1, 0, fp);
        [outputchan, error] = ITCGetOutputChan(epoch, segment-1, fp);
        fprintf(1, '\tSegment %d: channel = %d pre points = %d stm points = %d tail points = %d amp = %d mean = %d\n', segment, outputchan, prepts, stmpts, tailpts, amp, meanamp);
    end
    fprintf(1, '\n');
    
end