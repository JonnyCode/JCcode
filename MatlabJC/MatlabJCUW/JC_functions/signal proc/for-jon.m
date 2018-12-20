TargetCovar = [5.8920,4.2989; 4.2989, 7.8981];
TargetMeans = [1.7, 3.1];

minerr = 1e6;
for var1 = 1:1:10 ;
    for var2 = 1:1:10 ;
        for cc = 0:.1:.9 ;
            for mean1 = -10:1:4 ;
                for mean2 = -10:1:4 ;
                    covar = cc*sqrt(var1*var2) ;
                    TestCov = [var1 covar; covar var2];
                    [NewCovar, NewMeans] = CorrectCovarForRect(TestCov, [mean1 mean2], [0 0]);
                    err = sum((NewCovar([1 3 4]) - TargetCovar([1 3 4])).^2) + sum((NewMeans - TargetMeans).^2);
                    if (err < minerr)
                        NewCovar ;
                        NewMeans ;
                        %fprintf(1, 'err = %d (%d %d %d %d %d)\n', err, var1, var2, covar, mean1, mean2);
                        minerr = err;
                        finalVar1 = var1 ;
                        finalVar2 = var2 ;
                        finalCov = covar ;
                        finalMean1 = mean1 ;
                        finalMean2 = mean2 ;
                    end
                end
            end
        end
    end
end
