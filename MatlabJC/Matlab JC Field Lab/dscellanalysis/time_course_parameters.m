%Calculation of Parameters from Temporal Receptive Field
function[TCParams] = time_course_parameters(timecourse, correction)

% Input: Time course matrix (rows - Tc values, columns - each cell) 
%        If correction is needed because of noisy time courses (1 or 0)
% Returns structure with parameters calculated from each cell's time course in this order:

% 1. Standard Deviation
% 2. Variance
% 3. Mean
% 4. Minimum Value
% 5. Minimum Value Time
% 6. Peak Value
% 7. Peak Value Time
% 8. Norm
% 9. Absolute Area Under Curve
%9a. Total Area Under Curve
%9b. Absolute Area under curve from min to max
%9c. Total Area under curve from min to max
%9d. Time to first peak
%9e. Time to second peak
%9f. Absolute Area after 1st peak
%9g. Absolute Area before 2nd peak
%9h. Gradient after 1st peak
%9i. Gradient before 2nd peak
%10. Minimum Value Time*Peak Value Time
%11. Minimum Value Time+Peak Value Time
%12. Peak Value Time - Minimum Value Time
%13. Peak Value Time / Minimum Value Time
%14. Minimum Value*Peak Value 
%15. Minimum Value+Peak Value 
%16. Peak Value - Minimum Value
%17. Peak Value / Minimum Value
%18. Abs(Minimum Value+Peak Value)
%19. Time off dip ends (not always accurate)
%20. Time on dip starts (not always accurate)
%21. Time on dip ends (not always accurate)
%22. Time off dip starts (not always accurate)
%22a.Time of zero crossing
%23. Length of on component (not always accurate)
%24. Length of off component (not always accurate)
%25. Gradient between max and min
%26. Minimum Gradient
%27. Minimum Gradient Time
%28. Peak Gradient
%29. Peak Gradient Time
%30. Minimum Gradient Time*Peak Gradient Time
%31. Minimum Gradient Time+Peak Gradient Time
%32. Peak Gradient Time - Minimum Gradient Time
%33. Peak Gradient Time / Minimum Gradient Time
%34. Minimum Gradient*Peak Gradient 
%35. Minimum Gradient+Peak Gradient 
%36. Peak Gradient - Minimum Gradient
%37. Peak Gradient / Minimum Gradient
%38. Abs(Minimum Gradient+Peak Gradient)
%39. Total Area Under Curve / Absolute Area Under Curve (Degree of Transience)




TCParams = struct;

TCParams.stddev = std(timecourse);
TCParams.var = var(timecourse);
TCParams.ave = mean(timecourse);
[TCParams.minval TCParams.mintim] =  min(timecourse);
[TCParams.maxval TCParams.maxtim] =  max(timecourse);


for i = 1:size(timecourse, 2) %or nonds
  TCParams.norm(1, i) = norm( timecourse(:,i)); %Calculate norm (magnitude) for all time courses
  TCParams.areaabs(1,i) = trapz(abs(1:(size(timecourse,1))), abs(timecourse(:,i))); %Calculate Area Under Curve for all time courses
  TCParams.area(1,i) = trapz(abs(1:(size(timecourse,1))), timecourse(:,i));
  if(TCParams.mintim(1,i) < TCParams.maxtim(1,i))
      TCParams.peakminarea(1,i) = trapz(abs(TCParams.mintim(1,i):1:TCParams.maxtim(1,i)), abs(timecourse(TCParams.mintim(1,i):1:TCParams.maxtim(1,i),i)));
      TCParams.peakmintotarea(1,i) = trapz(abs(TCParams.mintim(1,i):1:TCParams.maxtim(1,i)), timecourse(TCParams.mintim(1,i):1:TCParams.maxtim(1,i),i));
      TCParams.firstpeaktime(1,i) = TCParams.maxtim(1,i);
      TCParams.secondpeaktime(1,i)= TCParams.mintim(1,i);      
  else
      TCParams.peakminarea(1,i) = trapz(abs(TCParams.maxtim(1,i):1:TCParams.mintim(1,i)), abs(timecourse(TCParams.maxtim(1,i):1:TCParams.mintim(1,i),i)));
      TCParams.peakmintotarea(1,i) = trapz(abs(TCParams.maxtim(1,i):1:TCParams.mintim(1,i)), timecourse(TCParams.maxtim(1,i):1:TCParams.mintim(1,i),i));
      TCParams.firstpeaktime(1,i) = TCParams.mintim(1,i);
      TCParams.secondpeaktime(1,i)= TCParams.maxtim(1,i);        
  end
  TCParams.areaaft1stpeak(1,i) = trapz(abs(TCParams.firstpeaktime(1,i):1:size(timecourse,1)), abs(timecourse(TCParams.firstpeaktime(1,i):size(timecourse,1),i)));
  if(TCParams.secondpeaktime(1,i) == 1) 
      TCParams.areabef2ndpeak(1,i) = 0;
      TCParams.gradbef2ndpeak(1,i) = 0;
  else
   TCParams.areabef2ndpeak(1,i) = trapz(abs(1:1:TCParams.secondpeaktime(1,i)), abs(timecourse(1:TCParams.secondpeaktime(1,i),i)));
   TCParams.gradbef2ndpeak(1,i) = timecourse(TCParams.secondpeaktime(1,i),i) - timecourse(TCParams.secondpeaktime(1,i)-1,i);
  end
  TCParams.gradaft1stpeak(1,i) = timecourse(TCParams.firstpeaktime(1,i)+1,i) - timecourse(TCParams.firstpeaktime(1,i),i);
end 

if (correction)  
 %TCParams.zerocrossing(1,:) = TCParams.mintstart(1,:);
TCParams.firstpeaktime(1,:) = TCParams.mintim(1,:);

for i = 1:length(TCParams.mintim)    
    for j = TCParams.mintim(1,i):-1:2
        if(timecourse(j,i) < 0 &&  timecourse(j-1,i) >= 0)
            TCParams.mintstart(1,i) = j;
            break;
        end
    end  
%    [spvalcorrection sptimcorrection] = max(timecourse(1:TCParams.mintstart(1,i),i));
    %TCParams.secondpeaktime(1,i) = sptimcorrection;
    %TCParams.maxtim(1,i) = sptimcorrection;
   % TCParams.maxval(1,i) = spvalcorrection;
end

for i = 1:length(TCParams.mintim)
     TCParams.areaaft1stpeak(1,i) = trapz(abs(TCParams.firstpeaktime(1,i):1:size(timecourse,1)), abs(timecourse(TCParams.firstpeaktime(1,i):size(timecourse,1),i)));
   %  TCParams.areabef2ndpeak(1,i) = trapz(abs(1:1:TCParams.secondpeaktime(1,i)), abs(timecourse(1:TCParams.secondpeaktime(1,i),i)));     
   %   TCParams.peakminarea(1,i) = trapz(abs(TCParams.secondpeaktime(1,i):1:TCParams.firstpeaktime(1,i)), abs(timecourse(TCParams.secondpeaktime(1,i):1:TCParams.firstpeaktime(1,i),i)));
    %  TCParams.peakmintotarea(1,i) = trapz(abs(TCParams.secondpeaktime(1,i):1:TCParams.firstpeaktime(1,i)), timecourse(TCParams.secondpeaktime(1,i):1:TCParams.firstpeaktime(1,i),i));
    %  TCParams.gradaft1stpeak(1,i) = timecourse(TCParams.firstpeaktime(1,i)+3,i) - timecourse(TCParams.firstpeaktime(1,i),i);
  %    TCParams.gradbef2ndpeak(1,i) = timecourse(TCParams.secondpeaktime(1,i),i) - timecourse(TCParams.secondpeaktime(1,i)-3,i);   

end

TCParams = rmfield(TCParams,'mintstart');
end  


TCParams.mintmaxtmult = TCParams.mintim.*TCParams.maxtim;
TCParams.mintmaxtadd = TCParams.mintim+TCParams.maxtim;
TCParams.mintmaxtsub = TCParams.maxtim - TCParams.mintim;
TCParams.mintmaxtdiv = TCParams.maxtim./TCParams.mintim;

TCParams.minmaxmult = TCParams.minval.*TCParams.maxval;
TCParams.minmaxadd = TCParams.minval+TCParams.maxval;
TCParams.minmaxsub = TCParams.maxval - TCParams.minval;
TCParams.minmaxdiv = TCParams.maxval./TCParams.minval;
TCParams.minmaxabsadd = abs(TCParams.minval) + abs(TCParams.maxval);

TCParams.mintstart = zeros(1,length(TCParams.mintim));
TCParams.mintend = 30*ones(1,length(TCParams.mintim));
TCParams.maxtstart = zeros(1,length(TCParams.mintim));
TCParams.maxtend = 30*ones(1,length(TCParams.mintim));

for i = 1:length(TCParams.mintim)
    for j = TCParams.mintim(1,i):-1:2
        if(timecourse(j,i) < 0 &&  timecourse(j-1,i) >= 0)
            TCParams.mintstart(1,i) = j;
            break;
        end
    end
    for k = TCParams.mintim(1,i):1:29
        if(timecourse(k,i) < 0 &&  timecourse(k+1,i) >= 0)
            TCParams.mintend(1,i) = k;
            break;
        end
    end   
    for m = TCParams.maxtim(1,i):-1:2
        if(timecourse(m,i) >= 0 &&  timecourse(m-1,i) < 0)
            TCParams.maxtstart(1,i) = m;
            break;
        end
    end
    for n = TCParams.maxtim(1,i):1:29
        if(timecourse(n,i) >= 0 &&  timecourse(n+1,i) < 0)
            TCParams.maxtend(1,i) = n;
            break;
        end
    end 
     for p = TCParams.firstpeaktime(1,i):-1:(TCParams.secondpeaktime(1,i)-1)
         if((timecourse(p,i) >= 0 &&  timecourse(p-1,i) < 0) || (timecourse(p,i) < 0 &&  timecourse(p-1,i) >= 0))
             TCParams.zerocrossing(1,i) = p;
             break;
         end
     end
    
end
    

TCParams.onlength = TCParams.maxtend - TCParams.maxtstart;
TCParams.offlength = TCParams.mintend - TCParams.mintstart; 

TCParams.maxmingrad = (TCParams.maxval - TCParams.minval)./(TCParams.maxtim - TCParams.mintim);


[TCParams.mingradval TCParams.mingradtim] =  min(diff(timecourse)); %times seem very consistent across a type
[TCParams.maxgradval TCParams.maxgradtim] =  max(diff(timecourse)); %times seem very consistent across a type


TCParams.mingradtmaxgradtmult = TCParams.mingradtim.*TCParams.maxgradtim;
TCParams.mingradtmaxgradtadd = TCParams.mingradtim+TCParams.maxgradtim;
TCParams.mingradtmaxgradtsub = TCParams.maxgradtim - TCParams.mingradtim;
TCParams.mingradtmaxgradtdiv = TCParams.maxgradtim./TCParams.mingradtim;

TCParams.mingradmaxgradmult = TCParams.mingradval.*TCParams.maxgradval;
TCParams.mingradmaxgradadd = TCParams.mingradval+TCParams.maxgradval;
TCParams.mingradmaxgradsub = TCParams.maxgradval - TCParams.mingradval;
TCParams.mingradmaxgraddiv = TCParams.maxgradval./TCParams.mingradval;
TCParams.mingradmaxgradabsadd = abs(TCParams.mingradval) + abs(TCParams.maxgradval);

TCParams.dot = TCParams.area ./ TCParams.areaabs; 
  

end