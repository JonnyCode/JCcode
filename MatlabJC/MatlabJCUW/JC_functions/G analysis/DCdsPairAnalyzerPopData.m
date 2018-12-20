% analysis of population statistics for dynamic clamp paired data
% JC 4/27/11

% paired t-test of +c+p vs. -c+p conditions
[h,ForIgor.psn141Bar1All] = ttest(ForIgor.snCorrsdc141Bar1All(:,1),ForIgor.snCorrsdc141Bar1All(:,2)) ;
[h,ForIgor.psn141Bar2All] = ttest(ForIgor.snCorrsdc141Bar2All(:,1),ForIgor.snCorrsdc141Bar2All(:,2)) ;
[h,ForIgor.psn141Bar3All] = ttest(ForIgor.snCorrsdc141Bar3All(:,1),ForIgor.snCorrsdc141Bar3All(:,2)) ;

[h,ForIgor.pst141Bar1All] = ttest(ForIgor.stCorrsdc141Bar1All(:,1),ForIgor.stCorrsdc141Bar1All(:,2)) ;
[h,ForIgor.pst141Bar2All] = ttest(ForIgor.stCorrsdc141Bar2All(:,1),ForIgor.stCorrsdc141Bar2All(:,2)) ;
[h,ForIgor.pst141Bar3All] = ttest(ForIgor.stCorrsdc141Bar3All(:,1),ForIgor.stCorrsdc141Bar3All(:,2)) ;

% mean
ForIgor.snCorrsdc141Bar1Mean = mean(ForIgor.snCorrsdc141Bar1All) ;
ForIgor.snCorrsdc141Bar2Mean = mean(ForIgor.snCorrsdc141Bar2All) ;
ForIgor.snCorrsdc141Bar3Mean = mean(ForIgor.snCorrsdc141Bar3All) ;

ForIgor.stCorrsdc141Bar1Mean = mean(ForIgor.stCorrsdc141Bar1All) ;
ForIgor.stCorrsdc141Bar2Mean = mean(ForIgor.stCorrsdc141Bar2All) ;
ForIgor.stCorrsdc141Bar3Mean = mean(ForIgor.stCorrsdc141Bar3All) ;

% sem
ForIgor.snCorrsdc141Bar1SEM = std(ForIgor.snCorrsdc141Bar1All)/sqrt(size(ForIgor.snCorrsdc141Bar1All,1)) ;
ForIgor.snCorrsdc141Bar2SEM = std(ForIgor.snCorrsdc141Bar2All)/sqrt(size(ForIgor.snCorrsdc141Bar2All,1)) ;
ForIgor.snCorrsdc141Bar3SEM = std(ForIgor.snCorrsdc141Bar3All)/sqrt(size(ForIgor.snCorrsdc141Bar3All,1)) ;

ForIgor.stCorrsdc141Bar1SEM = std(ForIgor.stCorrsdc141Bar1All)/sqrt(size(ForIgor.stCorrsdc141Bar1All,1)) ;
ForIgor.stCorrsdc141Bar2SEM = std(ForIgor.stCorrsdc141Bar2All)/sqrt(size(ForIgor.stCorrsdc141Bar2All,1)) ;
ForIgor.stCorrsdc141Bar3SEM = std(ForIgor.stCorrsdc141Bar3All)/sqrt(size(ForIgor.stCorrsdc141Bar3All,1)) ;


% paired t-test of +c+p vs. -c+p conditions
[h,ForIgor.psn139Bar1All] = ttest(ForIgor.snCorrsdc139Bar1All(:,1),ForIgor.snCorrsdc139Bar1All(:,2)) ;
[h,ForIgor.psn139Bar2All] = ttest(ForIgor.snCorrsdc139Bar2All(:,1),ForIgor.snCorrsdc139Bar2All(:,2)) ;
[h,ForIgor.psn139Bar3All] = ttest(ForIgor.snCorrsdc139Bar3All(:,1),ForIgor.snCorrsdc139Bar3All(:,2)) ;

[h,ForIgor.pst139Bar1All] = ttest(ForIgor.stCorrsdc139Bar1All(:,1),ForIgor.stCorrsdc139Bar1All(:,2)) ;
[h,ForIgor.pst139Bar2All] = ttest(ForIgor.stCorrsdc139Bar2All(:,1),ForIgor.stCorrsdc139Bar2All(:,2)) ;
[h,ForIgor.pst139Bar3All] = ttest(ForIgor.stCorrsdc139Bar3All(:,1),ForIgor.stCorrsdc139Bar3All(:,2)) ;

% mean
ForIgor.snCorrsdc139Bar1Mean = mean(ForIgor.snCorrsdc139Bar1All) ;
ForIgor.snCorrsdc139Bar2Mean = mean(ForIgor.snCorrsdc139Bar2All) ;
ForIgor.snCorrsdc139Bar3Mean = mean(ForIgor.snCorrsdc139Bar3All) ;

ForIgor.stCorrsdc139Bar1Mean = mean(ForIgor.stCorrsdc139Bar1All) ;
ForIgor.stCorrsdc139Bar2Mean = mean(ForIgor.stCorrsdc139Bar2All) ;
ForIgor.stCorrsdc139Bar3Mean = mean(ForIgor.stCorrsdc139Bar3All) ;

% sem
ForIgor.snCorrsdc139Bar1SEM = std(ForIgor.snCorrsdc139Bar1All)/sqrt(size(ForIgor.snCorrsdc139Bar1All,1)) ;
ForIgor.snCorrsdc139Bar2SEM = std(ForIgor.snCorrsdc139Bar2All)/sqrt(size(ForIgor.snCorrsdc139Bar2All,1)) ;
ForIgor.snCorrsdc139Bar3SEM = std(ForIgor.snCorrsdc139Bar3All)/sqrt(size(ForIgor.snCorrsdc139Bar3All,1)) ;

ForIgor.stCorrsdc139Bar1SEM = std(ForIgor.stCorrsdc139Bar1All)/sqrt(size(ForIgor.stCorrsdc139Bar1All,1)) ;
ForIgor.stCorrsdc139Bar2SEM = std(ForIgor.stCorrsdc139Bar2All)/sqrt(size(ForIgor.stCorrsdc139Bar2All,1)) ;
ForIgor.stCorrsdc139Bar3SEM = std(ForIgor.stCorrsdc139Bar3All)/sqrt(size(ForIgor.stCorrsdc139Bar3All,1)) ;


% paired t-test of +c+p vs. -c+p conditions
[h,ForIgor.psn157Bar1All] = ttest(ForIgor.snCorrsdc157Bar1All(:,1),ForIgor.snCorrsdc157Bar1All(:,2)) ;
[h,ForIgor.psn157Bar2All] = ttest(ForIgor.snCorrsdc157Bar2All(:,1),ForIgor.snCorrsdc157Bar2All(:,2)) ;
[h,ForIgor.psn157Bar3All] = ttest(ForIgor.snCorrsdc157Bar3All(:,1),ForIgor.snCorrsdc157Bar3All(:,2)) ;

[h,ForIgor.pst157Bar1All] = ttest(ForIgor.stCorrsdc157Bar1All(:,1),ForIgor.stCorrsdc157Bar1All(:,2)) ;
[h,ForIgor.pst157Bar2All] = ttest(ForIgor.stCorrsdc157Bar2All(:,1),ForIgor.stCorrsdc157Bar2All(:,2)) ;
[h,ForIgor.pst157Bar3All] = ttest(ForIgor.stCorrsdc157Bar3All(:,1),ForIgor.stCorrsdc157Bar3All(:,2)) ;

% mean
ForIgor.snCorrsdc157Bar1Mean = mean(ForIgor.snCorrsdc157Bar1All) ;
ForIgor.snCorrsdc157Bar2Mean = mean(ForIgor.snCorrsdc157Bar2All) ;
ForIgor.snCorrsdc157Bar3Mean = mean(ForIgor.snCorrsdc157Bar3All) ;

ForIgor.stCorrsdc157Bar1Mean = mean(ForIgor.stCorrsdc157Bar1All) ;
ForIgor.stCorrsdc157Bar2Mean = mean(ForIgor.stCorrsdc157Bar2All) ;
ForIgor.stCorrsdc157Bar3Mean = mean(ForIgor.stCorrsdc157Bar3All) ;

% sem
ForIgor.snCorrsdc157Bar1SEM = std(ForIgor.snCorrsdc157Bar1All)/sqrt(size(ForIgor.snCorrsdc157Bar1All,1)) ;
ForIgor.snCorrsdc157Bar2SEM = std(ForIgor.snCorrsdc157Bar2All)/sqrt(size(ForIgor.snCorrsdc157Bar2All,1)) ;
ForIgor.snCorrsdc157Bar3SEM = std(ForIgor.snCorrsdc157Bar3All)/sqrt(size(ForIgor.snCorrsdc157Bar3All,1)) ;

ForIgor.stCorrsdc157Bar1SEM = std(ForIgor.stCorrsdc157Bar1All)/sqrt(size(ForIgor.stCorrsdc157Bar1All,1)) ;
ForIgor.stCorrsdc157Bar2SEM = std(ForIgor.stCorrsdc157Bar2All)/sqrt(size(ForIgor.stCorrsdc157Bar2All,1)) ;
ForIgor.stCorrsdc157Bar3SEM = std(ForIgor.stCorrsdc157Bar3All)/sqrt(size(ForIgor.stCorrsdc157Bar3All,1)) ;

