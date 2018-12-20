% qualitative look at effect of background

% 10^-7 at around 200ml/min
Input(1).cellname = '120423_1' ;
Input(1).label = '10^-^7 around 200ml/min' ;
Input(1).control = [4] ;
Input(1).exp = [5,19] ;

Input(2).cellname = '120423_1' ; 
Input(2).label = '10^-^6 around 200ml/min' ;
Input(2).control = [94] ;
Input(2).exp = [95,114] ;

%pure at 12ml/min
Input(3).cellname = '120430_1' ; 
Input(3).label = 'pure at 12ml/min' ;
Input(3).control = [5] ;
Input(3).exp = [6,10] ;

% pure at 2-3 ml/min
Input(4).cellname = '120430_1' ; 
Input(4).label = 'pure at 2-3 ml/min' ;
Input(4).control = [5] ;
Input(4).exp = [11,15] ;

Input(5).cellname = '120516_1' ; 
Input(5).label = 'pure at 2-3 ml/min' ;
Input(5).control = [8] ;
Input(5).exp = [15,20] ;

Input(6).cellname = '120518_2' ; 
Input(6).label = 'pure at 2-3 ml/min' ;
Input(6).control = [6] ;
Input(6).exp = [13,18] ;

% pure at 0ml/min
Input(7).cellname = '120516_1' ; 
Input(7).label = 'pure at 0ml/min' ;
Input(7).control = [8] ;
Input(7).exp = [9,14] ;

Input(8).cellname = '120518_2' ; 
Input(8).label = 'pure at 0ml/min' ;
Input(8).control = [6] ;
Input(8).exp = [7,12] ;

Input(9).cellname = '120604_1' ; 
Input(9).label = 'pure at 0ml/min' ;
Input(9).control = [69] ;
Input(9).exp = [70,75] ;

% pure at 1ml/min
Input(10).cellname = '120604_1' ; 
Input(10).label = 'pure at 1 ml/min' ;
Input(10).control = [69] ;
Input(10).exp = [76,81] ;

% 10^-4 at 5-5.5ml/min
Input(11).cellname = '120716_1' ; 
Input(11).label = '10^-^4 at 5-5.5ml/min' ;
Input(11).control = [28] ;
Input(11).exp = [29,34] ;

Input(12).cellname = '120813_1' ; 
Input(12).label = '10^-^4 at 5-5.5ml/min' ;
Input(12).control = [35] ;
Input(12).exp = [36,41] ;

Input(13).cellname = '120918_1' ; 
Input(13).label = '10^-^4 at 5-5.5ml/min' ;
Input(13).control = [13] ;
Input(13).exp = [14,19] ;

Input(14).cellname = '120924_1' ; 
Input(14).label = '10^-^4 at 5-5.5ml/min' ;
Input(14).control = [23] ;
Input(14).exp = [24,29] ;

Input(15).cellname = '120927_1' ; 
Input(15).label = '10^-^4 at 5-5.5ml/min' ;
Input(15).control = [15] ;
Input(15).exp = [16,26] ;

Input(16).cellname = '121001_1' ; 
Input(16).label = '10^-^4 at 5-5.5ml/min' ;
Input(16).control = [22] ;
Input(16).exp = [23,33] ;

Input(17).cellname = '121005_1' ; 
Input(17).label = '10^-^4 at 5-5.5ml/min' ;
Input(17).control = [13] ;
Input(17).exp = [14,24] ;

for A=[1:17] ;
    sampleRate = 10000 ;
    
    odorRspTrials = [Input(A).control,Input(A).exp] ;
    
    rootdir = ['Z:\Cafaro Data Backup\', Input(A).cellname(1:6),'Data'] ;
    
    for a=1:length(odorRspTrials) ;
        temp = load([rootdir,'\',Input(A).cellname,'\','voltage_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
        vData(a,:) = temp.voltage ;
    end
    
    figure
    subplot(2,1,1)
    plot([1:length(vData)]/sampleRate,vData) ;
    title([Input(A).cellname(1:end-2),'\',Input(A).cellname(end-1:end),' ',Input(A).label])
    xlabel('time')
    ylabel('voltage (mV)')
    legend('control','exp Early', 'exp Late')
    
    subplot(2,1,2)
    plot([1:length(vData)]/sampleRate,vData) ;
    xlim([1.5:3])
    xlabel('time')
    ylabel('voltage (mV)')
    
    cd 
    saveas(gcf,['fig',num2str(A)],'pdf')
    
    clear vData odorRspTrials
end
    
    
    
    
    


