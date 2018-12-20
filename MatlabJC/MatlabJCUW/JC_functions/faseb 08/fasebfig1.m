% figure 1 for FASEB poster '08

% reproduced quick plotter
CellName_str = '052908Bc1' ;
epochs_str = {'[333:8:405]','[583:2:601]','[903:2:921]','[433:8:505]','[603:2:621]','[843:2:861]','[332:8:405]','[582:2:601]','[902:2:921]','[432:8:505]','[602:2:621]','[842:2:861]'} ;
Color = {'b','r','g','k','y','c','m','b-.','r-.','g-.','k-.','y-.','c-.','m-.','b:','r:','g:','k:'} ;
fig_num = 1;

[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',CellName_str]);


for set = 1:length(epochs_str) ;
    epochNums{set} = str2num(epochs_str{set}) ;
    for trace = 1:length(epochNums{set}) ;
        [Data{set}(trace,:), error] = ITCReadEpoch(epochNums{set}(trace), 0, fp);
    end
    Trace_mean{set} = mean(Data{set},1) ;
    Trace_mean{set} = Trace_mean{set} - mean(Trace_mean{set}(1)) ;
    
    figure(fig_num)
    plot(Trace_mean{set},Color{set})
    hold on
end

legend(epochs_str)
title(CellName_str)

% change into conductances (assume Ecl = -61.6 and Eampa = 0)
% name each vector in a structure for igor

ForIgor.Epochs = epochs_str ;
ForIgor.ExcInc = Trace_mean{1}/-61.6 ;
ForIgor.ExcIncAPB = Trace_mean{2}/-61.6 ;
ForIgor.ExcIncWash = Trace_mean{3}/-61.6 ;
ForIgor.InhInc = Trace_mean{4}/61.6 ;
ForIgor.InhIncAPB = Trace_mean{5}/61.6 ;
ForIgor.InhIncWash = Trace_mean{6}/61.6 ;
ForIgor.ExcDec = Trace_mean{7}/-61.6 ;
ForIgor.ExcDecAPB = Trace_mean{8}/-61.6 ;
ForIgor.ExcDecWash = Trace_mean{9}/-61.6 ;
ForIgor.InhDec = Trace_mean{10}/61.6 ;
ForIgor.InhDecAPB = Trace_mean{11}/61.6 ;
ForIgor.InhDecWash = Trace_mean{12}/61.6 ;

% export this structure as a HDF5 which will be read into igor (will write
% in current directory)
exportStructToHDF5(ForIgor, 'StepResonse', '~/data_analysis/For Igor') ;
















