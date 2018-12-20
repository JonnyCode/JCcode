function [movScrambled,permStruct] = movScrambler(mov,SpatialScrambleFlag, TemporalScrambleFlag, DoubleScrambleFlag)

% JC 10/7/2016

% This function will scramble the stixel positions and/or frame
% positions,of mov(x,y,t).  Ostensibly, this will spatially and/or temporally whiten
% the movie without changing intensity distribution or the spatial or
% temporal power spectrum.

movScrambled = nan(size(mov)) ; 

xl = size(mov,1) ;
yl = size(mov,2) ;
fl = size(mov,3) ;

if SpatialScrambleFlag ;
    permStruct.xyperm = randperm(xl*yl) ; 
    [meshX,meshY]=meshgrid([1:xl],[1:yl]) ;
    s=1 ; % stixel count
    for x=1:xl ;
        for y=1:yl ;
            movScrambled(x,y,:)= mov(meshX(permStruct.xyperm(s)), meshY(permStruct.xyperm(s)),:) ;
            s=s+1 ;
        end
    end
end
        
if TemporalScrambleFlag ;
    permStruct.fperm = randperm(fl) ;
    if SpatialScrambleFlag ;
        mov=movScrambled ;
    end
    
    for f=1:fl ;
        movScrambled(:,:,f)= mov(:,:,permStruct.fperm(f)) ;
    end
end        
        
if DoubleScrambleFlag ; % if you want to scramble all stixels in space and time (different from running both above flags)
    permStruct.xyfperm = randperm(xl*yl*fl) ;  
    [meshX,meshY,meshF]=meshgrid([1:xl],[1:yl],[1:fl]) ;
    s=1 ; % stixel count
    for x=1:xl ;
        for y=1:yl ;
            for f=1:fl ;
                movScrambled(x,y,f)= mov(meshX(permStruct.xyfperm(s)), meshY(permStruct.xyfperm(s)),meshF(permStruct.xyfperm(s))) ;
                s=s+1 ;
            end
        end
    end
end
