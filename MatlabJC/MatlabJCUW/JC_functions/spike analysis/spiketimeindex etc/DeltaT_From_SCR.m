function [DeltaT, tli_spike, tlj_spike] = DeltaT_From_SCR(scr,cost,tli,tlj)

%
%  determines which spikes in two spike trains (tli and tlj) were paired by the spk_d
%  algorithim (for a given 'cost' of shifting spikes in time). 
%
%
% flawed approach called Simple_DeltaT_from_SCR. don't use!
% early version of this successful approach called 'New_Approach'
%
% GJM  3/06



%
% compute the subtraction matrix... the distance that you have saved by
% making a given pairing (relative to removing/replacing a spike pair and all
% others considered before it).
%

[scr_rows,scr_columns]=size(scr);          

for i=1:scr_rows;                                                                                           % for each row of the scr matrix
    for j=1:scr_columns;                                                                                    % for each column of the scr matrix 
        raw_subtraction(i,j)=scr(i,j)-((i-1)+(j-1));                                                  % compute the subtraction matrix (how different each value is than it would be if the remove/replace operation had been performed)
    end
end

scale_factor=cost*.1;                                                                            %                                           
subtraction=raw_subtraction./scale_factor;                                        % this line, and the line either side of it, reduces the chance of small decimals messing up the min command below
subtraction=round(subtraction);                                                              %

%
% determine the number of rows/columns in the scr (and subtraction) matrix, then
% start at bottom right corner of the subtraction matrix. from each
% position, move up or left (according to whatever movement enables the
% most negative subtraction value). if values up and left in subtraction
% matrix from your current position are equal, move diagonally (but only
% keep DeltaT value if change in value (in subtraction matrix) is <2)
%

row=scr_rows;                                                                                                % start at highest row number
column=scr_columns;                                                                                     % start at highest column number
pair=1;                             

last_pair_row=scr_rows+1;                                                                               % imagine there is another pair so that you can get the last 
last_pair_column=scr_columns+1;                                                                    % imagine there is another pair so that you can get the last 

while (row>1 & column>1)                                                                         
    [left_min]=subtraction(row,column-1);                                                                                                 % what is the value in the current row but over one column 
    [up_min]=subtraction(row-1,column);                                                                                                   % what is the value in the current column but up one row
    if left_min==up_min;                                                                                                                                   % if the values up and left from your current position are equal...     
        if (row<last_pair_row & column<last_pair_column & subtraction(row,column)~=subtraction(row-1,column-1));            % if the value in the diagonal (up/left) is not the same as in your current position.... and you are over at least one row and one column from where the last pair occurred    
            DeltaT(pair)=abs(tli(row-1)-tlj(column-1));                                                                                                                                     % compute the delta T for this pair
            tli_spike(pair)=row-1;                                                                                                                                                                            % this spike from spike train 'tli' forms one half of the pair
            tlj_spike(pair)=column-1;                                                                                                                                                                        % this spike from spike train 'tlj' forms the other half
            last_pair_row=row;                                                                                                                                                                                     % this pair is now the most recent formed.... can't pick another pair that contains this spike from tli
            last_pair_column=column;                                                                                                                                                                           % this pair is now the most recent formed.... can't pick another pair that contains this spike from tlj
            pair=pair+1;                                                                                                                                                                                                        % update index counter
            row=row-1;                                                                                                                                                                                                            % move up one row
            column=column-1;                                                                                                                                                                                                 % move over one column
        elseif (column==last_pair_column | row==last_pair_row | subtraction(row,column)==subtraction(row-1,column-1)) ;     % if one spike from the current spike pairing has already been asigned to another pair, or spike taken out and replaced (thus no difference in subtraction), you can't assign/determine Delta T
            if left_min<up_min                                                                                                                                                                                      % never going to be true...                                                  
                row=row-1;                                                                                                                                                                                                  % good thing, because this is the wrong movement!
            elseif left_min>up_min                                                                                                                                                                                % never going to be true...         
                column=column-1;                                                                                                                                                                                     % good thing, because this is the wrong movement!
            elseif left_min==up_min
                row=row-1;
                column=column-1;
            end               
        end           
    elseif up_min<left_min;
        if (row<last_pair_row & column<last_pair_column & up_min~=subtraction(row,column));
            DeltaT(pair)=abs(tli(row-1)-tlj(column-1));                     
            tli_spike(pair)=row-1;
            tlj_spike(pair)=column-1;
            last_pair_row=row;
            last_pair_column=column;
             pair=pair+1;
            row=row-1;
        elseif (row==last_pair_row | column==last_pair_column | up_min==subtraction(row,column)); 
            row=row-1 ;      % not sure that this is going to work in all circumstances
        end
    elseif left_min<up_min;            
        if (row<last_pair_row & column<last_pair_column & left_min~=subtraction(row,column));
            DeltaT(pair)=abs(tli(row-1)-tlj(column-1));                     
            tli_spike(pair)=row-1;
            tlj_spike(pair)=column-1;
            last_pair_row= row;
            last_pair_column= column;
             pair=pair+1;
            column=column-1;
        elseif (column==last_pair_column | row==last_pair_row | left_min==subtraction(row,column));
            column=column-1;       % not sure that this is going to work in all circumstances
        end
    end
end

              
            
            
            
            
            
            