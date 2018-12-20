function [raw_subtraction,DeltaT, tli_spike, tlj_spike] = DeltaT_from_SCR_Matrix(scr, tli, tlj,cost)

%
%
%
%
%

clear minimum* subtraction

[scr_rows,scr_columns]=size(scr);                                                            % how many rows/columns are in the scr matrix

for i=1:scr_rows; 
    for j=1:scr_columns;
        raw_subtraction(i,j)=scr(i,j)-((i-1)+(j-1));                                            % compute the subtraction matrix (how different each value is than it would be if the remove/replace operation had been performed)
    end
end

scale_factor=cost*.1;                                                                               
subtraction=raw_subtraction./scale_factor;                                      % this line, and the line either side of it, reduces the chance of small decimals messing up 'find' below
subtraction=round(subtraction);                                                             
% 
% if scr_rows>=scr_columns                                                                   % .....if this matrix has more rows than columns
%     row_count=2;                                                                                        % start in the first row of the subtraction matrix
%     column_count=2;                                                                                     % start in the first column of the subtraction matrix
%     for a=2:scr_rows                                                                                         % while the row count is less than or equal to the number of rows    
%         while column_count<=scr_columns;                                                     % while the column count is less than or equal to the number of columns
%             row_min=min(subtraction(row_count,:));                                              % what is the mimum value in this row of the subtraction matrix
%             column_min=min(subtraction(:,column_count));                                    %  what is the mimum value in this column of the subtraction matrix
%             if row_min<column_min;                                                                                                                              % if the min of this row is less than the min of this column...
%                 [minimum1(row_count), minimum2(column_count)]=find(subtraction==row_min,1);                      % save the coordinates where the min occurs
%             elseif row_min>column_min;                                                                                                                            % otherwise, if the min of this row is greater than the min of this column...
%                 [minimum1(row_count), minimum2(column_count)]=find(subtraction==column_min,1);                    % save the coordinates where the min occur                
%             elseif row_min==column_min
%                 minimum1(row_count)=row_count;
%                 minimum2(column_count)=column_count;
%             end 
%             row_count=row_count+1;                                                                                               % step down one row
%             column_count=column_count+1;                                                                                     % step over one column
%         end
%     end
% elseif scr_rows<scr_columns                                                              % ....otherwise, if this matrix has more columns than rows
%     row_count=2;                                                                                        % start in the first row of the subtraction matrix
%     column_count=2;                                                                                     % start in the first column of the subtraction matrix
%     for a=2:scr_columns                                                                                   % while the row count is less than or equal to the number of rows    
%         while row_count<=scr_rows;                                                                   % while the column count is less than or equal to the number of columns
%             row_min=min(subtraction(row_count,:));                                              % what is the mimum value in this row of the subtraction matrix
%             column_min=min(subtraction(:,column_count));                                   %  what is the mimum value in this column of the subtraction matrix
%             if row_min<column_min;                                                                                                                              % if the min of this row is less than the min of this column...
%                 [minimum1(row_count), minimum2(column_count)]=find(subtraction==row_min,1);                      % save the coordinates where the min occurs
%             elseif row_min>column_min;                                                                                                                            % otherwise, if the min of this row is greater than the min of this column...
%                 [minimum1(row_count), minimum2(column_count)]=find(subtraction==column_min,1);                    % save the coordinates where the min occurs
%             elseif row_min==column_min
%                 minimum1(row_count)=row_count;
%                 minimum2(column_count)=column_count;
%             end
%             row_count=row_count+1;                                                                                               % step down one row
%             column_count=column_count+1;                                                                                     % step over one column
%         end
%     end
% end
%    
% %
% %
% %
% %
% 


row_step=2;
column_step=2;
row_place=2;
column_place=2;
pair=1;
while row_place<scr_rows-1;
    while column_place<scr_columns-1;
        [r1,c1]=min(subtraction(row_step,:));
        [r2,c2]=min(subtraction(:,column_step));
        if r1<r2;
            row_place=row_step
            column_place=c1
            [nextrow_min, nextrow_min_place]=min(subtraction(row_place+1,:));
            [nextcolumn_min, nextcolumn_min_place]=min(subtraction(:,column_place+1));    
        elseif r1>r2;
            row_place=c2
            column_place=column_step
            [nextrow_min, nextrow_min_place]=min(subtraction(row_place+1,:));
            [nextcolumn_min, nextcolumn_min_place]=min(subtraction(:,column_place+1));
            if row_place+1==nextcolumn_min_place | column_place+1==nextrow_min_place;
                DeltaT(pair)=abs(tli(row_place-1)-tlj(column_place-1));
                tli_spike(pair)=row_place-1;
                tlj_spike(pair)=column_place-1;
                pair=pair+1;
                row_step=row_step+1;
                column_step=column_step+1;
                clear r1 r2 c1 c2 next*
            else
                row_step=row_step+1;
                column_step=column_step+1;
                clear r1 r2 c1 c2 next*
            end
        end
    end
end

         

        
            
            
        

        
       
        
        























% 
% 
% 
% 
% row_step=2;
% column_step=2;
% row_place=2;
% column_place=2;
% pair=1;
% 
% while row_place<scr_rows;
%     while column_place<scr_columns;
%         [r1,c1]=min(subtraction(row_step,:));
%         [r2,c2]=min(subtraction(:,column_step));
%         if r1<r2;
%             row_place=row_step
%             column_place=c1
%         elseif r1>r2;
%             row_place=c2
%             column_place=column_step
%             if subtraction(row_place-1,column_place-1)==min(subtraction(row_place-1,:)) | subtraction(row_place-1,column_place-1)==min(subtraction(:,column_place-1));
%                 if min(subtraction(row_place,column_place+1))~=min(subtraction(:,column_place+1));
%                     DeltaT(pair)=abs(tli(row_place-1)-tlj(column_place-1));
%                     tli_spike(pair)=row_place-1;
%                     tlj_spike(pair)=column_place-1;
%                    pair=pair+1
%                 end
%             end
%         end
%         row_step=row_step+1;
%         column_step=column_step+1;
%         clear r1 r2 c1 c2
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% row_step=2;
% column_step=2;
% pair=1;
% 
% while row_step<15%scr_rows;
%     while column_step<15%scr_columns;
%         [r1,c1]=min(subtraction(row_step,:))
%         [r2,c2]=min(subtraction(:,column_step))
%         if r1<r2;
%             row_place=row_step
%             column_place=c1
%         elseif r1>r2;
%             row_place=c2
%             column_place=column_step
%             if subtraction(row_place-1,column_place-1)==min(subtraction(row_place-1,:)) | subtraction(row_place-1,column_place-1)==min(subtraction(:,column_place-1))
%                 DeltaT(pair)=abs(tli(row_place-1)-tlj(column_place-1));
%                 tli_spike(pair)=row_place-1
%                 tlj_spike(pair)=column_place-1
%                 pair=pair+1
%             end
%         end
%         row_step=row_step+1;
%         column_step=column_step+1;
%         clear r1 r2 c1 c2
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %             elseif subtraction(row_place-1,column_place-1)==subtraction(row_place,column_place)-2;
% %                 row_step=row_step+1;
% %                 column_step=column_step+1;
% %             elseif subtraction(row_step-1,column_step-1)~=min(subtraction(row_step-1,:));
% %                 row_step=row_step+1;
% %             elseif subtraction(row_step-1,column_step-1)~=min(subtraction(:,column_step-1));
% %                 column_step=column_step+1;
% 
% 
% % after looking at row 2, row 2 (and getting a pair), needs to now look at row 3, column 3 (not +1 from where the pair occurred) 
% 
% 
% 
% 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % tli_spikes=minimum1(2:end)-1;
% % tlj_spikes=minimum2(2:end)-1;
% %       
% % [x1,y1,z1]=unique(tli_spikes)
% % [x2,y2,z2]=unique(tlj_spikes)
% % 
% % a1=length(y1)
% % a2=length(y2)
% % if a1>=a2
% %     a=a2
% % elseif a1<a2
% %     a=a1
% % end
% % 
% % counter=1
% % pair=1
% % 
% % while counter<a
% %         if tli_spikes(counter+1)~=tli_spikes(counter) & tlj_spikes(counter+1)~=tlj_spikes(counter);
% %             DeltaT(pair)=abs(tli(tli_spikes(y1(counter)))-tlj(tlj_spikes(y2(counter))));
% %             counter=counter+1;
% %             pair=pair+1
% %         elseif tli_spikes(counter+1)==tli_spikes(counter) | tlj_spikes(counter+1)==tlj_spikes(counter);
% %             
% %             DeltaT(pair)=abs(tli(tli_spikes(y1(counter)))-tlj(tlj_spikes(y2(counter))));
% %             counter=counter+2;                                                                                                                      % this won't work forever.... need to find out how many in a row are like this number
% %         end
% % end
% % 
% % 
% %         
% %                 
% %             