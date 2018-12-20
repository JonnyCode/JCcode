function [OutputDist,OutputDistX] = pdfTransform(InputDist,InputDistX, Transform)

% this function will transform a set of input probability density functions
% via a monotonic function into a corresponding set of output pdfs.
% JC 10/4/12
% InputDist should be a matrix with pdfs in rows
% InputDistX should be a row vector (x values of input pdfs)
% Transform should be a row vector the same length as InputDistX

%output distributions 
OutputDist = nan(2,length(InputDistX)) ;
OutputDistX = Transform ;
a=1 ;
b=1 ;
while a<=length(InputDistX) ;
    if a==1 
        i = find(Transform(a+1:end)~=Transform(a),1,'first') ;
        if isempty(i) ;
            i = length(Transform)-a+1 ;
        end
        if i>1 ;
            weightingFnct = [1,2*ones(1,i-1),1] ; 
        else
            weightingFnct = 1 ;
        end
        OutputDist(:,b) = (InputDist(:,a:a-1+i)*weightingFnct'*(InputDistX(2)-InputDistX(1))/Transform(a)) ; 
        a = a + i  ;
        
    elseif a~=length(InputDistX) ; 
        i = find(Transform(a+1:end)~=Transform(a),1,'first') ;
        if isempty(i) ;
            i = length(Transform)-a+1 ;
        end
        weightingFnct = [1,2*ones(1,i-1),1] ; 
        OutputDist(:,b) = (InputDist(:,a-1:a-1+i)*weightingFnct'*(InputDistX(2)-InputDistX(1))/(Transform(a)-Transform(a-1)))-OutputDist(:,b-1) ;
        a = a + i  ;
    end
    
    b=b+1 ;
end