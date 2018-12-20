function MI = MutualInfoInJointPdf(JointPdf) 


FundMat = ones(size(JointPdf))*.00001 ; % add tiny fraction to make sure no zeros
JointPdf = (JointPdf+FundMat) ;
JointPdf = JointPdf/sum(JointPdf(:)) ;

Py = repmat(sum(JointPdf,1),[size(JointPdf,1),1]) ;
Px = repmat(sum(JointPdf,2),[1,size(JointPdf,2)]) ;

PxyOverPxPy = JointPdf.*log2(JointPdf./(Px.*Py)) ;
MI = sum(PxyOverPxPy(:)) ; % sum over x,y
