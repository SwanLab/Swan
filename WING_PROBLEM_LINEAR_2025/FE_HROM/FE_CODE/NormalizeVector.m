function t1 = NormalizeVector(t1,ndim)

normt1 = reshape(t1,ndim,[]) ;
normt1 = sqrt(sum(normt1.*normt1))  ; 
normt1 = repmat(normt1,ndim,1);
normt1 = normt1(:);
t1 = t1./normt1 ;