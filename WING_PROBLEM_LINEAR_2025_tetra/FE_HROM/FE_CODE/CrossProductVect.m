function t1t2 = CrossProductVect(t1,t2,ndim)

t1 =reshape(t1,ndim,[]) ;
t2 =reshape(t2,ndim,[]) ;
t1t2 = (cross(t1',t2',2))';
t1t2 = t1t2(:) ;