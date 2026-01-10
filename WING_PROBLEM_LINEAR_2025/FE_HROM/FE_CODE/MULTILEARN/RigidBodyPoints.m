function DOFp = RigidBodyPoints(PointsRef,COOR)
% See ALEX_THESIS_mine.pdf
% Select DOFs to prescribe so that rigid body motions are avoided

ndim = size(COOR,2) ;
DOFp1 = Nod2DOF(PointsRef(1),ndim) ;
if ndim == 3


DOFp23 = Nod2DOF(PointsRef(2:3),ndim) ;
c1 = COOR(PointsRef(2),:)-COOR(PointsRef(1),:) ; 
c2 = COOR(PointsRef(3),:)-COOR(PointsRef(1),:) ; 
Omega1 = [0 -c1(3) c1(3) ; c1(3) 0 -c1(1); -c1(2) c1(1) 0] ; 
Omega2 = [0 -c2(3) c2(3) ; c2(3) 0 -c2(1); -c2(2) c2(1) 0] ; 
Omega = [Omega1 Omega2] ; 
[~,idx]=licols(Omega) ; 
DOFp23 = DOFp23(idx) ; 
DOFp = [DOFp1 ; DOFp23] ;

else
    DOFp1 = Nod2DOF(PointsRef(1),ndim) ;
    DOFp2 = Nod2DOF(PointsRef(2),ndim) ;
    
    DOFp = [DOFp1; DOFp2(1)] ; 

end