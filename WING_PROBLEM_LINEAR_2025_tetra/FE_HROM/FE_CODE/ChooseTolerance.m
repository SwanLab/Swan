function TOL = ChooseTolerance(CN,COOR) ;


TOL = 1e-3 ; 
i1 = 1 ; 
i2 = 2; 
nodes1 = CN(:,i1) ; 
nodes2 = CN(:,i2) ; 
VECT= COOR(nodes1,:)-COOR(nodes2,:) ; 
normV = sqrt(sum(VECT'.*VECT')) ;

minSIZEelem = min(normV); minSIZEelem = minSIZEelem(1) ; 

TOL = (TOL*minSIZEelem) ;