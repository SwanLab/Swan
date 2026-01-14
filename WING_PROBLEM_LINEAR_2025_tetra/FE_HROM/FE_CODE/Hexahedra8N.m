function  [Ne BeXi] = Hexahedra8N(xiV) ; 
% Shape functions and derivatives for 4-node quadrilateral 
xi = xiV(1) ; eta = xiV(2) ; zeta = xiV(3) ; 
% Matrix of shape functions
Ne14 =0.125*[(1-xi)*(1-eta)*(1-zeta), (1+xi)*(1-eta)*(1-zeta), (1+xi)*(1+eta)*(1-zeta), (1-xi)*(1+eta)*(1-zeta) ];
Ne58 =0.125*[(1-xi)*(1-eta)*(1+zeta), (1+xi)*(1-eta)*(1+zeta), (1+xi)*(1+eta)*(1+zeta), (1-xi)*(1+eta)*(1+zeta) ];
Ne = [Ne14,Ne58] ; 
% Matrix of the gradient of shape functions 
BeXi14 = 0.125*[ -(1-eta)*(1-zeta),  (1-eta)*(1-zeta),  (1+eta)*(1-zeta), -(1+eta)*(1-zeta) ; 
             -(1-xi)*(1-zeta) , -(1+xi)*(1-zeta) , (1+xi)*(1-zeta), (1-xi)*(1-zeta)  ;
             -(1-xi)*(1-eta) , -(1+xi)*(1-eta) , -(1+xi)*(1+eta) , -(1-xi)*(1+eta) ] ; 
 
BeXi58 = 0.125*[ -(1-eta)*(1+zeta),  (1-eta)*(1+zeta),  (1+eta)*(1+zeta), -(1+eta)*(1+zeta) ; 
             -(1-xi)*(1+zeta) , -(1+xi)*(1+zeta) , (1+xi)*(1+zeta), (1-xi)*(1+zeta)  ;
             (1-xi)*(1-eta) , (1+xi)*(1-eta) , (1+xi)*(1+eta) , (1-xi)*(1+eta) ] ;     
         
         
BeXi = [BeXi14,BeXi58] ;    