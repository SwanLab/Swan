function RotMatFaceINER = AxesPrincipalInertia(RotMatFace,COORglo,GeometricMassMatrix)
% Determination principal axes of a given surfae, defined by its global coordinates (referred to the centroid, COORglo), 
% a rotation matrix RotMatFace = [n,t1,t2].  GeometricMassMatrix is the
% geometric mass matrix = int(Nshape^T,Nshape).  
% JAHO, 31-Aug-2021 



% Determinatio n
COORrelA_global =  RotMatFace'*COORglo   ;

RigidBodyModes = ConstructBasisRigidBody(COORrelA_global') ;
% GEOMETRIC PROPERTIES
nmodes = size(RigidBodyModes,2) ;
GGG = zeros(nmodes,nmodes)  ;
ndim = 3; 
for idim = 1:ndim
    GGG = GGG + (RigidBodyModes(idim:ndim:end,:)'*GeometricMassMatrix*RigidBodyModes(idim:ndim:end,:) ) ;
end

%     \item If $\G$ is not diagonal, then we can determine the principal directions as follows. Firstly, to avoid round-off errors, we symmetrize $\G$ by making $\G \leftarrow 0.5 (\G^T + \G)$.  Then we compute the eigenvalues and eigenvectors  of $\G$, such that
%    \begin{equation}
%     \B^T \G \B = \A
%    \end{equation}
% where   $\A$ is a diagonal matrix, while $\B$ is a basis change (rotation) matrix. Given the coordinates  $\x_{loc}$ of a vector expressed in the system of coordinates attached to the surface, with arbitrary tangent vectors, $\B_{princ,loc}[1,\zero; \zero \B$ maps this vector onto the principal axes:  $\x_{princ} = \B_{princ,loc}  \x_{loc}$. Recall  we seek a matrix $\R_{face}$ such that $\x_{glo} = \R_{face} \x_{princ} $; it follows thus that
%
% \begin{equation}
%  \x_{glo} = (\R_{face}^{non}) \x_{loc}  = (\R_{face}^{non}) \B_{princ,loc}^T \x_{princ}
% \end{equation}
% Thus
%
% \begin{equation}
%   \boxed{\R_{face} =  (\R_{face}^{non}) \B_{princ,loc}^T}
% \end{equation}

MomentOfInertia = GGG(5:6,5:6) ;
MomentOfInertia=  0.5*(MomentOfInertia'+MomentOfInertia);

[B,aaa,bbb] = eig(MomentOfInertia)  ;

B_princ_loc = eye(3) ;
B_princ_loc(2:3,2:3) = B ;

RotMatFaceINER = RotMatFace*B_princ_loc' ;

RotMatFaceINER(:,3) = cross(RotMatFaceINER(:,1),RotMatFaceINER(:,2)); 
