function  [TRANSLATION, R_incre ]= RotationMatrixFrom3points(Xref,Xrot,TOL_LENGTHS)

% How to compute the rotation  ? We have three vectors in the reference configuration
% and three vectors in the rotated configuration

R_incre = [] ;
% Reference and rotated vectors on ---------------------------------
vREF = zeros(size(Xref,1),2) ;
vROT = vREF ;
% -----------------
TRANSLATION = Xrot(:,1)-Xref(:,1) ; 
% Reference
vREF(:,1) =  Xref(:,2)-Xref(:,1) ;
vREF(:,2) =  Xref(:,3)-Xref(:,1) ;
% Norm
nvREF = sqrt(sum(vREF.^2,1)) ;
% ------------------
% Rotated
vROT(:,1) =  Xrot(:,2)-Xrot(:,1) ;
vROT(:,2) =  Xrot(:,3)-Xrot(:,1) ;
% Norm
nvROT = sqrt(sum(vROT.^2,1)) ;
% ----
% Check that the length of the vectors is the same ---within
% reasonable bounds
diffLENGTHS = abs(nvROT-nvREF)./nvREF ;
if any(diffLENGTHS>TOL_LENGTHS)
   % disp('No matching')    
else    
    %  Unitary vectors
    for inodesLOC = 1:size(vROT,2)
        vROT(:,inodesLOC) = vROT(:,inodesLOC)./nvROT(inodesLOC) ;
    end
    for inodesLOC = 1:size(vREF,2)
        vREF(:,inodesLOC) = vREF(:,inodesLOC)./nvREF(inodesLOC) ;
    end
    % Determination of  rotation matrix REFERENCE configuration
    Rref = zeros(size(vREF,1)) ;
    Rref(:,1) = vREF(:,1) ;
    Rref(:,3) = cross(vREF(:,1),vREF(:,2)) ;
    Rref(:,3) = Rref(:,3)/norm(Rref(:,3)) ;
    Rref(:,2) =cross(Rref(:,3),Rref(:,1)) ;
    % Determination of  rotation matrix ROTATED configuration
    Rrot = zeros(size(vROT,1)) ;
    Rrot(:,1) = vROT(:,1) ;
    Rrot(:,3) = cross(vROT(:,1),vROT(:,2)) ;
    Rrot(:,3) = Rrot(:,3)/norm(Rrot(:,3)) ;
    Rrot(:,2) =cross(Rrot(:,3),Rrot(:,1)) ;
    
    % Total rotation vector    v_rot =   Rrot*v_0
    % Incremental    --->  v_rot = R_incre*(Rref*v_0)
    % Therefore,  Rrot = R_incre*Rref ; --> R_incre = Rrot*Rref' ;
    R_incre = Rrot*Rref' ;
    
end