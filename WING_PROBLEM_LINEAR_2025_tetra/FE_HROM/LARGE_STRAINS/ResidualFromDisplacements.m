function [GLSTRAINS,PK2STRESS,celastST,RESID,Fint,FgradST,PoneST,detFgrad] = ResidualFromDisplacements(OPERFE,d,MATPRO,DATA,FEXT)
% 2. Deformation gradient at all Gauss points
if nargin == 0
    load('tmp.mat')
end

% Stresses from displacements 
[GLSTRAINS,PK2STRESS,celastST,FgradST,PoneST,detFgrad] = StressesFromDisplacements(OPERFE,d,MATPRO,DATA) ; 

if isempty(PK2STRESS)
    PoneST = [] ; Fint = [] ; RESID = [] ;
else
   
    % 6.1. Internal forces
    Fint = InternalForces(OPERFE,PoneST,PK2STRESS,DATA) ;
    
    % 6.2. Residual
    RESID  = Fint- FEXT;
end



% ndof = size(OPERFE.Bst,2) ;
% FgradST =  OPERFE.Bst*d(1:ndof,:)  ; %+ OPERFE.IDENTITY_F ;
% if DATA.SMALL_STRAIN_KINEMATICS ==0
%     for idim = 1:DATA.MESH.ndim
%         LOCROWS = idim:DATA.MESH.ndim^2:length(FgradST) ;
%         FgradST(LOCROWS,:) = 1+FgradST(LOCROWS,:) ;
%     end
%     % 3. Green-Lagrante strains at all Gauss points
%     GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
% else
%     % In small strain kinematics, FgradST == gradU (because we haven't added the identity matrix )
%     [GLSTRAINS,FgradST] = StrainGreenLagrange_small(FgradST,DATA.MESH.ndim) ;   
% end
% % 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
% [PK2STRESS,celastST,detFgrad ]= PK2stress_Constitutive_Model(GLSTRAINS,MATPRO,DATA,FgradST) ;
% 
% if isempty(PK2STRESS)
%     PoneST = [] ; Fint = [] ; RESID = [] ;
% else
%        % 5. 1st Piola-Kirchhoff stresses at all Gauss Points
%     PoneST = PK1stress(PK2STRESS,FgradST,DATA.MESH.ndim) ;
%     % 6.1. Internal forces
%     Fint = InternalForces(OPERFE,PoneST,DATA) ;
%     
%     % 6.2. Residual
%     RESID  = Fint- FEXT;
% end


