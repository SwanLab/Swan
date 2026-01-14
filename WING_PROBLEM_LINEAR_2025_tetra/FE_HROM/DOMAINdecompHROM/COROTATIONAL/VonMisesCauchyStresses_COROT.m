function [VMSTRESS,CauchyStress ]= VonMisesCauchyStresses_COROT(PoneST,FgradST,ndim,detFgrad,DATA) ;
% Adaptation of VonMisesCauchyStresses.m, for corotational approach
% ---------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

% Determinant of deformation gradient
% -----------------------------------
if   DATA.SMALL_STRAIN_KINEMATICS ==0
    if isempty(detFgrad)
        detFgrad = Determinant_Fgrad(FgradST,DATA.MESH.ndim) ;
    end
    CauchyStress = CauchyStressFromPK1(PoneST,FgradST,detFgrad,DATA.MESH.ndim) ;
else
    CauchyStress = PoneST ; 
end
% Cauchy Stresses
% Von Mises Stresses
% ------------------
ndimFINE =  DATA.MESH.ndim;  % THIS IS THE ONLY CHANGE ! 
VMSTRESS = VonMisesStressCOMP(CauchyStress,ndimFINE,DATA) ;
