function [VMSTRESS,CauchyStress ]= VonMisesCauchyStresses(PoneST,FgradST,ndim,detFgrad,DATA) 
%--------------------------------------------------------------------------
% [VMSTRESS, CauchyStress] = VonMisesCauchyStresses(PoneST, FgradST, ndim, detFgrad, DATA)
%
% PURPOSE:
%   Computes the Cauchy stress tensor from the first Piola-Kirchhoff (PK1) 
%   stress tensor, and evaluates the corresponding Von Mises equivalent 
%   stress at each integration point.
%
% DESCRIPTION:
%   For finite strains, the Cauchy stress is computed as:
%
%       σ = (1 / det(F)) * F * P
%
%   where:
%     - σ: Cauchy stress tensor
%     - P: PK1 stress tensor (input as PoneST)
%     - F: Deformation gradient
%     - det(F): Determinant of F
%
%   For small strain kinematics, the Cauchy stress is assumed to be equal
%   to the PK1 tensor (P = σ), and no transformation is performed.
%
%   After computing σ, the Von Mises stress is extracted using the helper
%   function VonMisesStressCOMP.
%
% INPUT:
%   PoneST     : First Piola-Kirchhoff stress tensor (stacked vector format)
%   FgradST    : Deformation gradient (stacked format)
%   ndim       : Number of spatial dimensions (2 or 3)
%   detFgrad   : Determinant of F, optional (can be empty)
%   DATA       : Structure containing problem data and flags, including:
%                - DATA.SMALL_STRAIN_KINEMATICS
%                - DATA.MESH.posgp
%
% OUTPUT:
%   VMSTRESS      : Vector of Von Mises stress values (1 per Gauss point)
%   CauchyStress  : Cauchy stress tensor (same size as PoneST)
%
% REFERENCE:
%   - DOCS/05_IMPLEMENTATION_STATIC.pdf, page 26
%   - J.A. Hernández Ortega, 10-Nov-2020
%   - Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

% Computation of Von Mises Stresses (Cauchy Stresses)
% See DOCS/05_IMPLEMENTATION_STATIC.pdf, page 26
% cauchy = 1/detF *F'*P
% JAHO, 10-nov-2020 
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
ndimFINE = size(DATA.MESH.posgp,1) ; 
VMSTRESS = VonMisesStressCOMP(CauchyStress,ndimFINE,DATA) ;
