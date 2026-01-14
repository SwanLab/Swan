function CauchyStress = CauchyStressFromPK1(PoneST,FgradST,detF,ndim) 
%--------------------------------------------------------------------------
% function CauchyStress = CauchyStressFromPK1(PoneST, FgradST, detF, ndim)
%
% PURPOSE:
%   Computes the Cauchy stress tensor (in Voigt notation) from the first 
%   Piola-Kirchhoff stress tensor using the standard transformation:
%
%       σ = (1 / det(F)) * F * Pᵗ
%
%   where:
%       - F  is the deformation gradient
%       - P  is the first Piola-Kirchhoff stress tensor
%       - σ  is the Cauchy stress tensor
%       - det(F) is the determinant of F
%
%   This routine handles both 2D (plane strain/stress) and 3D cases, with
%   a fully vectorized implementation per Gauss point for all elements.
%
% FORMULATION DETAILS:
%   - The transformation assumes Voigt notation:
%       In 2D: [σ_xx, σ_yy, σ_xy]
%       In 3D: [σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]
%   - This function is symbolically derived using the script:
%         /SYMBOLIC/SymCauchyStress.m
%   - The result is consistent with finite deformation theory, and the 
%     output is intended to be compatible with subsequent von Mises stress 
%     computations and postprocessing routines.
%
% INPUTS:
%   PoneST   : First Piola-Kirchhoff stress tensor (stacked vector form)
%   FgradST  : Deformation gradient at Gauss points (stacked vector form)
%   detF     : Determinant of the deformation gradient at each Gauss point
%   ndim     : Problem dimensionality (2 or 3)
%
% OUTPUT:
%   CauchyStress : Cauchy stress tensor in Voigt form, as a column vector
%                  stacked by Gauss points (size: ngaus * nstrain × 1)
%
% EXAMPLE USAGE:
%   Used during postprocessing or for constitutive updates in finite
%   deformation mechanics (nonlinear FEM). Should be preceded by a proper 
%   evaluation of F and P at Gauss points.
%
% NOTES:
%   - The stress update is fully vectorized over all Gauss points.
%   - If `detF` is not precomputed, see `Determinant_Fgrad`.
%   - For small strain cases, Cauchy = PK1 = PK2 and this transformation
%     is bypassed in the main solver (see `VonMisesCauchyStresses.m`).
%
% REFERENCES:
%   - Implementation notes: DOCS/05_IMPLEMENTATION_STATIC.pdf, page 26
%   - Symbolic derivation: SYMBOLIC/SymCauchyStress.m
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, UPC-CIMNE
%   Created: 10-Dec-2020, Barcelona
%   Comments by ChatGPT, 13-May-2025
%--------------------------------------------------------------------------

% Cauchy stresses from PK1 stress tensor ---- See
%  DOCS/05_IMPLEMENTATION_STATIC.pdf, page 26
% cauchy = 1/detF *F'*P
% See symbolic function /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/SymCauchyStress.m
% JAHO, 10-dec-2020
% --------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

if ndim == 2
    nF = 4 ;
    nelem_ngaus = length(FgradST)/nF  ;
    
    FROWS = cell(1,nF) ;
    for icols  =1:nF
        FROWS{icols} = icols:nF:size(FgradST,1) ;
    end
    nstrain = 3;
    
    CauchyStress = zeros(nelem_ngaus*nstrain,1) ;
    
    CROWS = cell(1,nstrain) ;
    for icols  =1:nstrain
        CROWS{icols} = icols:nstrain:size(CauchyStress,1) ;
    end
    
    
    % Generated automatically by  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/SymCauchyStress.m
    
    CauchyStress(CROWS{1}) = (FgradST(FROWS{1}).*PoneST(FROWS{1}) + FgradST(FROWS{4}).*PoneST(FROWS{4}))./detF;
    CauchyStress(CROWS{2}) = (FgradST(FROWS{2}).*PoneST(FROWS{2}) + FgradST(FROWS{3}).*PoneST(FROWS{3}))./detF;
    CauchyStress(CROWS{3}) = (FgradST(FROWS{1}).*PoneST(FROWS{3}) + FgradST(FROWS{4}).*PoneST(FROWS{2}))./detF;
    
    
else
    
    
    nF = 9 ;
    nelem_ngaus = length(FgradST)/nF  ;
    
    FROWS = cell(1,nF) ;
    for icols  =1:nF
        FROWS{icols} = icols:nF:size(FgradST,1) ;
    end
    nstrain = 6;
    
    CauchyStress = zeros(nelem_ngaus*nstrain,1) ;
    
    CROWS = cell(1,nstrain) ;
    for icols  =1:nstrain
        CROWS{icols} = icols:nstrain:size(CauchyStress,1) ;
    end
    % Generated automatically by  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/SymCauchyStress.m
    
    CauchyStress(CROWS{1}) = (FgradST(FROWS{1}).*PoneST(FROWS{1}) + FgradST(FROWS{8}).*PoneST(FROWS{8}) + FgradST(FROWS{9}).*PoneST(FROWS{9}))./detF;
    CauchyStress(CROWS{2}) = (FgradST(FROWS{2}).*PoneST(FROWS{2}) + FgradST(FROWS{6}).*PoneST(FROWS{6}) + FgradST(FROWS{7}).*PoneST(FROWS{7}))./detF;
    CauchyStress(CROWS{3}) = (FgradST(FROWS{3}).*PoneST(FROWS{3}) + FgradST(FROWS{4}).*PoneST(FROWS{4}) + FgradST(FROWS{5}).*PoneST(FROWS{5}))./detF;
    CauchyStress(CROWS{4}) = (FgradST(FROWS{2}).*PoneST(FROWS{4}) + FgradST(FROWS{7}).*PoneST(FROWS{3}) + FgradST(FROWS{6}).*PoneST(FROWS{5}))./detF;
    CauchyStress(CROWS{5}) = (FgradST(FROWS{1}).*PoneST(FROWS{5}) + FgradST(FROWS{8}).*PoneST(FROWS{3}) + FgradST(FROWS{9}).*PoneST(FROWS{4}))./detF;
    CauchyStress(CROWS{6}) = (FgradST(FROWS{1}).*PoneST(FROWS{6}) + FgradST(FROWS{9}).*PoneST(FROWS{2}) + FgradST(FROWS{8}).*PoneST(FROWS{7}))./detF;
    
end