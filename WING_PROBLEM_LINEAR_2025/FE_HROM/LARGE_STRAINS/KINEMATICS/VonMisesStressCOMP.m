function VMSTRESS = VonMisesStressCOMP(CauchyStress,ndim,DATA)
% Computation of Von Mises stress, given the stacked vector of Cauchy
% Stresses
% See, DOCS/05_IMPLEMENTATION_STATIC.pdf, page 27
% JAHO, 10-Dec-2020
% Comments by ChatGPT-4 (12-May-2025)
%VonMisesStressCOMP Compute the von Mises equivalent stress from Cauchy stress tensors
%
%   This function computes the von Mises stress for 2D or 3D problems given the
%   stacked Cauchy stress tensor at Gauss points, using the standard definition
%   based on deviatoric stress invariants. It supports both plane strain and
%   3D stress representations with either 3, 4, or 6 strain components.
%
%   INPUT:
%     CauchyStress - Matrix of size (nstrain × ntimes) containing stress components
%                    in stacked vector form per Gauss point (column-wise per time step)
%     ndim         - Spatial dimension of the problem (2 or 3)
%     DATA         - Data structure that includes:
%                      - DATA.MESH.nstrain: number of stress components per Gauss point
%
%   OUTPUT:
%     VMSTRESS     - Matrix of von Mises stress values (size: ngaus × ntimes)
%
%   FORMULAS USED:
%   -------------------------------------------------------------------------
%     For 2D, nstrain = 3 (σxx, σyy, σxy):
%         σ_vm = sqrt(σxx² + σyy² − σxx·σyy + 3·σxy²)
%
%     For 2D, nstrain = 4 (σxx, σyy, σxy, σzz) → plane strain:
%         σ_vm = sqrt(½[(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)²] + 3·σxy²)
%
%     For 3D, nstrain = 6 (σxx, σyy, σzz, σxy, σxz, σyz):
%         σ_vm = sqrt(½[(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)² + 
%                       6(σxy² + σxz² + σyz²)])
%
%   NOTES:
%     - The stress input is assumed to be in Voigt notation and stacked column-wise.
%     - The function automatically handles the appropriate formula based on `ndim` and `nstrain`.
%
%   REFERENCES:
%     - See DOCS/05_IMPLEMENTATION_STATIC.pdf, page 27
%     - Classical elasticity and plasticity textbooks on von Mises yield criteria
%
%   AUTHOR:
%     Joaquín A. Hernández Ortega (JAHO), UPC BarcelonaTech, 10-Dec-2020
%
%   SEE ALSO:
%     StrainGreenLagrange, PK2stress_Constitutive_ModelVAR

% --------------------------------------------------------------------
if nargin == 0
    load('tmp3.mat')
end


if ndim == 2
    nstrain =  DATA.MESH.nstrain ;
    %  nelem_ngaus = length(CauchyStress)/nstrain  ;
    CROWS = cell(1,nstrain) ;
    for icols  =1:nstrain
        CROWS{icols} = icols:nstrain:size(CauchyStress,1) ;
    end
    
    if nstrain == 3
        VMSTRESS =   (CauchyStress(CROWS{1},:) - CauchyStress(CROWS{2},:) ).^2   ;
        VMSTRESS = VMSTRESS +  CauchyStress(CROWS{2},:).^2   ;
        VMSTRESS = VMSTRESS +  CauchyStress(CROWS{1},:).^2   ;
        VMSTRESS = VMSTRESS +  6*CauchyStress(CROWS{3},:).^2   ;
    elseif nstrain == 4
        VMSTRESS =   (CauchyStress(CROWS{1},:) - CauchyStress(CROWS{2},:)).^2   ;
        VMSTRESS = VMSTRESS +  (CauchyStress(CROWS{2},:) - CauchyStress(CROWS{4},:) ).^2   ;
        VMSTRESS = VMSTRESS +  (CauchyStress(CROWS{1},:) - CauchyStress(CROWS{4},:) ).^2   ;
        VMSTRESS = VMSTRESS +  6*CauchyStress(CROWS{3},:).^2   ;
        
    else
        error('Erroneous option')
    end
    
    
    
    
    VMSTRESS = sqrt(0.5*VMSTRESS) ;
    
    
else
    nstrain = 6 ;
 %   nelem_ngaus = length(CauchyStress)/nstrain  ;
    
    
    CROWS = cell(1,nstrain) ;
    for icols  =1:nstrain
        CROWS{icols} = icols:nstrain:size(CauchyStress,1) ;
    end
    
    VMSTRESS =  (CauchyStress(CROWS{1},:) - CauchyStress(CROWS{2},:) ).^2   ;
    VMSTRESS = VMSTRESS +  (CauchyStress(CROWS{2},:) - CauchyStress(CROWS{3},:) ).^2   ;
    VMSTRESS = VMSTRESS +  (CauchyStress(CROWS{3},:) - CauchyStress(CROWS{1},:) ).^2   ;
    VMSTRESS = VMSTRESS +  6*(CauchyStress(CROWS{4},:).^2 + CauchyStress(CROWS{5},:).^2 +  CauchyStress(CROWS{6},:).^2   )  ;
    
    VMSTRESS = sqrt(0.5*VMSTRESS) ;
    
end
