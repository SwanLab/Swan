function CECMoutput = CECM_internal_forces(OPERFE,PhiDEF,BasisPone,DATA,MESH,DATA_ECM)
%--------------------------------------------------------------------------
% function CECMoutput = CECM_internal_forces(OPERFE, PhiDEF, BasisPone, DATA, MESH, DATA_ECM)
%
% PURPOSE:
%   Computes the reduced internal force vector using the Continuous Empirical
%   Cubature Method (CECM), based on a given reduced deformation (PhiDEF).
%
%   This function acts as a key online step in projection-based nonlinear
%   model order reduction, where the internal forces are approximated via 
%   a reduced stress basis (`BasisPone`) and empirical cubature points 
%   derived from the full-order integration domain.
%
% INPUTS:
% -------
%   - OPERFE     : Structure with offline operators (Bst, Nst, posgp, wSTs)
%                  containing the reduced strain-displacement operators and
%                  quadrature data.
%
%   - PhiDEF     : Reduced deformation vector (rDOFs x 1)
%
%   - BasisPone  : Matrix whose columns are the reduced stress modes used to
%                  reconstruct the internal forces in a Galerkin projection.
%                  It typically comes from snapshot-based SVD offline.
%
%   - DATA       : Structure with global data and optional parameters.
%                  Includes options for scaling, tolerances, etc.
%
%   - MESH       : Mesh structure with fields such as COOR, posgp, ngausE.
%
%   - DATA_ECM   : Configuration parameters for cubature generation. Includes
%                  tolerance, max points, and flags for using DECM points.
%
% OUTPUT:
% -------
%   - CECMoutput : Structure with the result of the cubature procedure,
%                  including selected integration points and reconstructed 
%                  internal force vector.
%
% METHOD OVERVIEW:
% ----------------
%   1. Compute reduced stress tensor:      `BstRED = Bst * PhiDEF`
%   2. Assemble internal force integrands: `A = BasisF_from_BasisStress_PK1(...)`
%   3. Interpolate node positions in deformed configuration: `xFE`
%   4. Call `ContinuousECMgen2023` to generate cubature points for reduced integration.
%   5. Store resulting indices and optionally compute a posteriori errors.
%
% NOTES:
% ------
%   - The integrand matrix `A` is of size [ngaus × r] and defines the contributions
%     of the stress basis to the internal forces at Gauss points.
%
%   - If `DATA_ECM.UseDECMpoints == 0`, an error indicator is computed via `ErrorCalcLocal2023`
%
%   - The `CECMoutput.DECM_indexes_points` field contains the selected integration
%     points to be used in online evaluations.
%
% DEPENDENCIES:
%   - BasisF_from_BasisStress_PK1
%   - ContinuousECMgen2023
%   - ErrorCalcLocal2023
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, 2023
%--------------------------------------------------------------------------


if nargin == 0
    load('tmp.mat')
end

% 
% %if DATAoffline.IncludeSingularValueBasisStressesPone == 1
%     BasisPone = bsxfun(@times,BasisPone',DATA.BasisSTRESS_SINGULAR_VALUES)' ;
% %end

BstRED= OPERFE.Bst*PhiDEF;
A = BasisF_from_BasisStress_PK1(BstRED,BasisPone,DATA)  ;

xNODES = MESH.COOR';
xFE = OPERFE.Nst*xNODES(:) ;

ndim = size(MESH.COOR,2) ;
xFE = reshape(xFE,ndim,[])' ;
%xFE = xFE(1:ndim:end,1) ;
%MESH.COOR = MESH.COOR(:,1)  ;
wFE = OPERFE.wSTs ;
MESH.posgp = OPERFE.posgp ; % 10-May-2023



delete('CECMoutputF.txt')
diary('CECMoutputF.txt')

% --------------------------------
% Evaluation integrands
% --------------------------------
 
DATAloc.ExactIntegral = A'*wFE ;
MESH.ngausE = DATA.MESH.ngaus_RHS ;

% Continuous Empirical Cubature Method
[CECMoutput,DATA_AUX]= ContinuousECMgen2023(A,xFE,wFE,DATAloc,MESH,DATA_ECM) ;

CECMoutput.DECM_indexes_points   = DATA_AUX.indexPoints_DECM ; 
 


if DATA_ECM.UseDECMpoints ==0 
  ErrorCalcLocal2023(DATA,CECMoutput,DATA_AUX,A,DATAloc)  ;
end