function celasLARGEmat =CelasLARGEmat_allgauss(celastST,FgradST,ndim) 
% ---------------------------------------------------------------------------------------------------
% FUNCTION: CelasLARGEmat_allgauss
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   Computes the **material tangent stiffness contribution** (also called material or constitutive
%   "celas" matrix) to the total tangent matrix used in large strain finite element analysis.
%   The result accounts for the **nonlinear material response** at each Gauss point using the
%   deformation gradient F and the material tangent operator C.
%
%   This routine returns the block-diagonal structure of the material part of the stiffness matrix
%   in a format compatible with later assembly (e.g., via `ConvertBlockDiag`).
%
% USAGE:
%   celasLARGEmat = CelasLARGEmat_allgauss(celastST,FgradST,ndim)
%
% INPUTS:
%   - celastST  : Elastic (material) tangent tensor at each Gauss point (Voigt form),
%                 size = [nstrain x ngauss]
%   - FgradST   : Deformation gradient (Voigt-like form), size = [nF x ngauss]
%   - ndim      : Spatial dimension (2 or 3)
%
% OUTPUTS:
%   - celasLARGEmat : Matrix containing all Gauss-point contributions to the material part of the
%                     tangent matrix, stored in unassembled form (prior to B^T * celas * B).
%
% IMPLEMENTATION DETAILS:
%   - In **2D**, Voigt notations are used:
%         FgradST rows: [F11; F22; F12; F21]
%         celastST rows: [C11; C22; C12]
%   - In **3D**, Voigt notations:
%         FgradST rows: [F11; F22; F33; F12; F21; F13; F31; F23; F32]
%         celastST rows: [C11; C22; C33; C12; C13; C23]
%
%   - The constitutive update uses analytical formulas derived via symbolic manipulation
%     (see `Kmaterial_T.m`) and uses full non-symmetric formulation to accommodate
%     large strain kinematics.
%
%   - Symmetry of celas is enforced manually for performance and compactness.
%
% CONTEXT:
%   This function is a core component of the total tangent matrix assembly in large deformation
%   settings, e.g., for:
%     - Neo-Hookean or other hyperelastic models
%     - Finite element solvers using Updated or Total Lagrangian frameworks
%     - EIFEM and HROM approaches with geometric nonlinearity
%
% REFERENCES:
%   - "README_RigidBodyMotions.pdf", page 19
%   - Symbolic derivations: `Kmaterial_T.m`
%   - Application context: EIFEM nonlinear ROM framework (UPC/CIMNE)
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Date: 06-Jan-2024, Barcelona
%   Comments by ChatGPT4, 12-May-2025
% SEE ALSO:
%   - KstiffLargeStrains.m
%   - ConvertBlockDiag.m
%   - CelasLARGEgeo_allgauss.m
%
% ---------------------------------------------------------------------------------------------------

% Assembly material celasLARGE Matrix
% See % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/
% README_RigidBodyMotions.pdf, page 19
% See also /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/Kmaterial_T.m
%
if nargin == 0
    load('tmp.mat')
end





if ndim == 2
    nF = 4 ;
    nstrain = 3;
    nelem_ngaus = length(FgradST)/nF  ;
    celasLARGEmat = zeros(nF*nelem_ngaus,nF)  ;
    
    FROWS = cell(1,nF) ;
    for icols  =1:nF
        FROWS{icols} = icols:nF:size(FgradST,1) ;
    end
    
    CROWS = cell(1,nstrain) ;
    for icols  =1:nstrain
        CROWS{icols} = icols:nstrain:size(celastST,1) ;
    end
    
    % Generated automatically by /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/SYMBOLIC/Kmaterial_T.m
    celasLARGEmat(FROWS{1},1) = FgradST(FROWS{1}).*(celastST(CROWS{1},1).*FgradST(FROWS{1}) + celastST(CROWS{3},1).*FgradST(FROWS{3})) + FgradST(FROWS{3}).*(celastST(CROWS{1},3).*FgradST(FROWS{1}) + celastST(CROWS{3},3).*FgradST(FROWS{3}));
    celasLARGEmat(FROWS{1},2) = FgradST(FROWS{2}).*(celastST(CROWS{1},2).*FgradST(FROWS{1}) + celastST(CROWS{3},2).*FgradST(FROWS{3})) + FgradST(FROWS{4}).*(celastST(CROWS{1},3).*FgradST(FROWS{1}) + celastST(CROWS{3},3).*FgradST(FROWS{3}));
    celasLARGEmat(FROWS{1},3) = FgradST(FROWS{1}).*(celastST(CROWS{1},3).*FgradST(FROWS{1}) + celastST(CROWS{3},3).*FgradST(FROWS{3})) + FgradST(FROWS{3}).*(celastST(CROWS{1},2).*FgradST(FROWS{1}) + celastST(CROWS{3},2).*FgradST(FROWS{3}));
    celasLARGEmat(FROWS{1},4) = FgradST(FROWS{4}).*(celastST(CROWS{1},1).*FgradST(FROWS{1}) + celastST(CROWS{3},1).*FgradST(FROWS{3})) + FgradST(FROWS{2}).*(celastST(CROWS{1},3).*FgradST(FROWS{1}) + celastST(CROWS{3},3).*FgradST(FROWS{3}));
    celasLARGEmat(FROWS{2},2) = FgradST(FROWS{2}).*(celastST(CROWS{2},2).*FgradST(FROWS{2}) + celastST(CROWS{3},2).*FgradST(FROWS{4})) + FgradST(FROWS{4}).*(celastST(CROWS{2},3).*FgradST(FROWS{2}) + celastST(CROWS{3},3).*FgradST(FROWS{4}));
    celasLARGEmat(FROWS{2},3) = FgradST(FROWS{1}).*(celastST(CROWS{2},3).*FgradST(FROWS{2}) + celastST(CROWS{3},3).*FgradST(FROWS{4})) + FgradST(FROWS{3}).*(celastST(CROWS{2},2).*FgradST(FROWS{2}) + celastST(CROWS{3},2).*FgradST(FROWS{4}));
    celasLARGEmat(FROWS{2},4) = FgradST(FROWS{4}).*(celastST(CROWS{2},1).*FgradST(FROWS{2}) + celastST(CROWS{3},1).*FgradST(FROWS{4})) + FgradST(FROWS{2}).*(celastST(CROWS{2},3).*FgradST(FROWS{2}) + celastST(CROWS{3},3).*FgradST(FROWS{4}));
    celasLARGEmat(FROWS{3},3) = FgradST(FROWS{1}).*(celastST(CROWS{2},3).*FgradST(FROWS{3}) + celastST(CROWS{3},3).*FgradST(FROWS{1})) + FgradST(FROWS{3}).*(celastST(CROWS{2},2).*FgradST(FROWS{3}) + celastST(CROWS{3},2).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{3},4) = FgradST(FROWS{4}).*(celastST(CROWS{2},1).*FgradST(FROWS{3}) + celastST(CROWS{3},1).*FgradST(FROWS{1})) + FgradST(FROWS{2}).*(celastST(CROWS{2},3).*FgradST(FROWS{3}) + celastST(CROWS{3},3).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{4},4) = FgradST(FROWS{4}).*(celastST(CROWS{1},1).*FgradST(FROWS{4}) + celastST(CROWS{3},1).*FgradST(FROWS{2})) + FgradST(FROWS{2}).*(celastST(CROWS{1},3).*FgradST(FROWS{4}) + celastST(CROWS{3},3).*FgradST(FROWS{2}));
    celasLARGEmat(FROWS{2},1) = celasLARGEmat(FROWS{1},2)  ;
    celasLARGEmat(FROWS{3},1) = celasLARGEmat(FROWS{1},3)  ;
    celasLARGEmat(FROWS{3},2) = celasLARGEmat(FROWS{2},3)  ;
    celasLARGEmat(FROWS{4},1) = celasLARGEmat(FROWS{1},4)  ;
    celasLARGEmat(FROWS{4},2) = celasLARGEmat(FROWS{2},4)  ;
    celasLARGEmat(FROWS{4},3) = celasLARGEmat(FROWS{3},4)  ;
    
    
    
    
else
    nF = 9 ;
    nstrain = 6;
    nelem_ngaus = length(FgradST)/nF  ;
    celasLARGEmat = zeros(nF*nelem_ngaus,nF)  ;
    
    FROWS = cell(1,nF) ;
    for icols  =1:nF
        FROWS{icols} = icols:nF:size(FgradST,1) ;
    end
    
    CROWS = cell(1,nstrain) ;
    for icols  =1:nstrain
        CROWS{icols} = icols:nstrain:size(celastST,1) ;
    end
    
    
    celasLARGEmat(FROWS{1},1) = FgradST(FROWS{1}).*(celastST(CROWS{1},1).*FgradST(FROWS{1}) + celastST(CROWS{5},1).*FgradST(FROWS{5}) + celastST(CROWS{6},1).*FgradST(FROWS{6})) + FgradST(FROWS{5}).*(celastST(CROWS{1},5).*FgradST(FROWS{1}) + celastST(CROWS{5},5).*FgradST(FROWS{5}) + celastST(CROWS{6},5).*FgradST(FROWS{6})) + FgradST(FROWS{6}).*(celastST(CROWS{1},6).*FgradST(FROWS{1}) + celastST(CROWS{5},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{1},2) = FgradST(FROWS{2}).*(celastST(CROWS{1},2).*FgradST(FROWS{1}) + celastST(CROWS{5},2).*FgradST(FROWS{5}) + celastST(CROWS{6},2).*FgradST(FROWS{6})) + FgradST(FROWS{4}).*(celastST(CROWS{1},4).*FgradST(FROWS{1}) + celastST(CROWS{5},4).*FgradST(FROWS{5}) + celastST(CROWS{6},4).*FgradST(FROWS{6})) + FgradST(FROWS{9}).*(celastST(CROWS{1},6).*FgradST(FROWS{1}) + celastST(CROWS{5},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{1},3) = FgradST(FROWS{3}).*(celastST(CROWS{1},3).*FgradST(FROWS{1}) + celastST(CROWS{5},3).*FgradST(FROWS{5}) + celastST(CROWS{6},3).*FgradST(FROWS{6})) + FgradST(FROWS{7}).*(celastST(CROWS{1},4).*FgradST(FROWS{1}) + celastST(CROWS{5},4).*FgradST(FROWS{5}) + celastST(CROWS{6},4).*FgradST(FROWS{6})) + FgradST(FROWS{8}).*(celastST(CROWS{1},5).*FgradST(FROWS{1}) + celastST(CROWS{5},5).*FgradST(FROWS{5}) + celastST(CROWS{6},5).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{1},4) = FgradST(FROWS{4}).*(celastST(CROWS{1},3).*FgradST(FROWS{1}) + celastST(CROWS{5},3).*FgradST(FROWS{5}) + celastST(CROWS{6},3).*FgradST(FROWS{6})) + FgradST(FROWS{2}).*(celastST(CROWS{1},4).*FgradST(FROWS{1}) + celastST(CROWS{5},4).*FgradST(FROWS{5}) + celastST(CROWS{6},4).*FgradST(FROWS{6})) + FgradST(FROWS{9}).*(celastST(CROWS{1},5).*FgradST(FROWS{1}) + celastST(CROWS{5},5).*FgradST(FROWS{5}) + celastST(CROWS{6},5).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{1},5) = FgradST(FROWS{5}).*(celastST(CROWS{1},3).*FgradST(FROWS{1}) + celastST(CROWS{5},3).*FgradST(FROWS{5}) + celastST(CROWS{6},3).*FgradST(FROWS{6})) + FgradST(FROWS{1}).*(celastST(CROWS{1},5).*FgradST(FROWS{1}) + celastST(CROWS{5},5).*FgradST(FROWS{5}) + celastST(CROWS{6},5).*FgradST(FROWS{6})) + FgradST(FROWS{6}).*(celastST(CROWS{1},4).*FgradST(FROWS{1}) + celastST(CROWS{5},4).*FgradST(FROWS{5}) + celastST(CROWS{6},4).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{1},6) = FgradST(FROWS{6}).*(celastST(CROWS{1},2).*FgradST(FROWS{1}) + celastST(CROWS{5},2).*FgradST(FROWS{5}) + celastST(CROWS{6},2).*FgradST(FROWS{6})) + FgradST(FROWS{5}).*(celastST(CROWS{1},4).*FgradST(FROWS{1}) + celastST(CROWS{5},4).*FgradST(FROWS{5}) + celastST(CROWS{6},4).*FgradST(FROWS{6})) + FgradST(FROWS{1}).*(celastST(CROWS{1},6).*FgradST(FROWS{1}) + celastST(CROWS{5},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{1},7) = FgradST(FROWS{7}).*(celastST(CROWS{1},2).*FgradST(FROWS{1}) + celastST(CROWS{5},2).*FgradST(FROWS{5}) + celastST(CROWS{6},2).*FgradST(FROWS{6})) + FgradST(FROWS{3}).*(celastST(CROWS{1},4).*FgradST(FROWS{1}) + celastST(CROWS{5},4).*FgradST(FROWS{5}) + celastST(CROWS{6},4).*FgradST(FROWS{6})) + FgradST(FROWS{8}).*(celastST(CROWS{1},6).*FgradST(FROWS{1}) + celastST(CROWS{5},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{1},8) = FgradST(FROWS{8}).*(celastST(CROWS{1},1).*FgradST(FROWS{1}) + celastST(CROWS{5},1).*FgradST(FROWS{5}) + celastST(CROWS{6},1).*FgradST(FROWS{6})) + FgradST(FROWS{3}).*(celastST(CROWS{1},5).*FgradST(FROWS{1}) + celastST(CROWS{5},5).*FgradST(FROWS{5}) + celastST(CROWS{6},5).*FgradST(FROWS{6})) + FgradST(FROWS{7}).*(celastST(CROWS{1},6).*FgradST(FROWS{1}) + celastST(CROWS{5},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{1},9) = FgradST(FROWS{9}).*(celastST(CROWS{1},1).*FgradST(FROWS{1}) + celastST(CROWS{5},1).*FgradST(FROWS{5}) + celastST(CROWS{6},1).*FgradST(FROWS{6})) + FgradST(FROWS{4}).*(celastST(CROWS{1},5).*FgradST(FROWS{1}) + celastST(CROWS{5},5).*FgradST(FROWS{5}) + celastST(CROWS{6},5).*FgradST(FROWS{6})) + FgradST(FROWS{2}).*(celastST(CROWS{1},6).*FgradST(FROWS{1}) + celastST(CROWS{5},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{6}));
    celasLARGEmat(FROWS{2},2) = FgradST(FROWS{2}).*(celastST(CROWS{2},2).*FgradST(FROWS{2}) + celastST(CROWS{4},2).*FgradST(FROWS{4}) + celastST(CROWS{6},2).*FgradST(FROWS{9})) + FgradST(FROWS{4}).*(celastST(CROWS{2},4).*FgradST(FROWS{2}) + celastST(CROWS{4},4).*FgradST(FROWS{4}) + celastST(CROWS{6},4).*FgradST(FROWS{9})) + FgradST(FROWS{9}).*(celastST(CROWS{2},6).*FgradST(FROWS{2}) + celastST(CROWS{4},6).*FgradST(FROWS{4}) + celastST(CROWS{6},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{2},3) = FgradST(FROWS{3}).*(celastST(CROWS{2},3).*FgradST(FROWS{2}) + celastST(CROWS{4},3).*FgradST(FROWS{4}) + celastST(CROWS{6},3).*FgradST(FROWS{9})) + FgradST(FROWS{7}).*(celastST(CROWS{2},4).*FgradST(FROWS{2}) + celastST(CROWS{4},4).*FgradST(FROWS{4}) + celastST(CROWS{6},4).*FgradST(FROWS{9})) + FgradST(FROWS{8}).*(celastST(CROWS{2},5).*FgradST(FROWS{2}) + celastST(CROWS{4},5).*FgradST(FROWS{4}) + celastST(CROWS{6},5).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{2},4) = FgradST(FROWS{4}).*(celastST(CROWS{2},3).*FgradST(FROWS{2}) + celastST(CROWS{4},3).*FgradST(FROWS{4}) + celastST(CROWS{6},3).*FgradST(FROWS{9})) + FgradST(FROWS{2}).*(celastST(CROWS{2},4).*FgradST(FROWS{2}) + celastST(CROWS{4},4).*FgradST(FROWS{4}) + celastST(CROWS{6},4).*FgradST(FROWS{9})) + FgradST(FROWS{9}).*(celastST(CROWS{2},5).*FgradST(FROWS{2}) + celastST(CROWS{4},5).*FgradST(FROWS{4}) + celastST(CROWS{6},5).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{2},5) = FgradST(FROWS{5}).*(celastST(CROWS{2},3).*FgradST(FROWS{2}) + celastST(CROWS{4},3).*FgradST(FROWS{4}) + celastST(CROWS{6},3).*FgradST(FROWS{9})) + FgradST(FROWS{1}).*(celastST(CROWS{2},5).*FgradST(FROWS{2}) + celastST(CROWS{4},5).*FgradST(FROWS{4}) + celastST(CROWS{6},5).*FgradST(FROWS{9})) + FgradST(FROWS{6}).*(celastST(CROWS{2},4).*FgradST(FROWS{2}) + celastST(CROWS{4},4).*FgradST(FROWS{4}) + celastST(CROWS{6},4).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{2},6) = FgradST(FROWS{6}).*(celastST(CROWS{2},2).*FgradST(FROWS{2}) + celastST(CROWS{4},2).*FgradST(FROWS{4}) + celastST(CROWS{6},2).*FgradST(FROWS{9})) + FgradST(FROWS{5}).*(celastST(CROWS{2},4).*FgradST(FROWS{2}) + celastST(CROWS{4},4).*FgradST(FROWS{4}) + celastST(CROWS{6},4).*FgradST(FROWS{9})) + FgradST(FROWS{1}).*(celastST(CROWS{2},6).*FgradST(FROWS{2}) + celastST(CROWS{4},6).*FgradST(FROWS{4}) + celastST(CROWS{6},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{2},7) = FgradST(FROWS{7}).*(celastST(CROWS{2},2).*FgradST(FROWS{2}) + celastST(CROWS{4},2).*FgradST(FROWS{4}) + celastST(CROWS{6},2).*FgradST(FROWS{9})) + FgradST(FROWS{3}).*(celastST(CROWS{2},4).*FgradST(FROWS{2}) + celastST(CROWS{4},4).*FgradST(FROWS{4}) + celastST(CROWS{6},4).*FgradST(FROWS{9})) + FgradST(FROWS{8}).*(celastST(CROWS{2},6).*FgradST(FROWS{2}) + celastST(CROWS{4},6).*FgradST(FROWS{4}) + celastST(CROWS{6},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{2},8) = FgradST(FROWS{8}).*(celastST(CROWS{2},1).*FgradST(FROWS{2}) + celastST(CROWS{4},1).*FgradST(FROWS{4}) + celastST(CROWS{6},1).*FgradST(FROWS{9})) + FgradST(FROWS{3}).*(celastST(CROWS{2},5).*FgradST(FROWS{2}) + celastST(CROWS{4},5).*FgradST(FROWS{4}) + celastST(CROWS{6},5).*FgradST(FROWS{9})) + FgradST(FROWS{7}).*(celastST(CROWS{2},6).*FgradST(FROWS{2}) + celastST(CROWS{4},6).*FgradST(FROWS{4}) + celastST(CROWS{6},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{2},9) = FgradST(FROWS{9}).*(celastST(CROWS{2},1).*FgradST(FROWS{2}) + celastST(CROWS{4},1).*FgradST(FROWS{4}) + celastST(CROWS{6},1).*FgradST(FROWS{9})) + FgradST(FROWS{4}).*(celastST(CROWS{2},5).*FgradST(FROWS{2}) + celastST(CROWS{4},5).*FgradST(FROWS{4}) + celastST(CROWS{6},5).*FgradST(FROWS{9})) + FgradST(FROWS{2}).*(celastST(CROWS{2},6).*FgradST(FROWS{2}) + celastST(CROWS{4},6).*FgradST(FROWS{4}) + celastST(CROWS{6},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{3},3) = FgradST(FROWS{3}).*(celastST(CROWS{3},3).*FgradST(FROWS{3}) + celastST(CROWS{4},3).*FgradST(FROWS{7}) + celastST(CROWS{5},3).*FgradST(FROWS{8})) + FgradST(FROWS{7}).*(celastST(CROWS{3},4).*FgradST(FROWS{3}) + celastST(CROWS{4},4).*FgradST(FROWS{7}) + celastST(CROWS{5},4).*FgradST(FROWS{8})) + FgradST(FROWS{8}).*(celastST(CROWS{3},5).*FgradST(FROWS{3}) + celastST(CROWS{4},5).*FgradST(FROWS{7}) + celastST(CROWS{5},5).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{3},4) = FgradST(FROWS{4}).*(celastST(CROWS{3},3).*FgradST(FROWS{3}) + celastST(CROWS{4},3).*FgradST(FROWS{7}) + celastST(CROWS{5},3).*FgradST(FROWS{8})) + FgradST(FROWS{2}).*(celastST(CROWS{3},4).*FgradST(FROWS{3}) + celastST(CROWS{4},4).*FgradST(FROWS{7}) + celastST(CROWS{5},4).*FgradST(FROWS{8})) + FgradST(FROWS{9}).*(celastST(CROWS{3},5).*FgradST(FROWS{3}) + celastST(CROWS{4},5).*FgradST(FROWS{7}) + celastST(CROWS{5},5).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{3},5) = FgradST(FROWS{5}).*(celastST(CROWS{3},3).*FgradST(FROWS{3}) + celastST(CROWS{4},3).*FgradST(FROWS{7}) + celastST(CROWS{5},3).*FgradST(FROWS{8})) + FgradST(FROWS{1}).*(celastST(CROWS{3},5).*FgradST(FROWS{3}) + celastST(CROWS{4},5).*FgradST(FROWS{7}) + celastST(CROWS{5},5).*FgradST(FROWS{8})) + FgradST(FROWS{6}).*(celastST(CROWS{3},4).*FgradST(FROWS{3}) + celastST(CROWS{4},4).*FgradST(FROWS{7}) + celastST(CROWS{5},4).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{3},6) = FgradST(FROWS{6}).*(celastST(CROWS{3},2).*FgradST(FROWS{3}) + celastST(CROWS{4},2).*FgradST(FROWS{7}) + celastST(CROWS{5},2).*FgradST(FROWS{8})) + FgradST(FROWS{5}).*(celastST(CROWS{3},4).*FgradST(FROWS{3}) + celastST(CROWS{4},4).*FgradST(FROWS{7}) + celastST(CROWS{5},4).*FgradST(FROWS{8})) + FgradST(FROWS{1}).*(celastST(CROWS{3},6).*FgradST(FROWS{3}) + celastST(CROWS{4},6).*FgradST(FROWS{7}) + celastST(CROWS{5},6).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{3},7) = FgradST(FROWS{7}).*(celastST(CROWS{3},2).*FgradST(FROWS{3}) + celastST(CROWS{4},2).*FgradST(FROWS{7}) + celastST(CROWS{5},2).*FgradST(FROWS{8})) + FgradST(FROWS{3}).*(celastST(CROWS{3},4).*FgradST(FROWS{3}) + celastST(CROWS{4},4).*FgradST(FROWS{7}) + celastST(CROWS{5},4).*FgradST(FROWS{8})) + FgradST(FROWS{8}).*(celastST(CROWS{3},6).*FgradST(FROWS{3}) + celastST(CROWS{4},6).*FgradST(FROWS{7}) + celastST(CROWS{5},6).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{3},8) = FgradST(FROWS{8}).*(celastST(CROWS{3},1).*FgradST(FROWS{3}) + celastST(CROWS{4},1).*FgradST(FROWS{7}) + celastST(CROWS{5},1).*FgradST(FROWS{8})) + FgradST(FROWS{3}).*(celastST(CROWS{3},5).*FgradST(FROWS{3}) + celastST(CROWS{4},5).*FgradST(FROWS{7}) + celastST(CROWS{5},5).*FgradST(FROWS{8})) + FgradST(FROWS{7}).*(celastST(CROWS{3},6).*FgradST(FROWS{3}) + celastST(CROWS{4},6).*FgradST(FROWS{7}) + celastST(CROWS{5},6).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{3},9) = FgradST(FROWS{9}).*(celastST(CROWS{3},1).*FgradST(FROWS{3}) + celastST(CROWS{4},1).*FgradST(FROWS{7}) + celastST(CROWS{5},1).*FgradST(FROWS{8})) + FgradST(FROWS{4}).*(celastST(CROWS{3},5).*FgradST(FROWS{3}) + celastST(CROWS{4},5).*FgradST(FROWS{7}) + celastST(CROWS{5},5).*FgradST(FROWS{8})) + FgradST(FROWS{2}).*(celastST(CROWS{3},6).*FgradST(FROWS{3}) + celastST(CROWS{4},6).*FgradST(FROWS{7}) + celastST(CROWS{5},6).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{4},4) = FgradST(FROWS{4}).*(celastST(CROWS{3},3).*FgradST(FROWS{4}) + celastST(CROWS{4},3).*FgradST(FROWS{2}) + celastST(CROWS{5},3).*FgradST(FROWS{9})) + FgradST(FROWS{2}).*(celastST(CROWS{3},4).*FgradST(FROWS{4}) + celastST(CROWS{4},4).*FgradST(FROWS{2}) + celastST(CROWS{5},4).*FgradST(FROWS{9})) + FgradST(FROWS{9}).*(celastST(CROWS{3},5).*FgradST(FROWS{4}) + celastST(CROWS{4},5).*FgradST(FROWS{2}) + celastST(CROWS{5},5).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{4},5) = FgradST(FROWS{5}).*(celastST(CROWS{3},3).*FgradST(FROWS{4}) + celastST(CROWS{4},3).*FgradST(FROWS{2}) + celastST(CROWS{5},3).*FgradST(FROWS{9})) + FgradST(FROWS{1}).*(celastST(CROWS{3},5).*FgradST(FROWS{4}) + celastST(CROWS{4},5).*FgradST(FROWS{2}) + celastST(CROWS{5},5).*FgradST(FROWS{9})) + FgradST(FROWS{6}).*(celastST(CROWS{3},4).*FgradST(FROWS{4}) + celastST(CROWS{4},4).*FgradST(FROWS{2}) + celastST(CROWS{5},4).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{4},6) = FgradST(FROWS{6}).*(celastST(CROWS{3},2).*FgradST(FROWS{4}) + celastST(CROWS{4},2).*FgradST(FROWS{2}) + celastST(CROWS{5},2).*FgradST(FROWS{9})) + FgradST(FROWS{5}).*(celastST(CROWS{3},4).*FgradST(FROWS{4}) + celastST(CROWS{4},4).*FgradST(FROWS{2}) + celastST(CROWS{5},4).*FgradST(FROWS{9})) + FgradST(FROWS{1}).*(celastST(CROWS{3},6).*FgradST(FROWS{4}) + celastST(CROWS{4},6).*FgradST(FROWS{2}) + celastST(CROWS{5},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{4},7) = FgradST(FROWS{7}).*(celastST(CROWS{3},2).*FgradST(FROWS{4}) + celastST(CROWS{4},2).*FgradST(FROWS{2}) + celastST(CROWS{5},2).*FgradST(FROWS{9})) + FgradST(FROWS{3}).*(celastST(CROWS{3},4).*FgradST(FROWS{4}) + celastST(CROWS{4},4).*FgradST(FROWS{2}) + celastST(CROWS{5},4).*FgradST(FROWS{9})) + FgradST(FROWS{8}).*(celastST(CROWS{3},6).*FgradST(FROWS{4}) + celastST(CROWS{4},6).*FgradST(FROWS{2}) + celastST(CROWS{5},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{4},8) = FgradST(FROWS{8}).*(celastST(CROWS{3},1).*FgradST(FROWS{4}) + celastST(CROWS{4},1).*FgradST(FROWS{2}) + celastST(CROWS{5},1).*FgradST(FROWS{9})) + FgradST(FROWS{3}).*(celastST(CROWS{3},5).*FgradST(FROWS{4}) + celastST(CROWS{4},5).*FgradST(FROWS{2}) + celastST(CROWS{5},5).*FgradST(FROWS{9})) + FgradST(FROWS{7}).*(celastST(CROWS{3},6).*FgradST(FROWS{4}) + celastST(CROWS{4},6).*FgradST(FROWS{2}) + celastST(CROWS{5},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{4},9) = FgradST(FROWS{9}).*(celastST(CROWS{3},1).*FgradST(FROWS{4}) + celastST(CROWS{4},1).*FgradST(FROWS{2}) + celastST(CROWS{5},1).*FgradST(FROWS{9})) + FgradST(FROWS{4}).*(celastST(CROWS{3},5).*FgradST(FROWS{4}) + celastST(CROWS{4},5).*FgradST(FROWS{2}) + celastST(CROWS{5},5).*FgradST(FROWS{9})) + FgradST(FROWS{2}).*(celastST(CROWS{3},6).*FgradST(FROWS{4}) + celastST(CROWS{4},6).*FgradST(FROWS{2}) + celastST(CROWS{5},6).*FgradST(FROWS{9}));
    celasLARGEmat(FROWS{5},5) = FgradST(FROWS{5}).*(celastST(CROWS{3},3).*FgradST(FROWS{5}) + celastST(CROWS{4},3).*FgradST(FROWS{6}) + celastST(CROWS{5},3).*FgradST(FROWS{1})) + FgradST(FROWS{1}).*(celastST(CROWS{3},5).*FgradST(FROWS{5}) + celastST(CROWS{4},5).*FgradST(FROWS{6}) + celastST(CROWS{5},5).*FgradST(FROWS{1})) + FgradST(FROWS{6}).*(celastST(CROWS{3},4).*FgradST(FROWS{5}) + celastST(CROWS{4},4).*FgradST(FROWS{6}) + celastST(CROWS{5},4).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{5},6) = FgradST(FROWS{6}).*(celastST(CROWS{3},2).*FgradST(FROWS{5}) + celastST(CROWS{4},2).*FgradST(FROWS{6}) + celastST(CROWS{5},2).*FgradST(FROWS{1})) + FgradST(FROWS{5}).*(celastST(CROWS{3},4).*FgradST(FROWS{5}) + celastST(CROWS{4},4).*FgradST(FROWS{6}) + celastST(CROWS{5},4).*FgradST(FROWS{1})) + FgradST(FROWS{1}).*(celastST(CROWS{3},6).*FgradST(FROWS{5}) + celastST(CROWS{4},6).*FgradST(FROWS{6}) + celastST(CROWS{5},6).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{5},7) = FgradST(FROWS{7}).*(celastST(CROWS{3},2).*FgradST(FROWS{5}) + celastST(CROWS{4},2).*FgradST(FROWS{6}) + celastST(CROWS{5},2).*FgradST(FROWS{1})) + FgradST(FROWS{3}).*(celastST(CROWS{3},4).*FgradST(FROWS{5}) + celastST(CROWS{4},4).*FgradST(FROWS{6}) + celastST(CROWS{5},4).*FgradST(FROWS{1})) + FgradST(FROWS{8}).*(celastST(CROWS{3},6).*FgradST(FROWS{5}) + celastST(CROWS{4},6).*FgradST(FROWS{6}) + celastST(CROWS{5},6).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{5},8) = FgradST(FROWS{8}).*(celastST(CROWS{3},1).*FgradST(FROWS{5}) + celastST(CROWS{4},1).*FgradST(FROWS{6}) + celastST(CROWS{5},1).*FgradST(FROWS{1})) + FgradST(FROWS{3}).*(celastST(CROWS{3},5).*FgradST(FROWS{5}) + celastST(CROWS{4},5).*FgradST(FROWS{6}) + celastST(CROWS{5},5).*FgradST(FROWS{1})) + FgradST(FROWS{7}).*(celastST(CROWS{3},6).*FgradST(FROWS{5}) + celastST(CROWS{4},6).*FgradST(FROWS{6}) + celastST(CROWS{5},6).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{5},9) = FgradST(FROWS{9}).*(celastST(CROWS{3},1).*FgradST(FROWS{5}) + celastST(CROWS{4},1).*FgradST(FROWS{6}) + celastST(CROWS{5},1).*FgradST(FROWS{1})) + FgradST(FROWS{4}).*(celastST(CROWS{3},5).*FgradST(FROWS{5}) + celastST(CROWS{4},5).*FgradST(FROWS{6}) + celastST(CROWS{5},5).*FgradST(FROWS{1})) + FgradST(FROWS{2}).*(celastST(CROWS{3},6).*FgradST(FROWS{5}) + celastST(CROWS{4},6).*FgradST(FROWS{6}) + celastST(CROWS{5},6).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{6},6) = FgradST(FROWS{6}).*(celastST(CROWS{2},2).*FgradST(FROWS{6}) + celastST(CROWS{4},2).*FgradST(FROWS{5}) + celastST(CROWS{6},2).*FgradST(FROWS{1})) + FgradST(FROWS{5}).*(celastST(CROWS{2},4).*FgradST(FROWS{6}) + celastST(CROWS{4},4).*FgradST(FROWS{5}) + celastST(CROWS{6},4).*FgradST(FROWS{1})) + FgradST(FROWS{1}).*(celastST(CROWS{2},6).*FgradST(FROWS{6}) + celastST(CROWS{4},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{6},7) = FgradST(FROWS{7}).*(celastST(CROWS{2},2).*FgradST(FROWS{6}) + celastST(CROWS{4},2).*FgradST(FROWS{5}) + celastST(CROWS{6},2).*FgradST(FROWS{1})) + FgradST(FROWS{3}).*(celastST(CROWS{2},4).*FgradST(FROWS{6}) + celastST(CROWS{4},4).*FgradST(FROWS{5}) + celastST(CROWS{6},4).*FgradST(FROWS{1})) + FgradST(FROWS{8}).*(celastST(CROWS{2},6).*FgradST(FROWS{6}) + celastST(CROWS{4},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{6},8) = FgradST(FROWS{8}).*(celastST(CROWS{2},1).*FgradST(FROWS{6}) + celastST(CROWS{4},1).*FgradST(FROWS{5}) + celastST(CROWS{6},1).*FgradST(FROWS{1})) + FgradST(FROWS{3}).*(celastST(CROWS{2},5).*FgradST(FROWS{6}) + celastST(CROWS{4},5).*FgradST(FROWS{5}) + celastST(CROWS{6},5).*FgradST(FROWS{1})) + FgradST(FROWS{7}).*(celastST(CROWS{2},6).*FgradST(FROWS{6}) + celastST(CROWS{4},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{6},9) = FgradST(FROWS{9}).*(celastST(CROWS{2},1).*FgradST(FROWS{6}) + celastST(CROWS{4},1).*FgradST(FROWS{5}) + celastST(CROWS{6},1).*FgradST(FROWS{1})) + FgradST(FROWS{4}).*(celastST(CROWS{2},5).*FgradST(FROWS{6}) + celastST(CROWS{4},5).*FgradST(FROWS{5}) + celastST(CROWS{6},5).*FgradST(FROWS{1})) + FgradST(FROWS{2}).*(celastST(CROWS{2},6).*FgradST(FROWS{6}) + celastST(CROWS{4},6).*FgradST(FROWS{5}) + celastST(CROWS{6},6).*FgradST(FROWS{1}));
    celasLARGEmat(FROWS{7},7) = FgradST(FROWS{7}).*(celastST(CROWS{2},2).*FgradST(FROWS{7}) + celastST(CROWS{4},2).*FgradST(FROWS{3}) + celastST(CROWS{6},2).*FgradST(FROWS{8})) + FgradST(FROWS{3}).*(celastST(CROWS{2},4).*FgradST(FROWS{7}) + celastST(CROWS{4},4).*FgradST(FROWS{3}) + celastST(CROWS{6},4).*FgradST(FROWS{8})) + FgradST(FROWS{8}).*(celastST(CROWS{2},6).*FgradST(FROWS{7}) + celastST(CROWS{4},6).*FgradST(FROWS{3}) + celastST(CROWS{6},6).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{7},8) = FgradST(FROWS{8}).*(celastST(CROWS{2},1).*FgradST(FROWS{7}) + celastST(CROWS{4},1).*FgradST(FROWS{3}) + celastST(CROWS{6},1).*FgradST(FROWS{8})) + FgradST(FROWS{3}).*(celastST(CROWS{2},5).*FgradST(FROWS{7}) + celastST(CROWS{4},5).*FgradST(FROWS{3}) + celastST(CROWS{6},5).*FgradST(FROWS{8})) + FgradST(FROWS{7}).*(celastST(CROWS{2},6).*FgradST(FROWS{7}) + celastST(CROWS{4},6).*FgradST(FROWS{3}) + celastST(CROWS{6},6).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{7},9) = FgradST(FROWS{9}).*(celastST(CROWS{2},1).*FgradST(FROWS{7}) + celastST(CROWS{4},1).*FgradST(FROWS{3}) + celastST(CROWS{6},1).*FgradST(FROWS{8})) + FgradST(FROWS{4}).*(celastST(CROWS{2},5).*FgradST(FROWS{7}) + celastST(CROWS{4},5).*FgradST(FROWS{3}) + celastST(CROWS{6},5).*FgradST(FROWS{8})) + FgradST(FROWS{2}).*(celastST(CROWS{2},6).*FgradST(FROWS{7}) + celastST(CROWS{4},6).*FgradST(FROWS{3}) + celastST(CROWS{6},6).*FgradST(FROWS{8}));
    celasLARGEmat(FROWS{8},8) = FgradST(FROWS{8}).*(celastST(CROWS{1},1).*FgradST(FROWS{8}) + celastST(CROWS{5},1).*FgradST(FROWS{3}) + celastST(CROWS{6},1).*FgradST(FROWS{7})) + FgradST(FROWS{3}).*(celastST(CROWS{1},5).*FgradST(FROWS{8}) + celastST(CROWS{5},5).*FgradST(FROWS{3}) + celastST(CROWS{6},5).*FgradST(FROWS{7})) + FgradST(FROWS{7}).*(celastST(CROWS{1},6).*FgradST(FROWS{8}) + celastST(CROWS{5},6).*FgradST(FROWS{3}) + celastST(CROWS{6},6).*FgradST(FROWS{7}));
    celasLARGEmat(FROWS{8},9) = FgradST(FROWS{9}).*(celastST(CROWS{1},1).*FgradST(FROWS{8}) + celastST(CROWS{5},1).*FgradST(FROWS{3}) + celastST(CROWS{6},1).*FgradST(FROWS{7})) + FgradST(FROWS{4}).*(celastST(CROWS{1},5).*FgradST(FROWS{8}) + celastST(CROWS{5},5).*FgradST(FROWS{3}) + celastST(CROWS{6},5).*FgradST(FROWS{7})) + FgradST(FROWS{2}).*(celastST(CROWS{1},6).*FgradST(FROWS{8}) + celastST(CROWS{5},6).*FgradST(FROWS{3}) + celastST(CROWS{6},6).*FgradST(FROWS{7}));
    celasLARGEmat(FROWS{9},9) = FgradST(FROWS{9}).*(celastST(CROWS{1},1).*FgradST(FROWS{9}) + celastST(CROWS{5},1).*FgradST(FROWS{4}) + celastST(CROWS{6},1).*FgradST(FROWS{2})) + FgradST(FROWS{4}).*(celastST(CROWS{1},5).*FgradST(FROWS{9}) + celastST(CROWS{5},5).*FgradST(FROWS{4}) + celastST(CROWS{6},5).*FgradST(FROWS{2})) + FgradST(FROWS{2}).*(celastST(CROWS{1},6).*FgradST(FROWS{9}) + celastST(CROWS{5},6).*FgradST(FROWS{4}) + celastST(CROWS{6},6).*FgradST(FROWS{2}));
    celasLARGEmat(FROWS{2},1) = celasLARGEmat(FROWS{1},2) ;
    celasLARGEmat(FROWS{3},1) = celasLARGEmat(FROWS{1},3) ;
    celasLARGEmat(FROWS{3},2) = celasLARGEmat(FROWS{2},3) ;
    celasLARGEmat(FROWS{4},1) = celasLARGEmat(FROWS{1},4) ;
    celasLARGEmat(FROWS{4},2) = celasLARGEmat(FROWS{2},4) ;
    celasLARGEmat(FROWS{4},3) = celasLARGEmat(FROWS{3},4) ;
    celasLARGEmat(FROWS{5},1) = celasLARGEmat(FROWS{1},5) ;
    celasLARGEmat(FROWS{5},2) = celasLARGEmat(FROWS{2},5) ;
    celasLARGEmat(FROWS{5},3) = celasLARGEmat(FROWS{3},5) ;
    celasLARGEmat(FROWS{5},4) = celasLARGEmat(FROWS{4},5) ;
    celasLARGEmat(FROWS{6},1) = celasLARGEmat(FROWS{1},6) ;
    celasLARGEmat(FROWS{6},2) = celasLARGEmat(FROWS{2},6) ;
    celasLARGEmat(FROWS{6},3) = celasLARGEmat(FROWS{3},6) ;
    celasLARGEmat(FROWS{6},4) = celasLARGEmat(FROWS{4},6) ;
    celasLARGEmat(FROWS{6},5) = celasLARGEmat(FROWS{5},6) ;
    celasLARGEmat(FROWS{7},1) = celasLARGEmat(FROWS{1},7) ;
    celasLARGEmat(FROWS{7},2) = celasLARGEmat(FROWS{2},7) ;
    celasLARGEmat(FROWS{7},3) = celasLARGEmat(FROWS{3},7) ;
    celasLARGEmat(FROWS{7},4) = celasLARGEmat(FROWS{4},7) ;
    celasLARGEmat(FROWS{7},5) = celasLARGEmat(FROWS{5},7) ;
    celasLARGEmat(FROWS{7},6) = celasLARGEmat(FROWS{6},7) ;
    celasLARGEmat(FROWS{8},1) = celasLARGEmat(FROWS{1},8) ;
    celasLARGEmat(FROWS{8},2) = celasLARGEmat(FROWS{2},8) ;
    celasLARGEmat(FROWS{8},3) = celasLARGEmat(FROWS{3},8) ;
    celasLARGEmat(FROWS{8},4) = celasLARGEmat(FROWS{4},8) ;
    celasLARGEmat(FROWS{8},5) = celasLARGEmat(FROWS{5},8) ;
    celasLARGEmat(FROWS{8},6) = celasLARGEmat(FROWS{6},8) ;
    celasLARGEmat(FROWS{8},7) = celasLARGEmat(FROWS{7},8) ;
    celasLARGEmat(FROWS{9},1) = celasLARGEmat(FROWS{1},9) ;
    celasLARGEmat(FROWS{9},2) = celasLARGEmat(FROWS{2},9) ;
    celasLARGEmat(FROWS{9},3) = celasLARGEmat(FROWS{3},9) ;
    celasLARGEmat(FROWS{9},4) = celasLARGEmat(FROWS{4},9) ;
    celasLARGEmat(FROWS{9},5) = celasLARGEmat(FROWS{5},9) ;
    celasLARGEmat(FROWS{9},6) = celasLARGEmat(FROWS{6},9) ;
    celasLARGEmat(FROWS{9},7) = celasLARGEmat(FROWS{7},9) ;
    celasLARGEmat(FROWS{9},8) = celasLARGEmat(FROWS{8},9) ;
    
    
    
end