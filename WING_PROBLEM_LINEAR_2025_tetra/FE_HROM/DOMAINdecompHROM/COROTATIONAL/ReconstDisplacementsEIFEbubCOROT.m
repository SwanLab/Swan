function [DISP,strainC,dClocINCREe_time] = ReconstDisplacementsEIFEbubCOROT(XfREF,dCeTIME,EIFE_prop,QrotINIe,DATA,lambdaLENe,QrotTIME,ndim,eIND)
%  Adaptation of ReconstDisplacementsEIFEbub.m to co-rotational approach
% See % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 6-Nov-2024, UPC, Terrassa
% Comments by ChatGPT-4
%ReconstDisplacementsEIFEbubCOROT Reconstruct fine-scale displacements using co-rotational EIFE with bubble enrichment
%
%   This function reconstructs the fine-scale displacement field of an EIFEM domain,
%   enriched with bubble functions and corrected for large rotations using a 
%   co-rotational transformation. It combines rotational motion, rigid body modes,
%   and strain-induced deformation to build the final nodal displacements over time.
%
%   INPUT:
%     XfREF         - Reference coordinates of the fine-scale mesh (parent domain)
%     dCeTIME       - Time-varying coefficients of coarse-scale DOFs (boundary + bubble)
%     EIFE_prop     - EIFE properties of the parent domain (operators, meshes, etc.)
%     QrotINIe      - Initial rotation matrix of the domain (local → global)
%     DATA          - Structure with simulation and reconstruction options
%     lambdaLENe    - Scaling factor matrix for this domain
%     QrotTIME      - Cell array with rotation matrices per time step
%     ndim          - Spatial dimension (e.g., 2 or 3)
%     eIND          - Indices in QrotTIME associated to this domain (row slicing)
%
%   OUTPUT:
%     DISP              - Matrix of fine-scale displacements at all nodes and time steps
%     strainC           - Time evolution of strain-mode coefficients (including bubbles)
%     dClocINCREe_time  - Time evolution of corrected coarse-scale DOFs (local reference frame)
%
%   PROCEDURE:
%   -----------------------------------------------------------------------------------------
%   1. For each time step, retrieve the rotation matrix QrotE and compute:
%      - The rotational displacement of the fine mesh nodes: Δu_rot = λ(QrotE - QrotINI)·Xref
%   2. Separate dCe into:
%      - Boundary DOFs (dCbE)
%      - Bubble DOFs (dCbubE)
%   3. Correct boundary DOFs by subtracting the rigid rotation and projecting into local frame
%   4. Form the complete incremental displacement vector (coarse local): ΔdC_loc = [ΔdCbE; dCbubE]
%   5. Reconstruct:
%      - Rigid body component: Φ_RB·(Pdown_RB·ΔdC_loc)
%      - Strain component:     Φ_DEF·(Pdown_DEF·ΔdC_loc)
%   6. Rotate fine-scale field back to global frame and add rotational part from Step 1
%
%   REMARKS:
%   - This function assumes the EIFE domain has been preprocessed with the necessary
%     reconstruction operators (basis and projection coefficients).
%   - The reconstructed displacement includes contributions from both geometric rotations
%     and physical deformation (strain + bubble).
%   - All vector-matrix operations follow the notation used in EIFEM_largeROT.tex
%
%   AUTHOR:
%     Joaquín A. Hernández Ortega (JAHO), UPC BarcelonaTech, 6-Nov-2024
%     Adapted from ReconstDisplacementsEIFEbub.m to support large rotation recovery.
%
%   SEE ALSO:
%     ReconstStressesEIFE_exactBUBcorot, PostProcess_EIFE_vectBUBcorot,
%     EIFE_operatorsBUB.m, DOMAINdecompHROM

% -----------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
% deformational part (amplitude modes)

% BOUNDARY DOFS AND COARSE-SCALE DOFS
% Defined in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/EIFE_operatorsBUB.m
% EIFEoper.INFO.DOFsBUB = DOFsBUB; % BUBBLE DOFS
% EIFEoper.INFO.DOFsB= DOFsB; % boundary DOFS
DOFsBUB= EIFE_prop.INFO.DOFsBUB;
DOFsB= EIFE_prop.INFO.DOFsB;
IndexesStrainModes= EIFE_prop.INFO.IndPhiDEF;
IndexesBubbleModes= EIFE_prop.INFO.IndGammaBUB;


DISP = zeros(size(EIFE_prop.RECONSTRUCTION.DEF_DISP.BASIS,1),size(dCeTIME,2)) ;
nstrainTOTAL = length(IndexesStrainModes) + length(IndexesBubbleModes) ; 
strainC = zeros(nstrainTOTAL,size(dCeTIME,2)) ;
dClocINCREe_time = zeros(length(DOFsBUB)+length(DOFsB),size(dCeTIME,2)) ;


for istep = 1:size(dCeTIME,2)
    dCe = dCeTIME(:,istep) ;
    
    % ------------------------
    % Rotational displacement
    % ------------------------
    % See
    % /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXT/EIFEM_largeROT.tex
    %  \dFqNe{I}{e}(t)  =  \lambdaLENe{e}  ({\QrotE{e}(t)  -\QrotINIe{e}}) \XfREFnE{I}{e}, \hspace{0.5cm}  I =1,2 \ldots \nnodeF
    if ~isempty( QrotTIME{istep})
        QrotE = QrotTIME{istep}(eIND,:) ;
    else
        QrotE = eye(ndim) ;
    end
    dFqE = lambdaLENe*(QrotE-QrotINIe)*XfREF';
    % ------------------------------
    % INCREMENTAL DISPLACEMENT
    % -----------------------------
    % COARSE-SCALE
    
    % 1) separate   $\dCe{e}$ into its boundary and bubble DOFs:
    %     \dCe{e} = \coldos{\dCbE{e}}{\dCbubE{e}}
    
    dCbE = dCe(DOFsB) ;
    dCbE = reshape(dCbE,ndim,[]) ;
    dCbubE = dCe(DOFsBUB) ;
    
    %  2) Compute the rotational displacement of the interface boundary DOFs;
    %  \dCqBnE{i}{e}(t) & =  \lambdaLENe{e}   (   \QrotE{e}(t)  -\QrotINIe{e}) \XcREFnE{i}{e},
    XcREFallBe = EIFE_prop.MESH.XcREFallBe ;
    dCqBe = lambdaLENe*(QrotE-QrotINIe)*XcREFallBe';
    
    %  3) Determine the incremental, coarse-scale nodal displacement (boundary DOFs) in the local reference system:
    % \Delta \dCbLOCnE{i}{e}(t) =  \QrotE{e}(t)^T  (\dCnE{i}{e}(t) - \dCqBnE{i}{e}(t) ),   \hspace{0.5cm} i = 1,2 \ldots \nnodeCb
    Delta_dCbLOCe = QrotE'*(dCbE-dCqBe) ;
    
    % 4) Construct the incremental vector of coarse-scale displacements, including bubble DOFs as
    %  \dClocINCREe{e}(t) = \coldos{\Delta \dCbLOCe{e}(t)}{\dCbubE{e}(t)}
    dClocINCREe = [Delta_dCbLOCe(:);dCbubE ] ;
    dClocINCREe_time(:,istep) = dClocINCREe ; 
    
    % 5) Fine-scale local displ.   $\dFlocINCREe{e}(t)$
    % \dFlocINCREe{e}(t) =  \PhiDEFallE{e} (\PdownsDEF^{(e)}  \dClocINCREe{e}(t))  +  \Rrb^{(e)}  (\PdownsRB^{(e)}   \dClocINCREe{e}(t))
    %  Here $\PhiDEFallE{e}$, $\PdownsDEF^{(e)}$, $\Rrb^{(e)}$ and $\PdownsRB^{(e)}$ are intrinsic variables of the parent domain.
    
    % RIGID BODY
    % -------------
    % %  See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/EIFE_operatorsBUB.m
    % %  EIFEoper.RECONSTRUCTION.RB_DISP.coeff = PdownsRB ;
    % EIFEoper.RECONSTRUCTION.RB_DISP.BASIS = MODES.PhiRB ;
    aRBe = EIFE_prop.RECONSTRUCTION.RB_DISP.coeff*dClocINCREe ;
    dFe = EIFE_prop.RECONSTRUCTION.RB_DISP.BASIS*aRBe ;
    % DEFORMATIONAL
    %--------------
    strainCe = EIFE_prop.RECONSTRUCTION.DEF_DISP.coeff*dClocINCREe ;
    strainC(:,istep) = strainCe ; 
    dFe = dFe +  EIFE_prop.RECONSTRUCTION.DEF_DISP.BASIS*strainCe ;
    % EIFEoper.RECONSTRUCTION.DEF_DISP.coeff = PdownsDEF ;
    % EIFEoper.RECONSTRUCTION.DEF_DISP.BASIS = PhiDEF ;
    
    
    % 6) Once we have $\dFlocINCREe{e}(t)$, we compute its components in the global system,
    %by multiplying   by $\QrotE{e}$, and sum it up to the total displacement.
    dFe = QrotE*reshape(dFe,ndim,[]) ;
    
    dFe = dFqE + dFe ;
    
    
    DISP(:,istep) = dFe(:) ;
    
end

%