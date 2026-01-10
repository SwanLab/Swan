function  [PsiSEf,PhiDEF,GammaBUB] =   BubbleModes_BasicSEfromDEF(PsiRBf,Mintf,DUMMY,DATALOC,...
    DATAoffline,MESH,BasisUdeform,Mdom,Kstiff,MATPRO,OPERFE,MintfINV_chol,PhiRB,MdomCHOL)
 %==========================================================================
% function  [PsiSEf,PhiDEF,GammaBUB] = BubbleModes_BasicSEfromDEF(...)
%
% Constructs a pair of basis matrices from deformational modes:
%
%   - PsiSEf: self-equilibrated stress modes (trace-free on internal DOFs)
%   - GammaBUB: complementary bubble modes, orthogonal to PsiSEf at interface
%
% These modes are fundamental to the EIFEM framework, enabling enriched
% model order reduction (MOR) with both static admissibility and deformation
% completeness.
%
%--------------------------------------------------------------------------
% THEORETICAL BACKGROUND:
% Based on the decomposition of the deformational space:
%       span([PhiDEF, GammaBUB]) = span([PhiDEF, PhiDEFcomp])
% where:
%   - PhiDEF      : basic elastic deformational modes
%   - PhiDEFcomp  : complement modes from partitioned SVD
%   - PsiSEf      : self-equilibrated interface stress modes, computed via
%                   weighted SVD of internal forces induced by PhiDEF
%   - GammaBUB    : modes orthogonal to PsiSEf in a zero-work or kinematic sense
%
% This function:
%   1. Constructs PsiSEf = WSVD(K * PhiDEF) restricted to interface DOFs
%   2. Verifies self-equilibration via norm ratios and optionally enforces it
%   3. Computes angles between stress and strain modes to ensure stability
%   4. Constructs bubble modes via `BubbleModes_EIFEM` if PhiDEFcomp is non-empty
%
%--------------------------------------------------------------------------
% REFERENCES:
% - J.A. Hernández Ortega, "Computation of bubble modes", Appendix (EIFEM_largeROTfinal.pdf)
% - BubbleModes_BasicSE.m, DeformModesNONL_INVAR.m
%
%--------------------------------------------------------------------------
% INPUTS:
% - PsiRBf      : rigid body equilibrium modes (for purging)
% - Mintf       : interface mass matrix (for WSVD projections)
% - DUMMY       : unused argument (legacy compatibility)
% - DATALOC     : local configuration structure
% - DATAoffline : offline data structure, incl. options (e.g., WSVD toggle)
% - MESH        : FE mesh structure with boundary DOF info
% - BasisUdeform: struct with .PhiDEFbs (basic) and .PhiDEFcomp (complement)
% - Mdom, Kstiff, MATPRO, OPERFE: domain matrices and operators
% - MintfINV_chol, PhiRB, MdomCHOL: optional Cholesky decompositions
%
%--------------------------------------------------------------------------
% OUTPUTS:
% - PsiSEf      : self-equilibrated static modes
% - PhiDEF      : basic elastic deformational modes (copied from input)
% - GammaBUB    : bubble modes orthogonal to PsiSEf, if computed
%
%--------------------------------------------------------------------------
% AUTHOR:
% Joaquín Alberto Hernández Ortega (UPC, Balmes 185, Barcelona)
%
% DATE:
% Original version: April 8th, 2024
% Comments and formatting updated: May 2025
%==========================================================================

% -------------------------------------------
if nargin == 0
    load('tmp.mat')
end




disp(['-----------------------------------------------------------'])
disp(['Determination of self-equilibrated modes'])
disp(['-----------------------------------------------------------'])

PsISEall = Kstiff*BasisUdeform.PhiDEFbs ;
b = MESH.faceDOFSall ;  % Boundary DOFs
s = setdiff(1:size(PsISEall,1),b) ; % remaining DOFs

PsISEb = PsISEall(b,:) ;
PsISEs = PsISEall(s,:) ;

ratioNORMpsi = norm(PsISEs, 'fro') / norm(PsISEb, 'fro');
fprintf('Norm PsISEs / Norm PsISEb = %.3e\n', ratioNORMpsi);
make_SELF_EQUILIBRATED = 0 ;
TOL_consider_selfequilibrated = 1e-6;
DATALOC = DefaultField(DATALOC,'ConversionINDEX_AUXDOM',[]) ; 
if ratioNORMpsi > TOL_consider_selfequilibrated
    fprintf('ERROR: norm(PsISEs) / norm(PsISEb) = %.3e exceeds tolerance of %.1e\n', ...
        ratioNORMpsi, TOL_consider_selfequilibrated);
    if  ~isempty(DATALOC.ConversionINDEX_AUXDOM)
        disp('This might be due to the fact that the deformational modes are obtained from different domains in this case')
        make_SELF_EQUILIBRATED = 1;
    else
        error('Reactive forces are not self-equilibrated according to the prescribed tolerance');
    end
end

%DATALOCww.Mchol = Mintfinv_chol ;


DATAoffline = DefaultField(DATAoffline,'USE_CHOLESKY_DECOMPOSITION',1) ; % JAHO 22-Apr-2024


if DATAoffline.USE_CHOLESKY_DECOMPOSITION == 1
    
    if isempty(MintfINV_chol)
        DATALOCww = []  ;
        MintfINV = inv(Mintf) ;
        [ PsiSEf,Sbs,~,MintfINV_chol] = WSVDT( PsISEb,MintfINV,DATALOCww) ;
    else
        DATALOCww.Mchol = MintfINV_chol ;
        %MintfINV = inv(Mintf) ;
        [ PsiSEf,Sbs,~,~] = WSVDT( PsISEb,[],DATALOCww) ;
    end
    
else
    DATALOCww.Mchol =[] ;
    %MintfINV = inv(Mintf) ;
    [ PsiSEf,Sbs,~,~] = WSVDT( PsISEb,[],DATALOCww) ;
    
    
end

% sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
if make_SELF_EQUILIBRATED == 1
    fprintf('Forcing the condition of self-equilibrated...\n');
    PsiSEf = SprojDEF_operator(PhiRB(b,:),speye(size(Mintf)),PsiSEf) ;
    
end





DATALOC.LEGEND_MODES_SE = 'Basic' ;
PlotModesSE_SubdomainLevel(DATALOC,PsiRBf,PsiSEf,MESH) ;

disp(['RECALL THAT THE APPROACH USED IS THAT IN WHICH THE  ONLY Self-equilibrated modes IN THE FORMULATION ARE BasiC Modes (ELASTIC)'])

PhiDEF = BasisUdeform.PhiDEFbs ;
% ------------------------------------------------------------------

CHECK_MATRIX_H_BASIC_MODES  = 1 ;
disp('------------------------------------------------------------------------------------------------')
if  CHECK_MATRIX_H_BASIC_MODES == 1
    % Principal angles formed by PhiDEFb and PsiSEf
    PhiDEFb = SVDT(PhiDEF(b,:)) ;
    PsiSEf_orth = SVDT(PsiSEf) ;
    [UUU,PrincipalAngles_cosines,VVV] = SVDT(PsiSEf_orth'*PhiDEFb) ;
    disp(['Cosine principal angles  formed by basic deformational and self-equilibrated modes'])
    PrincipalAngles_cosines
    TOL_limit_cosines = 1e-4 ;
    if PrincipalAngles_cosines(end) < TOL_limit_cosines
        warning(['Smallest cosine principal angle is below '],num2str(TOL_limit_cosines))
        disp(['This might be conducive to ill-posed problems'])
        pause
        
    end
    
    
    disp('------------------------------------------------------------------------------------------------')
    
    
end


CHECK_PRINCIPAL_ANGLES_ORTHOGONAL_COMPL_MODES = 0;
if CHECK_PRINCIPAL_ANGLES_ORTHOGONAL_COMPL_MODES == 1
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES    BY PsiSEf = PsISEfBS    AND   BasisUdeform.PhiDEFcomp(b,:)  (= BasisUdeform.PhiDEFbs(b,:))')
    disp('-------------------------------------------------')
    [PhiDEFcompB,SSS,VVV] =WSVDT(BasisUdeform.PhiDEFcomp(b,:),Mintf,[]) ;
    [uu,ss,vv] = SVDT(PsiSEf'*PhiDEFcompB) ;
    beta = real(acosd(ss))
end

%
%AsnapREACse = cell(size(AsnapREAC)) ;
b = MESH.faceDOFSall ;  % Boundary DOFs
% Self-equilibrated part (in principle, not necessary if one trains without external forces, just prescribed displacements)

if isempty(BasisUdeform.PhiDEFcomp)
    GammaBUB = [] ;
    
else
    
   GammaBUB = BubbleModes_EIFEM(PsiSEf,PhiDEF,BasisUdeform,b,PhiRB,Mdom,DATALOC,MESH,MdomCHOL,DATAoffline,MintfINV_chol) ;   
    
    
    
end
