function [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
    = Determine_qinf_qsup_TRANSLATION(CASES,NAMEsnap_base,DATAoffline,NAME_BASE)
%==========================================================================
% Determine_qinf_qsup_TRANSLATION is a Modification of Determine_qinf_qsup_1SVDd
% It proposes an anstaz for the nodal vector than involves a moving "bumb"
% SEe /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/18_movingLOAD.mlx
% https://chatgpt.com/share/691626f4-119c-8013-bae6-eaed95c65c42
% JAHO, 13-Nov-2025, Thursday, Balmes 185, Barcelona
% --------------------------------------------------


% -----------------------------------------------------
% File: Determine_qinf_qsup_1SVDd.m
%
% Purpose
% -------
%   Constructs the reduced-order basis and generalized coordinates (qINF, qSUP)
%   from multiple sets of displacement snapshots, properly handling
%   non-homogeneous Dirichlet boundary conditions.
%
%   This function is a refined version of `Determine_qinf_qsup_1SVD`, adapted
%   to problems where some prescribed displacements are **non-zero**. It ensures
%   the consistent extraction of linear (Φₗᵢₙ) and nonlinear (Φₙₒₙ) modes,
%   suitable for manifold-based Reduced Order Models (ROMs).
%
% Context
% -------
%   Used in nonlinear manifold-ROM workflows for both small- and
%   large-strain problems, where the reduced kinematic representation is
%   defined by:
%
%       u(x, t) ≈ Φₗᵢₙ·qₗᵢₙ(t) + Φₙₒₙ·qₙₒₙ(t)
%
%   The extracted qINF (scalar) and qSUP (vector) generalized coordinates
%   characterize the parametric evolution of the nonlinear manifold.
%
% Function Signature
% ------------------
%   [PhiLIN, PhiNON, qINF, qSUP, UU, SS, VV, ...
%    MESH, DATA, DOFl, DOFr, OPERFE, MATPRO, DISP_CONDITIONS, OTHER_output] = ...
%       Determine_qinf_qsup_1SVDd(CASES, NAMEsnap_base, DATAoffline)
%
% Inputs
% ------
%   CASES         : Array of integer IDs of training folders to process.
%   NAMEsnap_base : Base path to the snapshot directories.
%   DATAoffline   : Struct with offline parameters and tolerances:
%                     • errorDISP : SVD truncation tolerance for displacements.
%
% Outputs
% -------
%   PhiLIN          : [nDOF×1] Linear mode (first left singular vector).
%   PhiNON          : [nDOF×r] Nonlinear basis (remaining left singular vectors).
%   qINF            : [1×n_snap] Reduced coordinate on ΦLIN.
%   qSUP            : [r×n_snap] Reduced coordinates on ΦNON.
%   UU, SS, VV      : SVD of qSUP (used for further reduction or fitting).
%   MESH, DATA      : FE mesh and snapshot metadata (from first case).
%   DOFl, DOFr      : Indices of free and constrained DOFs.
%   OPERFE          : Finite element operators.
%   MATPRO          : Material properties.
%   DISP_CONDITIONS : Dirichlet boundary condition structure.
%   OTHER_output    : Auxiliary info (Dirichlet data, plotting basis, etc.).
%
% Method
% -------
%   1. Load displacement snapshots for all CASES.
%   2. Concatenate snapshots and reconstruct U*S*Vᵀ blocks.
%   3. Compute an overall SVD (for both constrained and unconstrained DOFs)
%      to obtain a visualization basis (OTHER_output.Phi_To_Plot).
%   4. Restrict each snapshot to the unconstrained DOFs (DOFl).
%   5. Apply weighted SVD (SRSVD) to extract Φ and singular values.
%   6. Split Φ = [Φₗᵢₙ Φₙₒₙ].
%   7. Compute reduced coordinates:
%         qINF = Φₗᵢₙᵀ·U,   qSUP = Φₙₒₙᵀ·U
%   8. Remove duplicates in qINF (MATLAB 'unique') and truncate qSUP.
%   9. Perform SVD on qSUP for nonlinear fitting (UU, SS, VV).
%
% Features & Notes
% ----------------
%   • Handles **non-homogeneous Dirichlet BCs** by consistently using the
%     unconstrained subspace (DOFl) for SVD extraction.
%   • Produces diagnostic plot of qINF vs snapshot index.
%   • Designed to feed manifold-learning or spline-regression routines
%     (e.g., Bsplines_EP_1stmodeDMG, RBF mappings).
%   • Optional tolerance control via DATAoffline.errorDISP.
%
% Reference Path
% --------------
%   /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/
%
% Author
% ------
%   Joaquín A. Hernández Ortega (JAHO)
%   CIMNE / Universitat Politècnica de Catalunya (UPC)
%
% History (versioned)
% -------------------
%   2025-11-07 — Barcelona, Spain — Comments section updated automatically by ChatGPT; no logic modified. (JAHO)
%   2025-07-27 — Borač, Serbia — Added support for non-homogeneous Dirichlet BCs. (JAHO)
%   2025-07-15 — Original implementation for homogeneous BCs. (JAHO)
%==========================================================================

%--------------------------------------------------------------------------



%  if    DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
%      % We are interested in a global basis matrix, also including the
%      % constrained DOFs
%      SNAPdisp_ALL =cell(1,length(CASES)) ;
%  end

if nargin == 0
    load('tmp.mat')
    close all
end


% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
SNAPdisp_orthog =cell(1,length(CASES)) ;
for iproj = 1:length(CASES)
    NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    if iproj == 1
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'MATPRO','MESH','OPERFE','DISP_CONDITIONS',...
            'OTHER_output') ;
        DATA = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ;
        DOFl = DISP_CONDITIONS.DOFl ;
        DOFr = DISP_CONDITIONS.DOFr ;
        OTHER_output.DISP_CONDITIONS = DISP_CONDITIONS ;
    end
    
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    DISP_LOC = cell(1,length(NAME_SNAP_loc)) ;
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
        
        %DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        
        % Or the whole matrix ....
        DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U',SNAP_cluster.DISP.S)' ;
        DISP_LOC{iloc} =  DISP_LOC{iloc}*SNAP_cluster.DISP.V(2:end,:)' ;
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    
    
    
end

% Let us begin by extracting, for instance, the x coordinates of the top
% fiber
SNAPdisp = cell2mat(SNAPdisp) ;
iSURFtop = 3;
NodesTOP = MESH.NODES_FACES{iSURFtop}  ;
ndim = size(MESH.COOR,2) ;
DofsTOP = small2large(NodesTOP,ndim) ;
CoorTOP = MESH.COOR(NodesTOP,1) ;
disTOP = SNAPdisp(DofsTOP(2:ndim:end),:) ;

figure(109)
hold on
xlabel('x')
ylabel('displacement')
title('Displacement (y) top surface versus x')
freq = 9 ;
TIme_select = 1:freq:size(disTOP,2)

for itime = 1:length(TIme_select)
    plot(CoorTOP,disTOP(:,TIme_select(itime)),'DisplayName',['FE, istep = ',num2str(TIme_select(itime))]) ;
end

[maxDISP,indXmax] = max(abs(disTOP)) ;
qMASTER=  CoorTOP(indXmax) ; % This play the role of latent variable (master variable)
% Let us begin by fitting the maximum displacements at each time step
figure(285)
hold on
xlabel('qMASTER')
ylabel('Amplitude RBF')
plot(qMASTER,maxDISP,'DisplayName','FE')
grid on

% Not worthy
%  [c0,c1,pshape] =  Find_Parameters_amplitude(qMASTER,maxDISP) ; 
%  
%  AMp_approx   = c0 + c1*(qMASTER.*(1-qMASTER)).^pshape ; 
% plot(qMASTER,AMp_approx,'DisplayName',['c0+c1(q(1-q))^p'])
% legend show




% But is it worthy to go ino this ?  Let us first see if it can handle one
% single function with a Gaussian Kernel
istep = 50 ;
x = CoorTOP ;
y = disTOP(:,istep) ;
x0 = qMASTER(istep) ;
A = -maxDISP(istep) ;
kernel = 'gaussian' ;


[sigma_opt, J_opt] = fit_rbf_width(x, y, x0, A, kernel ) ;

figure(109)
hold on
%   xlabel('x')
%  ylabel('dTOP')


for istep = 1:length(TIme_select)
    istepLOC = TIme_select(istep);
    % A = -maxDISP(istepLOC) ;
    x0 = qMASTER(istepLOC) ;
    %   y = disTOP(:,istepLOC) ;
    dAPPROX =  A*exp(-0.5*(x-x0).^2/sigma_opt^2) ;
    plot(x,dAPPROX,'DisplayName',['Gauss. approx = ',num2str(istepLOC)],'LineStyle','--')
    %    plot(x,y,'DisplayName','Exact ')
    
    %      % Entire vector of displacement
    %      x = MESH.COOR(:,1) ;
    %      dGAUSSIAN_y = A*exp(-0.5*(x-x0).^2/sigma_opt^2) ;
    %      y = MESH.COOR(:,2) ;
    %      yMIN = min(y) ;
    %      yMAX = max(y) ;
    %
    %      dLINEAR_y = (y-yMIN)/(yMAX-yMIN) ;
    %
    %      SNAPdispGAUSS(2:ndim:end,:) = dLINEAR_y(:).*dGAUSSIAN_y(:) ;
    
    
    
    
end
legend show

SNAPdispGAUSS = zeros(size(SNAPdisp)) ;

for istep = 1:size(SNAPdispGAUSS,2)
    % A = -maxDISP(istepLOC) ;
    x0 = qMASTER(istep) ;
    %   y = disTOP(:,istepLOC) ;
    dAPPROX =  A*exp(-0.5*(x-x0).^2/sigma_opt^2) ;
    %   plot(x,dAPPROX,'DisplayName',['Gauss. approx = ',num2str(istepLOC)],'LineStyle','--')
    %    plot(x,y,'DisplayName','Exact ')
    
    % Entire vector of displacement
    x = MESH.COOR(:,1) ;
    dGAUSSIAN_y = A*exp(-0.5*(x-x0).^2/sigma_opt^2) ;
    y = MESH.COOR(:,2) ;
    yMIN = min(y) ;
    yMAX = max(y) ;
    
    dLINEAR_y = (y-yMIN)/(yMAX-yMIN) ;
    
    SNAPdispGAUSS(2:ndim:end,istep) = dLINEAR_y(:).*dGAUSSIAN_y(:) ;
    
end
% PLOTTING RESULT IN GID
% ---------------------------------
NAME_MODES_FOLDER = [cd,filesep,'GIDPOST',filesep];
if  ~exist(NAME_MODES_FOLDER)
    mkdir(NAME_MODES_FOLDER)
end
NAME__DISP = [NAME_MODES_FOLDER,NAME_BASE,'DISPgaussian'] ;

% NameFileMesh = [NAME_MODES_DISP,'.msh'] ;
% NameFile_res = [NAME_MODES_DISP,'.res'] ;

t = 1:size(SNAPdispGAUSS,2)  ;
NAME_INPUT_DATA = [] ;
GidPostProcessDynamic_loc(MESH.COOR,MESH.CN,MESH.TypeElement,SNAPdispGAUSS,NAME__DISP,DATA.MESH.posgp,NAME__DISP,t);


% NOW WE APPLY THE SVD TO THE DIFFERENCE BETWEEN THE FE DISPLACEMENTS
% AND THE GAUSSIAN ANSTAZ
SNAPdisp  =  SNAPdisp(DOFl,:) ;
SNAPdispGAUSS  =  SNAPdispGAUSS(DOFl,:) ;

SNAPdisp_FLUCT = SNAPdisp-SNAPdispGAUSS ;

errorGAUSS  = norm(SNAPdisp_FLUCT,'fro') ; 
errorGAUSS = errorGAUSS/norm(SNAPdisp,'fro') ; 
disp(['Error made by Gaussian predictions (over 1) = ',num2str(errorGAUSS)])


TOL_BLOCK = DATAoffline.errorDISP  ;
% 
% % This is just for plotting purposes.
% % We apply the SVD to a matrix containing both unconstrained and
% % constrained DOFs
% DATAsvd=[];
% [OTHER_output.Phi_To_Plot,S,V] = SRSVD(SNAPdisp,TOL_BLOCK,DATAsvd) ;


% SVD of the unconstrained block of SNAPdisp

    
 




[PhiDEV,S,V] = SRSVD(SNAPdisp_FLUCT,TOL_BLOCK) ;
disp('***********************************************************')
disp(['Number of displacement modes for deviations with respect Gaussian =',num2str(size(PhiDEV,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])

% d = Decoder(q) = d_GAUSSIAN(q) + d_DEVIATIONS(q)
% where d_DEVIATIONS = PhiDEV*S*V(q)^T
% Derivative
%  d_Decoder/dq = der(d_GAUSSIAN(q))/dq +  PhiDEV*S*der(V^T(q))


error('To be continued')
 %  [DATA_evaluateTAU_and_DER, nREDcoor] = BsplinesLeastSquares_SHIFT(DATA_interp, qINF, VrightT, UU, SS)