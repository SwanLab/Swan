function [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
    = Determine_qinf_qsup_MACRO(CASES,NAMEsnap_base,DATAoffline)
% Determine_qinf_qsup_MACRO is a modification of Determine_qinf_qsup_1SVD
% In homogenization problems, rather than using as latent variable the amplitude of the first SVD mode,
% we use macroscopic strains themselves 
% JAHO, 5-Oct-2025, Sunday, 11:23, Secrets by Farga, Barcelona. 
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This is a modification of 
% function Determine_qinf_qsup_1SVD 
%  whose description is shown below. The modification is related with the
%  fact that the former version was not able to deal with non-homogeneous
%  Dirichlet boundary conditions
%  Modification made on July 27th 2025, Sunday, Boraĉ, Serbia 


% function [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,...
%           MATPRO,DISP_CONDITIONS,OTHER_output] = Determine_qinf_qsup_1SVD(CASES,NAMEsnap_base,DATAoffline)
%
% PURPOSE:
%   This function constructs a reduced-order basis from multiple sets of 
%   displacement snapshots, computes a linear/nonlinear decomposition via SVD,
%   and generates projection coefficients used for reduced-order modeling.
%   The linear component is treated separately, and the remaining modes
%   are used to define the nonlinear reduced subspace.
%
%   This is a key step in nonlinear manifold-based ROMs, where the 
%   generalized coordinates (qINF, qSUP) are extracted for each case.
%
% INPUTS:
%   - CASES          : Array of case identifiers to loop over (each corresponding to a folder with snapshots)
%   - NAMEsnap_base  : Base path to snapshot directories (e.g., 'DATA_snap_case_')
%   - DATAoffline    : Structure containing tolerances and offline settings for SVD and error threshold
%
% OUTPUTS:
%   - PhiLIN         : First left singular vector, used as linear mode
%   - PhiNON         : Remaining left singular vectors (nonlinear subspace)
%   - qINF           : Projections of snapshots onto PhiLIN (scalar generalized coordinate)
%   - qSUP           : Projections of snapshots onto PhiNON (nonlinear coordinates)
%   - UU, SS, VV     : SVD of qSUP used for dimensionality reduction or nonlinear mappings
%   - MESH           : Mesh structure from the first case
%   - DATA           : Complete input data structure used in snapshot generation
%   - DOFl, DOFr     : Left and right Dirichlet boundary conditions
%   - OPERFE         : Operators used for finite element assembly or ROM evaluation
%   - MATPRO         : Material properties structure
%   - DISP_CONDITIONS: Structure specifying displacement constraints
%   - OTHER_output   : Any other relevant auxiliary data loaded from snapshots
%
% METHOD:
%   1. Load snapshot matrices for each case in CASES.
%   2. Concatenate the displacement data from each location.
%   3. Build a global snapshot matrix and apply SVD to extract modes.
%   4. The first mode (PhiLIN) is stored as the linear component.
%   5. Remaining modes (PhiNON) define the nonlinear manifold.
%   6. Project each snapshot onto PhiLIN and PhiNON to compute qINF and qSUP.
%   7. Perform SVD of qSUP to extract a compact representation (UU, SS, VV).
%
% NOTES:
%   - Snapshots are assumed to be stored as U*S*V' matrices, possibly preweighted.
%   - If available, the function can later be extended to handle full-rank or truncated SVD variants.
%   - The output qINF is cleaned from duplicates using MATLAB’s 'unique' function.
%   - Diagnostic plot is produced to visualize the evolution of qINF.
%   JAHO, 15th July 2025
%   See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/02_1param_BEAM.mlx
%--------------------------------------------------------------------------


% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
SNAPdisp_orthog =cell(1,length(CASES)) ;  

%  if    DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
%      % We are interested in a global basis matrix, also including the
%      % constrained DOFs
%      SNAPdisp_ALL =cell(1,length(CASES)) ;
%  end

if nargin == 0
    load('tmp.mat')
end


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
        DISP_LOC{iloc} =  DISP_LOC{iloc}*SNAP_cluster.DISP.V' ;
        
    end
    SNAPdisp{iproj} = cell2mat(DISP_LOC);
    
%     if iproj == 1
%         % Nonlinear problems
%         % First nonzero colum in SNAPdisp{1} is considered the linear mode
%         isnap  =2 ;
%         PhiLIN = SNAPdisp{iproj}(:,isnap) ;
%         n_PhiLIN = norm(PhiLIN,'fro') ;
%         if n_PhiLIN >0
%             PhiLIN = PhiLIN/n_PhiLIN ;
%         else
%             error('Zero norm...Choose another snapshot')
%         end
%         
%     end
%     SNAPdisp_orthog{iproj} = SNAPdisp{iproj} - PhiLIN*(PhiLIN'*SNAPdisp{iproj}) ;
    
end

% Orthogonal complement and projection onto the nonlinear basis



TOL_BLOCK = DATAoffline.errorDISP*ones(length(SNAPdisp),1)' ;

% This is just for plotting purposes. 
% We apply the SVD to a matrix containing both unconstrained and
% constrained DOFs
DATAsvd=[];
[OTHER_output.Phi_To_Plot,S,V] = SRSVD(SNAPdisp,TOL_BLOCK,DATAsvd) ;


% SVD of the unconstrained block of SNAPdisp

for irpo = 1:length(SNAPdisp)
     SNAPdisp{irpo} =  SNAPdisp{irpo}(DOFl,:) ; 
end



DATAsvd=[];
[Phi,S,V] = SRSVD(SNAPdisp,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of displacement modes =',num2str(size(Phi,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])



PhiLIN = Phi(:,1) ; 
PhiNON = Phi(:,2:end) ; 

% TRAIN NEURAL NETWORK, COEFFICIENTS PROJECTION

% Now that we know the modes BasisU  = [PhiLIN,PhiNON], we can collect the
% data for training or nonlinear fitting function (either NNs, Radial Basis functios...)
qINF = cell(size(SNAPdisp)) ;
qSUP = cell(size(SNAPdisp)) ;

for iproj = 1:length(SNAPdisp)
    qINF{iproj} = PhiLIN'*SNAPdisp{iproj} ;
    qSUP{iproj} = PhiNON'*SNAPdisp{iproj} ;
end

qINF = cell2mat(qINF) ;
[qINF,aaa,bbb ]= unique(qINF) ;

figure(10)
hold on
xlabel('Number of snapshot')
ylabel('qINF')
plot(qINF,'Marker','*')
grid on

qSUP = cell2mat(qSUP) ;
qSUP = qSUP(:,aaa) ;

[UU,SS,VV ] = SVDT(qSUP) ;
%SS = SS/SS(1) ;