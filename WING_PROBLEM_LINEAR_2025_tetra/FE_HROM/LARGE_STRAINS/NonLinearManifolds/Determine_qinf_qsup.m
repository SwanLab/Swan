function [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
    = Determine_qinf_qsup(CASES,NAMEsnap_base,DATAoffline)
%--------------------------------------------------------------------------
% function [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,...
%           MATPRO,DISP_CONDITIONS,OTHER_output] = Determine_qinf_qsup(CASES,NAMEsnap_base,DATAoffline)
%
% PURPOSE:
%   Constructs a reduced basis for a nonlinear reduced-order model (ROM),
%   by separating a linear mode (PhiLIN) from nonlinear components (PhiNON).
%   The function projects the snapshot data onto this basis to obtain
%   generalized coordinates (qINF, qSUP), which are later used for 
%   nonlinear manifold modeling, neural networks, or hyperreduction.
%
%   Unlike Determine_qinf_qsup_1SVD, here the first mode is manually selected
%   as the linear one, and the nonlinear basis is computed from the 
%   orthogonal complement using a separate SVD on the deflated snapshots.
%
% INPUTS:
%   - CASES          : Array of case identifiers (each corresponding to a folder with snapshots)
%   - NAMEsnap_base  : Base directory name prefix for accessing snapshot data
%   - DATAoffline    : Structure with offline processing parameters (e.g., error thresholds)
%
% OUTPUTS:
%   - PhiLIN         : Manually selected linear displacement mode (normalized)
%   - PhiNON         : Orthonormal nonlinear basis from SVD of orthogonal complement
%   - qINF           : Projection of snapshots onto PhiLIN (scalar coordinate)
%   - qSUP           : Projection of snapshots onto PhiNON (nonlinear coordinates)
%   - UU, SS, VV     : SVD of qSUP used for reduced modeling or nonlinear interpolation
%   - MESH           : Finite element mesh structure (from the first snapshot case)
%   - DATA           : Full input data structure from the simulation setup
%   - DOFl, DOFr     : Degrees of freedom (left/right) subject to boundary conditions
%   - OPERFE         : Finite element operators used in offline/online assembly
%   - MATPRO         : Material properties structure
%   - DISP_CONDITIONS: Displacement boundary conditions and DOF indexing
%   - OTHER_output   : Auxiliary data structure from snapshot generation
%
% METHOD:
%   1. Load the snapshot matrices (U*S*V') for each case and reconstruct the displacement snapshots.
%   2. Concatenate snapshots into a single matrix per case.
%   3. Select the second snapshot vector as the linear mode (PhiLIN), normalize it, and store it.
%   4. For each case, subtract the projection onto PhiLIN to obtain the orthogonal component.
%   5. Apply SVD (via SRSVD) to the orthogonal snapshots to extract the nonlinear basis PhiNON.
%   6. Project original snapshots onto PhiLIN and PhiNON to obtain qINF and qSUP.
%   7. Perform a final SVD (via SVDT) on qSUP to enable further reduction or model fitting.
%
% NOTES:
%   - PhiLIN is selected as the second snapshot of the first case (assuming the first is zero or trivial).
%   - The SRSVD procedure allows error-controlled compression of nonlinear modes.
%   - qINF is deduplicated using the 'unique' function before processing qSUP.
%   - A diagnostic plot of qINF is shown to validate the spread of linear coordinates.
%   - The output is typically used to train surrogate models (NNs, RBFs, etc.) or define interpolation maps.
%  
%  JAHO, 15-July-2025
%  %   See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/02_1param_BEAM.mlx
%--------------------------------------------------------------------------


% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;
SNAPdisp_orthog =cell(1,length(CASES)) ;

%  if    DATAoffline.USE_ALL_DISP_MODES_TO_COMPUTE_ECMpoints == 1
%      % We are interested in a global basis matrix, also including the
%      % constrained DOFs
%      SNAPdisp_ALL =cell(1,length(CASES)) ;
%  end




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
    
    if iproj == 1
        % Nonlinear problems
        % First nonzero colum in SNAPdisp{1} is considered the linear mode
        isnap  =2 ;
        PhiLIN = SNAPdisp{iproj}(:,isnap) ;
        n_PhiLIN = norm(PhiLIN,'fro') ;
        if n_PhiLIN >1e-6
            PhiLIN = PhiLIN/n_PhiLIN ;
        else
            error('Zero norm...Choose another snapshot')
        end
        
    end
    SNAPdisp_orthog{iproj} = SNAPdisp{iproj} - PhiLIN*(PhiLIN'*SNAPdisp{iproj}) ;
    
end

% Orthogonal complement and projection onto the nonlinear basis



TOL_BLOCK = DATAoffline.errorDISP*ones(length(SNAPdisp),1)' ;

DATAsvd=[];
[PhiNON,S,V] = SRSVD(SNAPdisp_orthog,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of displacement modes, nonlinear part =',num2str(size(PhiNON,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])

% TRAIN NEURAL NETWORK, COEFFICIENTS PROJECTION

% Now that we know the modes BasisU  = [PhiLIN,PhiNON], we can collect the
% data for training or nonlinear fitting function (either NNs, Radial Basis functios...)
qINF = cell(size(SNAPdisp_orthog)) ;
qSUP = cell(size(SNAPdisp_orthog)) ;

for iproj = 1:length(SNAPdisp_orthog)
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