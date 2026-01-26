function [PhiLIN,PhiNON,qINF,qSUP,UU,SS,VV,MESH,DATA,DOFl,DOFr,OPERFE,MATPRO,DISP_CONDITIONS,OTHER_output]...
    = Determine_qinf_qsup_LINd(CASES,NAMEsnap_base,DATAoffline,DATA_interp)
%--------------------------------------------------------------------------
% Determine_qinf_qsup_LINd is a modification of Determine_qinf_qsup_1SVDd,
% described below
% Given a set of "projects", this version determine basis matrices by
% selecting the first snapshot of each project
% 26-Sept-2025, Molinos Marfagones,  Cartagena
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
if nargin == 0
    load('tmp.mat')
end

% Let us construct first the matrix the snapshots of each project
SNAPdisp =cell(1,length(CASES)) ;


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
        OTHER_output.MESH = MESH; 
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
    
    
end


% Linear modes
% -------------
DATA_interp = DefaultField(DATA_interp,'IND_PROJECT_LINEAR',1) ;
SNAPdispLIN = cell2mat(SNAPdisp(DATA_interp.IND_PROJECT_LINEAR)) ;
SNAPdispLINplot = SNAPdispLIN ; 
SNAPdispLIN = SNAPdispLIN(DOFl,:) ;
% -------------------------------------------------------------
TOLloc = 1e-3;
DATAlll.RELATIVE_SVD = 1;
[PhiLIN,SSSS] = SVDT(SNAPdispLIN,TOLloc,DATAlll) ;

[PhiLINplot,SSSS] = SVDT(SNAPdispLINplot,TOLloc,DATAlll) ;


% Orthogonal complement (Frobenius norm)
% --------------------------------------
CASES_compl = setdiff(CASES,DATA_interp.IND_PROJECT_LINEAR) ; 
SNAPdispORTH  = cell(1,length(CASES_compl)) ;
SNAPdispORTHplot  = cell(1,length(CASES_compl)) ;

for iprojLOC = 1:length(CASES_compl)
    iproj = CASES_compl(iprojLOC) ; 
    SNAPdispORTH{iprojLOC} = SNAPdisp{iproj}(DOFl,:) - PhiLIN*(PhiLIN'*SNAPdisp{iproj}(DOFl,:)) ; 
    SNAPdispORTHplot{iprojLOC} = SNAPdisp{iproj}(:,:) - PhiLINplot*(PhiLINplot'*SNAPdisp{iproj}(:,:)) ; 
end


TOL_BLOCK = DATAoffline.errorDISP*ones(length(SNAPdispORTH),1)' ;

DATAsvd=[];
[PhiNON,S,V] = SRSVD(SNAPdispORTH,TOL_BLOCK,DATAsvd) ;
disp('***********************************************************')
disp(['Number of linear displacement modes =',num2str(size(PhiLIN,2))])
disp('***********************************************************')
 
 
disp(['Number of nonlinear displacement modes =',num2str(size(PhiNON,2)), ' (for ERRORdisp = ',num2str(DATAoffline.errorDISP),')'])
disp('***********************************************************')
disp(['Singular Values =',num2str(S')])

[PhiNONplot,S,V] = SRSVD(SNAPdispORTHplot,TOL_BLOCK,DATAsvd) ;


OTHER_output.Phi_To_Plot = [PhiLINplot,PhiNONplot] ; 


qINF = cell(size(SNAPdispORTH)) ;
qSUP = cell(size(SNAPdispORTH)) ;

for iprojLOC = 1:length(CASES_compl)
    iproj = CASES_compl(iprojLOC) ; 
     qINF{iproj} = PhiLIN'*SNAPdisp{iproj}(DOFl,2:end) ;
    qSUP{iproj} = PhiNON'*SNAPdisp{iproj}(DOFl,2:end) ;
end

qINF = cell2mat(qINF) ;
[qINF,aaa,bbb ]= unique(qINF) ;
% We remove points at which the gradient is very small (at the origin)


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

end 