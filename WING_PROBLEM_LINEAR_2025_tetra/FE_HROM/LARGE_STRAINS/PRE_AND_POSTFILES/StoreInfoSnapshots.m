function [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: StoreInfoSnapshots
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   Stores selected field variables (e.g., displacements, strains, stresses) into a snapshot matrix
%   structure for a given time step and manages compression of these snapshots **cluster by cluster**.
%   Optionally performs SVD compression of the snapshot data for use in **reduced-order modeling (ROM)**
%   via Empirical Cubature or Proper Orthogonal Decomposition (POD).
%
%   This function is typically called after each time step in a nonlinear simulation.
%
% USAGE:
%   [SNAP, DATA, icluster] = StoreInfoSnapshots(istep, icluster, SNAP, VAR, DATA, CONVERGED)
%
% INPUTS:
%   - istep     : Current global time step index.
%   - icluster  : Current cluster index (group of time steps).
%   - SNAP      : Structure containing the accumulating snapshot matrices for current cluster.
%   - VAR       : Structure containing current state variables to be stored (e.g., DISP, PK2STRESS, etc.).
%   - DATA      : Structure containing the time-stepping and storage configuration.
%   - CONVERGED : Flag (1 = converged, 0 = failed) or structure (future-proofing).
%
% OUTPUTS:
%   - SNAP      : Updated snapshot structure (may be reset if a new cluster begins).
%   - DATA      : Updated data structure (includes saved file names, step logs, etc.).
%   - icluster  : Updated cluster index if the time step starts a new cluster.
%
% FUNCTIONALITY:
%   - Checks whether current `istep` belongs to current `icluster`.
%   - If so and `CONVERGED == 1`, stores all fields from `VAR` into `SNAP`.
%   - When the last time step in the cluster is reached:
%       * If `DATA.STORE.COMPRESS_WITH_SVD == 1`, compresses each field in `SNAP` using `RSVDT`.
%       * Saves the compressed `SNAP_cluster` or raw `SNAP` to disk using the specified MAT filename.
%   - If `istep` starts a **new cluster**, initializes `SNAP` with zero arrays of the appropriate size.
%
% SPECIAL CASES:
%   - If SVD compression leads to an empty decomposition (e.g., constant or zero snapshot), a dummy
%     `U = 0`, `S = 0`, `V = 0` structure is saved to maintain format consistency (patch added 10-Apr-2023).
%
% NOTES:
%   - Snapshot data are saved in files named in `DATA.STORE.NAME_MATFILE_STORE{icluster}`.
%   - The stored snapshots (compressed or not) are used for building ROMs via ECM, Gappy POD, or SVD/PCA.
%
% REFERENCES:
%   - Snapshot management and compression strategy adapted from:
%     /TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE  
%   Date: April 2023 (SVD patch), preamble added May 2025 (ChatGPT4)
%   
% DEPENDENCIES:
%   - RSVDT (Randomized SVD function)
%   - DefaultField
%   - fieldnames, bsxfun (for compact SVD reconstruction)
%
% ---------------------------------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
end

[ISCURRENT_CLUSTER,istepLOC] = ismember(istep,DATA.STORE.NSTEPS_CLUSTER{icluster}) ;

if ISCURRENT_CLUSTER
    
    %      if ~isfield(SNAP,'DISP')
    %          SNAP.DISP(:,istepLOC) = d;
    %      end
    %       if ~isfield(SNAP,'PK2STRESS')
    %          SNAP.PK2STRESS(:,istepLOC) = StwoST;
    %       end
    %       if ~isfield(SNAP,'GLSTRAINS')
    %          SNAP.GLSTRAINS(:,istepLOC) = EgreenlST;
    %       end
    
    if ~isstruct(CONVERGED) && CONVERGED == 1
        fff = fieldnames(SNAP) ;
        
        for iii = 1:length(fff)
            nrows = size(VAR.(fff{iii})) ;
            SNAP.(fff{iii})(1:nrows,istepLOC) = VAR.(fff{iii}) ;
        end
    end
    
    if istepLOC == length(DATA.STORE.NSTEPS_CLUSTER{icluster})  || CONVERGED == 0
        % STORE information previous cluster
        if DATA.STORE.COMPRESS_WITH_SVD  == 1
            disp(['Compressing snapshots cluster = ',num2str(icluster),'...'])
            abcd = tic ; 
            fff= fieldnames(SNAP) ;
            for iii = 1:length(fff)
                % [U,S,V,eSVD,Rsup] = RSVDT(A,e0,mu,R,DATA)
                DATASVD.RELATIVE_SVD = 1;
                
                [U,S,V] = RSVDT(SNAP.(fff{iii}),DATA.STORE.TOLERANCE_SVD_COMPRESSION.(fff{iii}),0,1,DATASVD) ;
                disp([fff{iii},' RANK = ',num2str(length(S))])
                
                %   SV = bsxfun(@times,V',S)';
                
                SNAP_cluster.(fff{iii}) = [] ;
                if ~isempty(S) > 0 
                SNAP_cluster.(fff{iii}).U = U ;
                SNAP_cluster.(fff{iii}).S = S ;
                SNAP_cluster.(fff{iii}).V = V ;
                else
                    % MOdification 10-APRIL-2023
                     SNAP_cluster.(fff{iii}).U = zeros(size(SNAP.(fff{iii}),1),1) ;
                SNAP_cluster.(fff{iii}).S = 0 ;
                SNAP_cluster.(fff{iii}).V = zeros(size(SNAP.(fff{iii}),2),1) ; ;
                    
                end
            end
            TimeStore = toc(abcd) ;
            disp(['Done in ',num2str(TimeStore)])
            STEP_LOC = DATA.STORE.NSTEPS_CLUSTER{icluster}(1:istepLOC);
            save(DATA.STORE.NAME_MATFILE_STORE{icluster},'STEP_LOC','SNAP_cluster') ;
        else
              STEP_LOC = DATA.STORE.NSTEPS_CLUSTER{icluster}(1:istepLOC);
            save(DATA.STORE.NAME_MATFILE_STORE{icluster},'STEP_LOC','SNAP') ;
         %   disp('Option not implemented')
        end
        % SAve information
        
      
    end
    
else
    
    fff= fieldnames(SNAP) ;
    % Change storing cluster
    icluster = icluster + 1;
    % Initialization new snapshot
    ntimestepsLOC = length(DATA.STORE.NSTEPS_CLUSTER{icluster}) ;
    for iii = 1:length(fff)
        ncomp = size(SNAP.(fff{iii}),1) ;
        SNAP.(fff{iii}) = zeros(ncomp,ntimestepsLOC) ;
    end
    
    [ISCURRENT_CLUSTER,istepLOC] = ismember(istep,DATA.STORE.NSTEPS_CLUSTER{icluster}) ;
    
    if CONVERGED == 1
        fff = fieldnames(SNAP) ;
        
        for iii = 1:length(fff)
            SNAP.(fff{iii})(:,istepLOC) = VAR.(fff{iii}) ;
        end
    end
end