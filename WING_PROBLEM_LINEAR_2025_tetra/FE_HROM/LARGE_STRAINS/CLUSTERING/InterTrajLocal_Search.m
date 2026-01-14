function  [indCANDIDATES_self,indCANDIDATES_rest] =  InterTrajLocal_Search(d,nd,SNAPdisp,DATAloc,ipoints,WHICH_TRAJ_IS_EACH_SNAP,...
    IndClusters,DISTANCES_POINTS)

if nargin == 0
    load('tmp1.mat')
    DATAloc = DefaultField(DATAloc,'InterTrajectorySearchWithTwoSelfPoints',1) ;
    DATAloc = DefaultField(DATAloc,'MinNumberNeighbors_interTRAJ',4) ;
    
    %  DATAloc.TolStopSearch_interTRAJ = 1e-5;
end

% ---------------------------------------------------
% Search over points belonging to the same trajectory
% ---------------------------------------------------
itraj = WHICH_TRAJ_IS_EACH_SNAP(ipoints) ;  % Current trajectory
% Indexes snapshots pertaining to this  trajectory
INDEXsamTRAJ = find(WHICH_TRAJ_IS_EACH_SNAP == itraj) ;
IndRest = setdiff(INDEXsamTRAJ,ipoints,'stable') ;
% Local search for candidates
TolStopSearch = DATAloc.TolStopSearch ;
if isempty(DATAloc.MaxNumberNeighbors)
    MaxNumberNeighbors = length(IndRest) ;
else
    MaxNumberNeighbors = min(DATAloc.MaxNumberNeighbors,length(IndRest)) ;
end
[indCANDIDATES_self,sortedDISTANCES,indCANDIDATESall ]=  LocalSearchNeigh_Clustering(ipoints,d,nd,SNAPdisp,IndRest,TolStopSearch,...
    MaxNumberNeighbors,DISTANCES_POINTS) ;

% ---------------------------------------------------
% Search over points belonging to the rest of trajectories
% ---------------------------------------------------
% Minimum inter-trajectory distance
IndRest = find(WHICH_TRAJ_IS_EACH_SNAP ~= itraj) ;
TolStopSearch = DATAloc.TolStopSearch_interTRAJ ;
if isempty(DATAloc.MaxNumberNeighbors_interTRAJ)
    MaxNumberNeighbors = length(IndRest) ;
else
    MaxNumberNeighbors = min(DATAloc.MaxNumberNeighbors_interTRAJ,length(IndRest)) ;
end

if DATAloc.InterTrajectorySearchWithTwoSelfPoints == 1
    error('This option turned out to be unreliable')
    % Include in the inter-trajectory search points of the same trajectory
    % as point "ipoints" (located at a distance greater than the closest inter-trajectory point)
    if ~ isempty(DISTANCES_POINTS)
        radius_MIN_DIST_loc  =  min(DISTANCES_POINTS(IndRest,ipoints));
    else
        error('MEthodology only available when DISTANCES_POINTS has been previously calculated')
    end
    indLocalSELF = find(sortedDISTANCES >= radius_MIN_DIST_loc) ;
    indLocalSELF = indLocalSELF(1:2) ;
    indGlobalSelf = indCANDIDATESall(indLocalSELF) ;
    
    if isempty(DATAloc.MinNumberNeighbors_interTRAJ)
        MinNumberNeighbors_interTRAJ =0 ;   
    else
        MinNumberNeighbors_interTRAJ = DATAloc.MinNumberNeighbors_interTRAJ ; 
    end
    
    
    indCANDIDATES_rest =  LocalSearchNeigh_ClusteringMIXED(ipoints,d,nd,SNAPdisp,IndRest,TolStopSearch,...
        MaxNumberNeighbors,DISTANCES_POINTS,indGlobalSelf,MinNumberNeighbors_interTRAJ) ;
    
    
else
    
    
    indCANDIDATES_rest =  LocalSearchNeigh_Clustering(ipoints,d,nd,SNAPdisp,IndRest,TolStopSearch,...
        MaxNumberNeighbors,DISTANCES_POINTS) ;
end

%%%%






%indCANDIDATES = [indCANDIDATES_self;indCANDIDATES_rest] ;


% DATAloc = DefaultField(DATAloc,'TolStopSearch_interTRAJ',1e-4) ;  % Search tolerance
% DATAloc = DefaultField(DATAloc,'MaxNumberNeighbors_interTRAJ',10) ;  %  Max. Number neig.
