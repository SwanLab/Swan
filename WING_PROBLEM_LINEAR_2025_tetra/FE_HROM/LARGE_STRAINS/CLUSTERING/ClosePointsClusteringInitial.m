function [SNAPSHOTS_PER_CLUSTER_after_overlap,CENT,IndClusters_all ]= ...
    ClosePointsClusteringInitial(SNAPdisp,DATAloc,InitialClusterDATA,trajectoryEACHsnap,DATAoffline)

if nargin == 0
    load('tmp.mat')
    %     DATAloc.MaxNumberNeighbors_interTRAJ = [] ;
    %     DATAloc.TolStopSearch_interTRAJ  =1e-3 ;
end

%  DATAloc.MaxNumberNeighbors = [];  % Maximum number of points for stopping searching process.
%                                                       % Empty == Until                                                    % achieving
%                                                       % convergence
%
% DATAloc.UseInterTrajectoryApproach = 1;
% DATAloc.TolStopSearch_interTRAJ = 1e-4;   % Inter-trajectory approach
% DATAloc.MaxNumberNeighbors_interTRAJ = 6;

DATAloc = DefaultField(DATAloc,'TolStopSearch',1e-5) ;   % Search tolerance ALL points or TRAJEC. POINTS
DATAloc = DefaultField(DATAloc,'MaxNumberNeighbors',[]) ;  % Max. Number neig. ALL points or TRAJEC. POINTS

DATAloc = DefaultField(DATAloc,'UseInterTrajectoryApproach',0) ;  % Enable inter-trajectory approach

DATAloc = DefaultField(DATAloc,'TolStopSearch_interTRAJ',1e-4) ;  % Search tolerance
DATAloc = DefaultField(DATAloc,'MaxNumberNeighbors_interTRAJ',[]) ;  %  Max. Number neig.

DATAloc = DefaultField(DATAloc,'InterTrajectorySearchWithTwoSelfPoints',0) ;

DATAloc = DefaultField(DATAloc,'PrecomputeDistancePoints',1) ;  %  Max. Number neig.

disp('-----------------------------------')
disp('Distance between all points')
disp('-----------------------------------')

if DATAloc.PrecomputeDistancePoints == 1
    DISTANCES_POINTS = zeros(size(SNAPdisp,2)) ;
    for ipoints = 1:size(SNAPdisp,2)
        d = SNAPdisp(:,ipoints) ;
        DIST = bsxfun(@minus,SNAPdisp,d) ;
        DISTANCES_POINTS(:,ipoints)  = sqrt(sum(DIST.^2,1)) ;
    end
end

% Euclidean norm snapshots
normSNAPdisp = sqrt(sum(SNAPdisp.^2,1)) ;
% ------------------------------------------------
% DETERMINATION OF INITIAL CLUSTERS (POINTS AROUND ZERO)
% ------------------------------------------------------------
% switch   DATAoffline.METHOD_PARTITION_SNAPSHOTS
%     case 'ClosestPoints_KMEANS'
%         IndClustersZERO = [] ;
%     otherwise

nstart = InitialClusterDATA.NumberOfSnapshotsEachTrajectory  ;  % Number of initial snapshots  (prescribed by the user)
INITIAL_snap = cell(size(trajectoryEACHsnap)) ;
for itraj = 1:length(trajectoryEACHsnap)
    INITIAL_snap{itraj} = trajectoryEACHsnap{itraj}(1:nstart) ;
end
INITIAL_snap = cell2mat(INITIAL_snap) ;

% --- Select all snapshots in a radius -->  max(normSNAPdisp(INITIAL_snap))
normINI = max(normSNAPdisp(INITIAL_snap)) ;
IndClustersZERO = find(normSNAPdisp<=normINI) ;  %  Indexes snapshots associated to the initial cluster
%end
%
% Indexes remaining snapshots = ADDITIONAL CLUSTERS (one point, one cluster = centroid )
IndClusters = setdiff(1:length(normSNAPdisp),IndClustersZERO) ;  %

%% Centroid matrix
CENT = SNAPdisp(:,IndClusters);
%if ~isempty(IndClustersZERO)
CENT = [zeros(size(CENT,1),1),CENT] ;  % The first column is the "zero" centroid
%end

% Indexes snapshots for each cluster
SNAPSHOTS_PER_CLUSTER_after_overlap = cell(1,size(CENT,2)) ;
%if  ~isempty(IndClustersZERO)
SNAPSHOTS_PER_CLUSTER_after_overlap{1} = IndClustersZERO' ;  % First cluster is the zero cluster
IndClusters_all = [1,IndClusters]  ; % Indexes corresponding to all clusters (t)
WHICH_CLUST_IS_EACH_SNAP(IndClusters)  =2:length(IndClusters)+1 ;

% else
%     IndClusters_all = IndClusters ;
%     WHICH_CLUST_IS_EACH_SNAP   =1:length(IndClusters)  ;
%
% end


[WHICH_TRAJ_IS_EACH_SNAP ]=  GetTrajectoryEach_Snapshot(trajectoryEACHsnap,size(SNAPdisp,2)) ;

PLOT_TRAJ_NUMBER = 0 ;
if PLOT_TRAJ_NUMBER == 1
    NumberOfPointsPerTraj.SELF = cell(size(trajectoryEACHsnap)) ;
    NumberOfPointsPerTraj.INTER = cell(size(trajectoryEACHsnap)) ;
end


for ipointsLOC = 1:length(IndClusters)
    
    ipoints = IndClusters(ipointsLOC) ;  % Snapshot index associated to local cluster ipointsLOC
    % itraj = WHICH_TRAJ_IS_EACH_SNAP(ipoints) ;
    
    d = SNAPdisp(:,ipoints) ;
    nd = norm(d) ;
    fprintf(1,'Neigh. of point = %d\n',ipoints) ;
    
    if DATAloc.UseInterTrajectoryApproach == 0 || length(trajectoryEACHsnap) == 1
        % ---------------------
        % REMAINING SNAPSHOTS
        % ---------------------
        IndRest = setdiff(1:size(SNAPdisp,2),ipoints,'stable') ;
        % Local search for candidates
        TolStopSearch = DATAloc.TolStopSearch ;
        if isempty(DATAloc.MaxNumberNeighbors)
            MaxNumberNeighbors = length(IndRest) ;
        else
            MaxNumberNeighbors = min(DATAloc.MaxNumberNeighbors,length(IndRest)) ;
        end
        indCANDIDATES =  LocalSearchNeigh_Clustering(ipoints,d,nd,SNAPdisp,IndRest,TolStopSearch,MaxNumberNeighbors,DISTANCES_POINTS) ;
    else
        [indCANDIDATES_self,indCANDIDATES_rest]  =...
            InterTrajLocal_Search(d,nd,SNAPdisp,DATAloc,ipoints,WHICH_TRAJ_IS_EACH_SNAP,...
            IndClusters,DISTANCES_POINTS)   ;
        indCANDIDATES = [indCANDIDATES_self;indCANDIDATES_rest]  ;
        if  PLOT_TRAJ_NUMBER == 1
            ITRAJ = WHICH_TRAJ_IS_EACH_SNAP(ipoints) ;
            NumberOfPointsPerTraj.SELF{ITRAJ}(end+1) = length(indCANDIDATES_self)  ;
            NumberOfPointsPerTraj.INTER{ITRAJ}(end+1) =   length(indCANDIDATES_rest)  ;
        end
    end
    
    
    SNAPSHOTS_PER_CLUSTER_after_overlap{ipointsLOC+1} = [indCANDIDATES(:)',ipoints]' ;
    
end
DATAoffline = DefaultField(DATAoffline,'RANDOM_INITIALIZATION_KMEANS',0)  ;
switch   DATAoffline.METHOD_PARTITION_SNAPSHOTS
    case 'ClosestPoints_KMEANS'
        % Re-defining both CENT (matrix of centroids) and SNAPSHOTS_PER_CLUSTER_after_overlap
        ClusterIni = SNAPSHOTS_PER_CLUSTER_after_overlap{1} ;
        if  DATAoffline.RANDOM_INITIALIZATION_KMEANS ==0
            rng default  % To avoid producing different result in each run
        end
        DistanceClustering = 'sqeuclidean' ;
        
        if isempty(DATAoffline.NCLUSTERS_BASIS_DISP)
            DATAoffline.NCLUSTERS_BASIS_DISP = length(IndClusters) ; 
        end
        
        [IND_loc_kmeans,CENT_kmeans,~,DISTANCE_CENT] = kmeans(SNAPdisp(:,IndClusters)',DATAoffline.NCLUSTERS_BASIS_DISP,...
            'Distance', DistanceClustering) ;
        % IND_loc_kmeans refer to  IndClusters ... We need the global
        % counterpart of  IND_loc_kmeans: IND_glo_kmeans
        %  IND_glo_kmeans =  zeros(size())  IndClusters(IND_loc_kmeans) ;
        
        % Centroides
        CENT = [zeros(size(CENT_kmeans,2),1),CENT_kmeans'] ; % Initial cluster has the centroid at  ZERO
        
        % Overlapping
        % -----------
        SNAPSHOT_per_clusterNEW = cell(DATAoffline.NCLUSTERS_BASIS_DISP+1,1) ; % New clusters
        SNAPSHOT_per_clusterNEW{1} = ClusterIni ;  % The first one is the ZERO cluster
        disp('Loop over KMEANS cluster')
        for iclusterLOC = 1:DATAoffline.NCLUSTERS_BASIS_DISP
            
            IndPointsKMEANS = find(IND_loc_kmeans == iclusterLOC) ;  % Indexes snapshots cluster iclusterLOC
            IndPointsKMEANS = IndClusters(IndPointsKMEANS) ;
            
            ALL_LOC_CLUSTER_OLD = WHICH_CLUST_IS_EACH_SNAP(IndPointsKMEANS) ; % Indexes of the corresponding noninitial-points
            ALL_LOC_CLUSTER_OLD = ALL_LOC_CLUSTER_OLD(ALL_LOC_CLUSTER_OLD>0) ;
            
            ALL_P =  SNAPSHOTS_PER_CLUSTER_after_overlap(ALL_LOC_CLUSTER_OLD) ;
            ALL_P = unique(cell2mat(ALL_P')) ;  % CLuster points + neighboring points
            
            SNAPSHOT_per_clusterNEW{iclusterLOC+1} = unique([IndPointsKMEANS(:);ALL_P(:)]) ;   % Final cluster
            INDEX_OVERLAP = length(SNAPSHOT_per_clusterNEW{iclusterLOC+1})/length(IndPointsKMEANS)*100-100 ;
            
            disp(['Kmeans cl =',num2str(iclusterLOC),'  OVERLAPPING = ',num2str(INDEX_OVERLAP),' %']);
            
            
            
        end
        
        
end

SNAPSHOTS_PER_CLUSTER_after_overlap = SNAPSHOT_per_clusterNEW ;


if PLOT_TRAJ_NUMBER == 1
    figure(167)
    hold on
    itraj = 16;
    title(['Nearest Points trajectory ',num2str(itraj)])
    xlabel('Index Point/Cluster')
    ylabel('Number Neigh. points')
    h1 =  plot(NumberOfPointsPerTraj.SELF{itraj})  ;
    h2 =  plot(NumberOfPointsPerTraj.INTER{itraj})  ;
    legend([h1,h2],{'Self','Inter'})
end

