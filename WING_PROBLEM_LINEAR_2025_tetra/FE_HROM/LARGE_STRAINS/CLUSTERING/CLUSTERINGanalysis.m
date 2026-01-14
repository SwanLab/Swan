function [BasisU_cluster,IND,CENT,SNAPdisp,TransitionCLUST,IND_input,DATAoffline,SNAPSHOTS_PER_CLUSTER_after_overlap,DELAUNAY_PARAM_SPACE ]= ...
    CLUSTERINGanalysis(SNAPdisp_traj,DATAoffline,MESH,TIME_STEPS,MATRIX_POINTS_SPACE_PARAMETER,BasisU_ALL,TESTING_PARAM,...
    INDEXES_TRAJECTORIES,DISP_CONDITIONS,NAME_INPUTS)

if nargin == 0
    load('tmp.mat')
    %   DATAoffline.NCLUSTERS_BASIS_DISP = 8;
    %     DATAoffline.PLOT_CLUSTERS  = 0 ;
    %     close all
end
%IND_STUDY = 1;

DATAoffline = DefaultField(DATAoffline,'EfficentMemoryOption',1) ; 


DATAoffline.CLUSTER.WITH_ALL_DOFS =   1; %PRELIMINARY_ANALYSIS = 1;   = 0;
DATAoffline = DefaultField(DATAoffline,'CLUSTER',[]) ;
DATAoffline = DefaultField(DATAoffline,'PRELIMINARY_ANALYSIS',0) ;
DATAoffline.CLUSTER.WITH_ALL_DOFS = DATAoffline.PRELIMINARY_ANALYSIS ;
if  DATAoffline.PRELIMINARY_ANALYSIS ==1
    disp(['******************************************************++'])
    disp(['PRELIMINARY (A PRIORI ) CLUSTERING ANALYSIS '])
    disp(['******************************************************++'])
    
end

DELAUNAY_PARAM_SPACE = [] ;

WITH_ALL_DOFS =   1;

DATAoffline = DefaultField(DATAoffline,'TESTING_TRAJ_IS_TRAINING_TRAJ',[]) ;
if ~isempty(DATAoffline.TESTING_TRAJ_IS_TRAINING_TRAJ)
    TESTING_SNAPSHOTS = SNAPdisp_traj{DATAoffline.TESTING_TRAJ_IS_TRAINING_TRAJ} ;
else
    TESTING_SNAPSHOTS = [] ;
end

% INDEX_TRAIN_TRAJ = cell(size(SNAPdisp_traj)) ;
% iacum = 1;
% for iii = 1:length()
% end

%[~,DATAoffline.Ind_Trajectories] = cellfun(@size,SNAPdisp_traj) ;
% IND_traj = cell(size(SNAPdisp_traj)) ;
% for itraj = 1:length(SNAPdisp_traj)
% end
trajectoryEACHsnap = [] ; 
if isstruct(SNAPdisp_traj)
    SNAPdisp =  (SNAPdisp_traj.DOFl.coeff) ;  % This method work on a reduced space 
    trajectoryEACHsnap = SNAPdisp_traj.DOFl.nrowsPROJECT ; 
else
    SNAPdisp= cell2mat(SNAPdisp_traj ) ; 
end

% for idim = 1:size(MESH.COOR,2)
%    IND = idim:size(MESH.COOR,2):size(SNAPdisp,1);
%
%    SNAPdisp(IND,:) = bsxfun(@plus,SNAPdisp(IND,:),MESH.COOR(:,idim));
%
% end

%DATAoffline = DefaultField(DATAoffline,'HierarchicalClustering',0); % No longer available
DATAoffline = DefaultField(DATAoffline,'UseSparseSubspaceClustering',0);
DATAoffline = DefaultField(DATAoffline,'METHOD_PARTITION_SNAPSHOTS','PartitionStateManifold');


CENTROIDS_REAL = [] ;DELAUNAY_PARAM_SPACE = [] ;
 DistanceCentroid = [] ; DistanceCentroid.NEIGHBORING_CLUSTERS =[] ;  
%if  DATAoffline.HierarchicalClustering == 0
switch DATAoffline.METHOD_PARTITION_SNAPSHOTS
    case 'PartitionStateManifold'
        if DATAoffline.UseSparseSubspaceClustering ==1
            error('This option has not been properly examined')
            [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline] = BasisMatricesClusterStandardSparse(DATAoffline,SNAPdisp,DISP_CONDITIONS) ;
        else
            [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline,CENTROIDS_REAL,...
                SNAPSHOTS_PER_CLUSTER_after_overlap] = BasisMatricesClusterStandard(DATAoffline,SNAPdisp,DISP_CONDITIONS) ;
            CentroidCoordinates = [] ;
            IND_input = IND ;
        end
    case 'PartitionParameterDomain'
        [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline,CENTROIDS_REAL,...
            SNAPSHOTS_PER_CLUSTER_after_overlap,DELAUNAY_PARAM_SPACE] = ...
            BasisMatricesParameterSpace(DATAoffline,SNAPdisp,DISP_CONDITIONS,MATRIX_POINTS_SPACE_PARAMETER) ;
        CentroidCoordinates = [] ;
        IND_input = IND ;
    case {'MaximumClusteringClosestPoints','ClosestPoints_KMEANS' }
        [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline,CENTROIDS_REAL,...
            SNAPSHOTS_PER_CLUSTER_after_overlap,DELAUNAY_PARAM_SPACE,DistanceCentroid] = ...
            BasisMatricesMaxClClosePoints(DATAoffline,SNAPdisp,DISP_CONDITIONS,MATRIX_POINTS_SPACE_PARAMETER,BasisU_ALL,...
            trajectoryEACHsnap) ;
        CentroidCoordinates = [] ;
        IND_input = IND ;
        DATAoffline.DistanceClustering = 'euclidean' ; 
        
    case 'SEQUENTIAL_unitrajectory_3'    
        disp('--------------------------------------------------')
        disp('Option only valid for one-trajectory problems ')
        disp('---------------------------------------------------')
         [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline,CENTROIDS_REAL,...
            SNAPSHOTS_PER_CLUSTER_after_overlap,DELAUNAY_PARAM_SPACE,DistanceCentroid] = ...
            SEQUENTIAL_unitrajectory_3_fun(DATAoffline,SNAPdisp,DISP_CONDITIONS,MATRIX_POINTS_SPACE_PARAMETER,BasisU_ALL,...
            trajectoryEACHsnap) ;
         CentroidCoordinates = [] ;
        IND_input = IND ;
        DATAoffline.DistanceClustering = 'euclidean' ; 
        
    otherwise
        error('Option not implemented')
end
% else
%         error('Abandon this approach... It is not useful...')
%     [BasisU_cluster,IND,CentroidIndex,CENT]  = ClusterViaHierarchicalKmeans(SNAPdisp,DATAoffline,DISP_CONDITIONS) ;
%     IND_input = zeros(size(SNAPdisp,2),1) ;
%     for index_cluster = 1:length(IND)
%         IND_input(IND{index_cluster}) = index_cluster ;
%     end
%
%     DATAoffline.CentroidIndex = CentroidIndex ;
%     errorCLUSTER = [] ;
%     DISTANCES = zeros(size(CENT,2),size(CENT,2));
%     for icluster = 1:size(CENT,2)
%         for jcluster = 1:size(CENT,2)
%             DISTANCES(icluster,jcluster) = norm(CENT(:,icluster)-CENT(:,jcluster)) ;
%         end
%     end
%     DISTANCE_CENT = [] ;
%     CENT = CENT' ;
% end

% Display table with distances

if ~isempty(DISTANCES)
    NameColums = TableDistances(errorCLUSTER,SNAPdisp,DISTANCES) ;
end
% Conncectivites, Intersection  subpaces
[TransMatrix,CONNECTIVITY_CLUSTERS,LOWERBOUND_matrix,UPPERBOUND_matrix] = ...
    SubspacesIntersectionClusters(BasisU_cluster,DATAoffline,SNAPdisp,BasisU_ALL,DistanceCentroid) ;
TransitionCLUST.TransMatrix = TransMatrix ;
TransitionCLUST.LOWERBOUND_matrix = LOWERBOUND_matrix ;
TransitionCLUST.UPPERBOUND_matrix = UPPERBOUND_matrix ;
 TransitionCLUST.NEIGHBORING_CLUSTERS = DistanceCentroid.NEIGHBORING_CLUSTERS; 
 DistanceCentroid.NEIGHBORING_CLUSTERS = [];
 TransitionCLUST.DistanceCentroid = DistanceCentroid; 



if DATAoffline.CLUSTER.WITH_ALL_DOFS  == 1
    % ONLY WHEN   DATAoffline.PRELIMINARY_ANALYSIS  =1
    % Plot parametric space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    DATAoffline = DefaultField(DATAoffline,'PLOT_CLUSTERS_parameter_space',1) ;
    % Plot Deformed shapes
    % PLITTING CLUSTERSS (DEFORMED SHAPES)
    DATAoffline = DefaultField(DATAoffline,'PLOT_CLUSTERS_DEFORMED_SHAPE',1) ;
    DATAoffline = DefaultField(DATAoffline,'PLOT_CLUSTERS_DEFORMED_SHAPE_ONLY_CENTROIDS',1) ;
    DATAoffline = DefaultField(DATAoffline,'PLOT_ONLY_deformed_shapes_GIVEN_SNAPSHOTS',[]) ;
    
    DATAoffline.COLOR_plot_def = 'r';
    DATAoffline.LINEWIDTH_plot_def = 3;
    
    if  ~isempty(CENTROIDS_REAL)
        DATAoffline =  PlotDeformedShapeClusters(DATAoffline,IND,MESH,CENTROIDS_REAL,SNAPdisp) ;
    end
    
    DATAoffline = DefaultField(DATAoffline,'PlotDistanceVersusParameterMapping',0);
    
    if DATAoffline.PlotDistanceVersusParameterMapping == 1
        DATAoffline =   PlotMappingSpaceParameterDistance(DATAoffline,MATRIX_POINTS_SPACE_PARAMETER,INDEXES_TRAJECTORIES,SNAPdisp)
    end
    
    % CHECKING TESTING TRAJECTORY
    IdealIndexCluster =  PlotParametricSpaceClusters(DATAoffline,MATRIX_POINTS_SPACE_PARAMETER,IND,DISTANCE_CENT,BasisU_cluster,...
        LOWERBOUND_matrix,TransMatrix,UPPERBOUND_matrix,CENT,BasisU_ALL,TESTING_SNAPSHOTS,TESTING_PARAM,MESH,NAME_INPUTS) ;
    %
    %
    %     NMAME = ['Ideal_',NAME_INPUTS,'.mat'] ;
    %     save(NMAME,'IdealIndexCluster')
    
    
    %%%% MONITORING TRAJECTORY OF THE DISPLACEMENT VECTOR WITHIN EACH cluster
    disp('***************************************************************************************************')
    disp(['Analyzing the trajectory of the displacement state within each cluster (for each loading state)'])
    disp('***************************************************************************************************')
    
    for itrajectory = 1:length(SNAPdisp_traj)
        error('Not finished yet-....')
        STEPS =  TIME_STEPS{itrajectory} ;
        disp(['Loading state  =  ',num2str(itrajectory)])
        disp('************************************************')
        
        isnap = 1;
        indCLUSTER = 1;
        d = SNAPdisp_traj{itrajectory}(:,isnap);  % current state
        nsteps = size(SNAPdisp_traj{itrajectory},2) ;
        % In the initial snapshot, the determination of the cluster by
        % means of the distance
        DISTANCEScent = zeros(size(CENT,1),1) ;
        for icentroids = 1:size(CENT,1)
            DISTANCEScent(icentroids) = norm(CENT(icentroids,:)'-d);
        end
        [~,kCLUSTER] = min(DISTANCEScent)  ;
        disp(['Initial cluster = ',num2str(kCLUSTER)]) ;
        disp(['Neighboring clusters = ',num2str(CONNECTIVITY_CLUSTERS{kCLUSTER})]) ;
        
        % Plotting aroung the vicinity of the current cluseter
        PLOT_NEIG = 1 ;
        if PLOT_NEIG == 1
            MAX_snap = isnap + 100 ;
            ClusterPlotLocalIni(kCLUSTER,isnap,MAX_snap,BasisU_cluster,SNAPdisp_traj{itrajectory},CONNECTIVITY_CLUSTERS) ;
        end
        %  isnapINI =
        isnap = isnap + 1;
        while  isnap <= length(STEPS)
            %
            d = SNAPdisp_traj{itrajectory}(:,isnap); % Full displacements
            qCL=  BasisU_cluster{kCLUSTER}'*SNAPdisp_traj{itrajectory}(:,isnap);  % Reduced displacements
            
            % Loop over neighboring clusters
            for ineigCLUSTER = 1:length(CONNECTIVITY_CLUSTERS{kCLUSTER})
                gCLUSTER= CONNECTIVITY_CLUSTERS{kCLUSTER}(ineigCLUSTER) ;
                % Intersection space
                qCL_int = IntersectionMatrices{kCLUSTER,gCLUSTER}'*qCL ;
                norm_qCL_int= sqrt(sum(qCL_int.^2,1));
                
                alphaINT = norm_qCL_int./norm_qCL ;
                
                hhh(ineigCLUSTER) = plot(REMAIN_steps,alphaINT,'--') ;
                LLL{ineigCLUSTER} =['clust = ',num2str(gCLUSTER)] ;
                
                
            end
            
        end
        
        
        
        
        
    end
    
end




