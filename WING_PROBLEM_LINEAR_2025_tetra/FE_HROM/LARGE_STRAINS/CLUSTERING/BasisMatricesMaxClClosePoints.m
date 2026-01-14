function [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline,...
    CENTROIDS_REAL,SNAPSHOTS_PER_CLUSTER_after_overlap,DELAUNAY_PARAM_SPACE,DistanceCentroid] =...
    BasisMatricesMaxClClosePoints(DATAoffline,SNAPdisp,DISP_CONDITIONS,MATRIX_POINTS_SPACE_PARAMETER,BasisU_ALL,...
    trajectoryEACHsnap)

if nargin == 0
    load('tmp1.mat')
end

format long g

IND = [] ; CENT = [] ; DISTANCE_CENT = []  ; CENTROIDS_REAL = [] ;  DISTANCES = [] ; SNAPSHOTS_PER_CLUSTER_after_overlap = [] ;
DELAUNAY_PARAM_SPACE = [] ;


DATAloc = DATAoffline.MaxClustCPoints ;


% OUTPUTS REQUIRED FOR CONSTRUCTING THE BASIS MATRICES
% -----------------------------------------------------
DATAoffline = DefaultField(DATAoffline,'InitialClusterDATA',[]) ;
InitialClusterDATA = DATAoffline.InitialClusterDATA ;
InitialClusterDATA = DefaultField(InitialClusterDATA,'NumberOfSnapshotsEachTrajectory',[]) ;

if isempty(InitialClusterDATA.NumberOfSnapshotsEachTrajectory)
    error('Better use the other option (or disable this message if you want to compare their performance)')
    [SNAPSHOTS_PER_CLUSTER_after_overlap,CENT ]= ClosePointsClusteringStandard(SNAPdisp,DATAloc) ;
else
    [SNAPSHOTS_PER_CLUSTER_after_overlap,CENT,IND] =...
        ClosePointsClusteringInitial(SNAPdisp,DATAloc,InitialClusterDATA,trajectoryEACHsnap,DATAoffline) ;
end


STUDY_ACCURACY_SVD_standard = 0 ;

if STUDY_ACCURACY_SVD_standard ==1    
    disp('Accuracy provided by the SVD for these clusters')    
    ERRORloc_standardSVD = zeros(size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;    
    for iclusters = 1:length(SNAPSHOTS_PER_CLUSTER_after_overlap)        
        SNAPloc = SNAPdisp(:,SNAPSHOTS_PER_CLUSTER_after_overlap{iclusters}) ;
        ERRORloc_standardSVD(iclusters) = norm(SNAPloc-BasisU_ALL*(BasisU_ALL'*SNAPloc),'fro')/norm(SNAPloc,'fro') ;        
    end    
    figure(2034)
    hold on
    xlabel('Number of point')
    ylabel('Log10(error)')
    plot(log10(ERRORloc_standardSVD))
    aaa = axis;    
    plot([aaa(1),aaa(2)],log10(DATAloc.TolStopSearch)*ones(1,2),'k','LineWidth',6)  
end


figure(902)
hold on
[NSNAP_per_cluster_overla,aaa ] = cellfun(@size,SNAPSHOTS_PER_CLUSTER_after_overlap) ;
hhh= bar(NSNAP_per_cluster_overla)
DATAoffline.InitialclusterINdex = [] ; % IND(DATAoffline.InitialclusterINdex )   ; %%%%s
%end
xlabel('Index centroids')
ylabel('Number of points per cluster')


% LOOP  OVER CLUSTERS TO COMPUTE the  BASIS MATRICES
% --------------------------------------------
BasisU_cluster = cell( size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;
errorCLUSTER = zeros(size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;
%DISTANCES = zeros(size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;

WITH_ALL_DOFS = DATAoffline.CLUSTER.WITH_ALL_DOFS  ;
%
if  WITH_ALL_DOFS == 1 ||  (DATAoffline.CompressClusterDisplacementStresses==1)
    DOFSincluded = 1:size(SNAPdisp,1) ;
else
    DOFSincluded = DISP_CONDITIONS.DOFl ;
    
end


DATAoffline = DefaultField(DATAoffline,'USE_STANDARD_SVD_FOR_CLUSTERS_PARAM',0) ;


for i=1:length(SNAPSHOTS_PER_CLUSTER_after_overlap)    
    % STEPSCLUSTER = find(IND == i) ;    
    STEPSCLUSTER = SNAPSHOTS_PER_CLUSTER_after_overlap{i}  ;    
    % [U,S,V,eSVD,Rsup] = RSVDT(A,e0,mu,R,DATA) ;
    e0 = DATAoffline.errorDISP_cluster ;
    mu = 0 ; R = 0 ; DATALOC.RELATIVE_SVD = 1;
    DATALOC.HIDE_OUTPUT = 1 ;
    if  DATAoffline.USE_STANDARD_SVD_FOR_CLUSTERS_PARAM == 1
        [BasisU_cluster{i},SS,VV,errorCLUSTER(i)] =   SVDT(SNAPdisp(DOFSincluded,STEPSCLUSTER),e0,DATALOC) ;
    else
        [BasisU_cluster{i},SS,VV,errorCLUSTER(i)] =   RSVDT(SNAPdisp(DOFSincluded,STEPSCLUSTER),e0,mu,R,DATALOC) ;
    end    
    disp(['SVD cluster =',num2str(i),'; Nmodes =',num2str(length(SS))]) ;    
    %     for j = 1:size(CENT,1)
    %         DISTANCES(i,j) = norm(CENT(i,:)-CENT(j,:));
    %     end   
    
end

figure(903)
hold on
xlabel('Cluster')
ylabel('Number of disp. modes')
[~,nnn_modes_cluster]  =cellfun(@size,BasisU_cluster) ;
bar(nnn_modes_cluster) ;

% figure(67)
% hold on
% xlabel('param 1') ; ylabel('param 2') ; zlabel('Number of modes')
% surf(CENTROID_param(:,1),CENTROID_param(:,2),nnn_modes_cluster) ;


DATAoffline = DefaultField(DATAoffline,'MinimumNumberOfModesPerCluster',2)  ;

if any(nnn_modes_cluster <DATAoffline.MinimumNumberOfModesPerCluster)
    error(['The number of modes per cluster should be higher than  =',num2str(DATAoffline.MinimumNumberOfModesPerCluster)])
end


DATAoffline = DefaultField(DATAoffline,'MaximumNumberNeighborsToSearchInterConnectivities',20) ;  
NEIGHBORING_CLUSTERS = ClusterInterConnectivities(SNAPSHOTS_PER_CLUSTER_after_overlap,CENT,DATAoffline) ; 



[aaa,bbb] = cellfun(@size,NEIGHBORING_CLUSTERS)  ;

figure(190)
hold on
xlabel('Cluster/Point')
ylabel('Number of neighboring clusters')
bar(bbb)


%%%% COMPUTING THE quantities needed to
% compute the DISTANCE BETWEEN POINTS IN A FAST FASHION
% Recall that
% IndClusters = find(normSNAPdisp>TOLzero) ;  % IndClusters(j) returns the SNAPSHOT INDEX corresponding to CENTROID j
% We need the inverse, that is, a vector such that IndClustersInverse(j) returns the
% CENTROID index associated with the j-th snapshot
% IndClustersInverse = zeros(npoints,1) ;
% IndClustersInverse(IndClusters)  =[1:length(IndClusters)]' ;


if DATAoffline.CompressClusterDisplacementStresses == 1
    DOFl = 1:size(CENT,1) ;
else
    DOFl =  DISP_CONDITIONS.DOFl  ;
end

% The following variables come into play when computing, in the online
% phase, the intercentroid distance
DistanceCentroid = DistanceCentroidsFUN(CENT,NEIGHBORING_CLUSTERS,BasisU_cluster,DATAoffline,DOFl) ; 






%
%
%
% DATAoffline =  DefaultField(DATAoffline,'ExploreIsoMap',0) ;
%
%
% if DATAoffline.ExploreIsoMap == 1
%
%     ExploreIsoMapForClustering(SNAPdisp,DATAoffline) ;
% end