function [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline,...
    CENTROIDS_REAL,SNAPSHOTS_PER_CLUSTER_after_overlap,DELAUNAY_PARAM_SPACE] =...
    BasisMatricesParameterSpace(DATAoffline,SNAPdisp,DISP_CONDITIONS,MATRIX_POINTS_SPACE_PARAMETER)

if nargin == 0
    load('tmp.mat')
end

format long g

IND = [] ; CENT = [] ; DISTANCE_CENT = []  ; CENTROIDS_REAL = [] ;  DISTANCES = [] ;

% Number of input parameters
ninp = size(MATRIX_POINTS_SPACE_PARAMETER,2) ;
nsnap =  size(MATRIX_POINTS_SPACE_PARAMETER,1) ;

DELAUNAY_PARAM_SPACE = [] ;

% Partition parameter domain (default UNIFORM)
DATAoffline = DefaultField(DATAoffline,'PartitionParamType','Delaunay') ; %'UNIFORM_cartesian'
switch DATAoffline.PartitionParamType
    case 'UNIFORM_cartesian'
        error('Option not tested yet')
        [INDEX_CLUSTERS_WITH_OVERLAPPING,INDEX_CLUSTERS_NO_OVERLAPPING] = CartesianPartParamSpace(DATAoffline) ;
    case 'Delaunay'
        [SNAPSHOTS_PER_CLUSTER_after_overlap,INDEX_CLUSTERS_NO_OVERLAPPING,DELAUNAY_PARAM_SPACE] = DelaunayPartParamSpace(DATAoffline,MATRIX_POINTS_SPACE_PARAMETER) ;
        
        
    otherwise
        error('Option not implemented')
end

%
%
%
% if DATAoffline.PROVIDE_INITIAL_CLUSTERS_MAXIMALLY_INDEPENDENT ==1
%     error('This option is not reliable')
%     INDICES_snapshots = DEIMfun_SELECT(SNAPdisp(:,INDEX_SELECT),DATAoffline) ;
%     [IND_loc,CENT,~,DISTANCE_CENT] = kmeans(SNAPdisp(:,INDEX_SELECT)',DATAoffline.NCLUSTERS_BASIS_DISP,...
%         'start',SNAPdisp(:,INDEX_SELECT(INDICES_snapshots))','Distance',DATAoffline.DistanceClustering) ;
% else
%
%     [IND_loc,CENT,~,DISTANCE_CENT] = kmeans(SNAPdisp(:,INDEX_SELECT)',DATAoffline.NCLUSTERS_BASIS_DISP,...
%         'Distance',DATAoffline.DistanceClustering) ;
%
% end
%
% % Inter-cluster connectivity based on DISTANCE_CENT
% [DATAoffline,IND_loc,SNAPSHOTS_PER_CLUSTER_after_overlap] = OverlappingClusters(DATAoffline,DISTANCE_CENT,IND_loc,INDEX_SELECT) ;
%
%
%
% %  switch    DATAoffline.DistanceClustering
% %    case 'cosine'
% IND =  ones(size(SNAPdisp,2),1) ;
% IND(INDEX_SELECT) = IND_loc ;  %INDEX_SELECT in the cosine approach only includes those snapshots with nonzero norm
%
% [aaa,IND_CENT] = min(DISTANCE_CENT) ;
%
% CENTROIDS_REAL = SNAPdisp(:,IND_CENT) ;
%
% % end
% NSNAP_per_cluster = zeros(size(CENT,1),1) ;
% for iii= 1:size(CENT,1)
%     NSNAP_per_cluster(iii) = length(find(IND_loc==iii)) ;
% end

figure(902)
hold on

[NSNAP_per_cluster_nooverla,sewegtsd ] = cellfun(@size,INDEX_CLUSTERS_NO_OVERLAPPING) ;

hh= bar(NSNAP_per_cluster_nooverla)
xlabel('Cluster')
ylabel('Number of snapshots (nonzero)')

[NSNAP_per_cluster_overla,aaa ] = cellfun(@size,SNAPSHOTS_PER_CLUSTER_after_overlap) ;

hhh= bar(NSNAP_per_cluster_overla)

legend([hh,hhh],{'NO overlap','Overlap'})

DATAoffline.InitialclusterINdex = [] ; % IND(DATAoffline.InitialclusterINdex )   ; %%%%s
%end





% LOOP  OVER CLUSTERS TO COMPUTE the  BASIS MATRICES
% --------------------------------------------
BasisU_cluster = cell( size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;
errorCLUSTER = zeros(size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;
%DISTANCES = zeros(size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;

WITH_ALL_DOFS = DATAoffline.CLUSTER.WITH_ALL_DOFS  ;
%
if  WITH_ALL_DOFS == 1
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


DATAoffline = DefaultField(DATAoffline,'MinimumNumberOfModesPerCluster',3)  ;

if any(nnn_modes_cluster <DATAoffline.MinimumNumberOfModesPerCluster)
    error(['The number of modes per cluster should be higher than  =',num2str(DATAoffline.MinimumNumberOfModesPerCluster)])
end







DATAoffline =  DefaultField(DATAoffline,'ExploreIsoMap',0) ;


if DATAoffline.ExploreIsoMap == 1
    
    ExploreIsoMapForClustering(SNAPdisp,DATAoffline) ;
end