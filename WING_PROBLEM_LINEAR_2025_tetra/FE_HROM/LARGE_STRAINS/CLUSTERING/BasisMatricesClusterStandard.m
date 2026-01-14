function [BasisU_cluster,DISTANCES,IND,CENT,DISTANCE_CENT,errorCLUSTER,DATAoffline,CENTROIDS_REAL,...
    SNAPSHOTS_PER_CLUSTER_after_overlap] = BasisMatricesClusterStandard(DATAoffline,SNAPdisp,DISP_CONDITIONS)

format long g
DATAoffline =  DefaultField(DATAoffline,'PROVIDE_INITIAL_CLUSTERS_MAXIMALLY_INDEPENDENT',0) ;
DATAoffline =  DefaultField(DATAoffline,'DistanceClustering','sqeuclidean') ;
DATAoffline =  DefaultField(DATAoffline,'RANDOM_INITIALIZATION_KMEANS',0) ;



DATAoffline.InitialclusterINdex = [] ;

% if DATAoffline.PROVIDE_INITIAL_CLUSTERS_MAXIMALLY_INDEPENDENT ==1
%
%     INDICES_snapshots = DEIMfun_SELECT(SNAPdisp,DATAoffline) ;
%     [IND,CENT,~,DISTANCE_CENT] = kmeans(SNAPdisp',DATAoffline.NCLUSTERS_BASIS_DISP,...
%         'start',SNAPdisp(:,INDICES_snapshots)','Distance',DATAoffline.DistanceClustering) ;
%
%     CENTROIDS_REAL = CENT ;
%
% else
if  DATAoffline.RANDOM_INITIALIZATION_KMEANS ==0
    rng default  % To avoid producing different result in each run
end
INDEX_SELECT  = 1:size(SNAPdisp,2) ;
switch    DATAoffline.DistanceClustering
    case 'cosine'
        % Remove zero points
        NNN = sum(SNAPdisp.^2,1) ;
        ddddd =  find(NNN<1e-16) ;
        INDEX_SELECT(ddddd) = [] ;
        %    SNAPdisp(:,ddddd) = [] ;
        NNNdiff = diff(ddddd) ; % Provisional, modify this
        DATAoffline.InitialclusterINdex = [2:NNNdiff(1):size(SNAPdisp,2)] ;
        % otherwise
        
end

if DATAoffline.PROVIDE_INITIAL_CLUSTERS_MAXIMALLY_INDEPENDENT ==1
    error('This option is not reliable')
    INDICES_snapshots = DEIMfun_SELECT(SNAPdisp(:,INDEX_SELECT),DATAoffline) ;
    [IND_loc,CENT,~,DISTANCE_CENT] = kmeans(SNAPdisp(:,INDEX_SELECT)',DATAoffline.NCLUSTERS_BASIS_DISP,...
        'start',SNAPdisp(:,INDEX_SELECT(INDICES_snapshots))','Distance',DATAoffline.DistanceClustering) ;
else
    
    [IND_loc,CENT,~,DISTANCE_CENT] = kmeans(SNAPdisp(:,INDEX_SELECT)',DATAoffline.NCLUSTERS_BASIS_DISP,...
        'Distance',DATAoffline.DistanceClustering) ;
    
end

% Inter-cluster connectivity based on DISTANCE_CENT
[DATAoffline,IND_loc,SNAPSHOTS_PER_CLUSTER_after_overlap] = OverlappingClusters(DATAoffline,DISTANCE_CENT,IND_loc,INDEX_SELECT) ; 

 

%  switch    DATAoffline.DistanceClustering
%    case 'cosine'
IND =  ones(size(SNAPdisp,2),1) ;
IND(INDEX_SELECT) = IND_loc ;  %INDEX_SELECT in the cosine approach only includes those snapshots with nonzero norm 

[aaa,IND_CENT] = min(DISTANCE_CENT) ;

CENTROIDS_REAL = SNAPdisp(:,IND_CENT) ;

% end
NSNAP_per_cluster = zeros(size(CENT,1),1) ;
for iii= 1:size(CENT,1)
    NSNAP_per_cluster(iii) = length(find(IND_loc==iii)) ;
end

figure(902)
hold on
hh= bar(NSNAP_per_cluster)
xlabel('Cluster')
ylabel('Number of snapshots (nonzero)')

[aaaa,NSNAP_per_cluster_overla ] = cellfun(@size,SNAPSHOTS_PER_CLUSTER_after_overlap) ; 

 hhh= bar(NSNAP_per_cluster_overla)

legend([hh,hhh],{'NO overlap','Overlap'})

DATAoffline.InitialclusterINdex = IND(DATAoffline.InitialclusterINdex )   ; %%%%s
%end





% LOOP  OVER CLUSTERS TO COMPUTE the  BASIS MATRICES
% --------------------------------------------
BasisU_cluster = cell(1,size(CENT,1)) ;
errorCLUSTER = zeros(1,size(CENT,1)) ;
DISTANCES = zeros(size(CENT,1)) ;

WITH_ALL_DOFS = DATAoffline.CLUSTER.WITH_ALL_DOFS  ;
%
if  WITH_ALL_DOFS == 1
    DOFSincluded = 1:size(SNAPdisp,1) ;
else
    DOFSincluded = DISP_CONDITIONS.DOFl ;
    
end

for i=1:size(CENT,1)
    
   % STEPSCLUSTER = find(IND == i) ;
    
    STEPSCLUSTER = SNAPSHOTS_PER_CLUSTER_after_overlap{i}  ; 
    
    % [U,S,V,eSVD,Rsup] = RSVDT(A,e0,mu,R,DATA) ;
    e0 = DATAoffline.errorDISP_cluster ;
    mu = 0 ; R = 0 ; DATALOC.RELATIVE_SVD = 1;
    DATALOC.HIDE_OUTPUT = 1 ;
    [BasisU_cluster{i},SS,VV,errorCLUSTER(i)] =   RSVDT(SNAPdisp(DOFSincluded,STEPSCLUSTER),e0,mu,R,DATALOC) ;
    
    disp(['SVD cluster =',num2str(i),'; Nmodes =',num2str(length(SS))]) ;
    
    for j = 1:size(CENT,1)
        DISTANCES(i,j) = norm(CENT(i,:)-CENT(j,:));
    end
    
    
end

figure(903)
hold on
xlabel('Cluster')
ylabel('Number of disp. modes')
[~,nnn_modes_cluster]  =cellfun(@size,BasisU_cluster) ;
bar(nnn_modes_cluster) ;

DATAoffline = DefaultField(DATAoffline,'MinimumNumberOfModesPerCluster',3)  ;

if any(nnn_modes_cluster <3)
    error(['The number of modes per cluster should be higher than  =',num2str(DATAoffline.MinimumNumberOfModesPerCluster)])
end



DATAoffline =  DefaultField(DATAoffline,'ExploreIsoMap',0) ;


if DATAoffline.ExploreIsoMap == 1
    
    ExploreIsoMapForClustering(SNAPdisp,DATAoffline) ;
end