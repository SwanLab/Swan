function [ECMdata_cluster,setCandidates,LIST_OF_CANDIDATES] = ECM_local_clusteringLOOP(nclusters,NEW_ORDER_clusters,BasisStwo_cluster,OPERFE,DATAoffline,BasisFint_cluster,...
    DATA,BasisU_cluster,BasisPone_cluster)

if nargin == 0
    load('tmp1.mat')
end

setCandidates = [] ;
ECMdata_cluster = cell(1,nclusters) ;

AAA = tic ; 

LIST_OF_CANDIDATES = cell(1,nclusters) ; 
for iclusterLOC = 1:nclusters
    LIST_OF_CANDIDATES{iclusterLOC} = setCandidates ; 
    icluster = NEW_ORDER_clusters(iclusterLOC) ;
    %     if icluster == 553
    %         disp('borrar esto....')
    %     end
    
    disp('*************************************************++')
    disp(['CLUSTER = ',num2str(icluster)])
    disp('*************************************************++')
    if isstruct(BasisStwo_cluster)
        BasisStwo= BasisStwo_cluster.BASIS*BasisStwo_cluster.coeff{icluster} ; ;   % Bssis matrix provided in a "multiscale" fashion
    else
        BasisStwo = BasisStwo_cluster{icluster} ;
    end
    [ECMdata_cluster{icluster},setCandidates ]=...
        ECM_clustersNC(OPERFE.wSTs,DATAoffline,BasisFint_cluster{icluster},DATA,BasisStwo,setCandidates) ;
    disp(['*****************************'])
    
    disp(['*********************************************************************'])
    if isstruct(BasisU_cluster)
        nmodesU = size(BasisU_cluster.coeff{icluster},2) ;
        nmodesS = size(BasisStwo_cluster.coeff{icluster},2) ;
    else
        nmodesU =size(BasisU_cluster{icluster},2) ;
        nmodesS = size(BasisStwo_cluster{icluster},2) ;
    end
    disp(['Number of displacement modes = ',num2str(nmodesU)])
    disp(['Number of PK2-stress modes = ',num2str(nmodesS)])
    disp(['Number of PK1-stress modes = ',num2str(size(BasisPone_cluster{icluster},2))])
end

AAA= toc(AAA) ; 
disp(['Time for the SAW-ECM: ',num2str(AAA),' s'])
