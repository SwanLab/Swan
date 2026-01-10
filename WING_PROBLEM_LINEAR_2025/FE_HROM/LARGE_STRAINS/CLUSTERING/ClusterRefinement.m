function [BasisU_cluster_new,INDICES_CLUSTERScell_new,CentroidIndex_NEW] = ClusterRefinement(A,BasisU_cluster,INDICES_CLUSTERScell,DATAoffline,CentroidIndex)


BasisU_cluster_new = {} ; 
INDICES_CLUSTERScell_new ={} ; 
INDICES_CENTROIDS_NEW = {}  ; 
CentroidIndex_NEW = [] ; 
iacum = 1; 
for icluster = 1:length(BasisU_cluster)
    if size(BasisU_cluster{icluster},2)  > DATAoffline.nminBASIS
        disp(['Cluster =',num2str(icluster),'; N. Basis Matrices = ',num2str(size(BasisU_cluster{icluster},2) )])
        % New snapshot matrix
        LOCind = INDICES_CLUSTERScell{icluster}; 
        Aloc1 = A(:,LOCind) ;
        DATAoffline_loc1 = DATAoffline ;
        DATAoffline_loc1.NCLUSTERS_BASIS_DISP = 2;  %
        [BasisU_cluster_loc1,INDICES_CLUSTERScell_loc1,IndLOC1,CentroidIndex1] = ClusteringAndBases(Aloc1,DATAoffline_loc1) ;
        % 4. Level 2.  SEcond refinement
        INDICES_CLUSTERScell{icluster} = cell(1,length(INDICES_CLUSTERScell_loc1)) ; 
        for iiiLOC = 1:length(INDICES_CLUSTERScell_loc1)
            INDICES_CLUSTERScell_new{iacum}  =  LOCind(INDICES_CLUSTERScell_loc1{iiiLOC}) ;  
            BasisU_cluster_new{iacum} = BasisU_cluster_loc1{iiiLOC} ; 
            CentroidIndex_NEW(iacum) = LOCind(CentroidIndex1(iiiLOC)) ; 
            iacum = iacum + 1; 
        end 
        
    else
        BasisU_cluster_new{iacum} = BasisU_cluster{icluster} ; 
        INDICES_CLUSTERScell_new{iacum} = INDICES_CLUSTERScell{icluster} ; 
         CentroidIndex_NEW(iacum)  = CentroidIndex(icluster) ; 
        iacum = iacum + 1; 
    end
    
end