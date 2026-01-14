function DistanceCentroid = DistanceCentroidsFUN(CENT,NEIGHBORING_CLUSTERS,BasisU_cluster,DATAoffline,DOFl)

if nargin == 0
    load('tmp.mat')
end


nclusters = size(CENT,2) ;
DistanceCentroid.C2 = zeros(nclusters,1) ;

if DATAoffline.EfficentMemoryOption==1
    
    % JAHO, introduced July 24th 2022 
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/09_HOMOG_2D_1T/MakingItEfficient.mlx
    
     % DistanceCentroid.Ct_BasisU = cell(nclusters,nclusters) ;
     NumberOfNeighbors = cellfun(@length,NEIGHBORING_CLUSTERS) ; 
     MaxNumberNeigh = max(NumberOfNeighbors) ; 
     DistanceCentroid.Ct_BasisU = cell(nclusters,MaxNumberNeigh+1) ;
    for icluster = 1:nclusters
        DistanceCentroid.C2(icluster)  = sum(CENT(DOFl,icluster).^2) ;  % W
        ipointsLOC = [icluster,NEIGHBORING_CLUSTERS{icluster}] ; % The first cluster is the current  cluster itself.
        for jclusterLOC = 1:length(ipointsLOC)
            jcluster = ipointsLOC(jclusterLOC) ;
            DistanceCentroid.Ct_BasisU{icluster,jclusterLOC}  =  CENT(DOFl,jcluster)'*BasisU_cluster{icluster} ;
        end
    end
    
    DistanceCentroid.NumberOfNeighbors = NumberOfNeighbors ; 
    
else
    
    DistanceCentroid.Ct_BasisU = cell(nclusters,nclusters) ;
    for icluster = 1:nclusters
        DistanceCentroid.C2(icluster)  = sum(CENT(DOFl,icluster).^2) ;  % W
        ipointsLOC = [icluster,NEIGHBORING_CLUSTERS{icluster}] ; % The first cluster is the current  cluster itself.
        for jclusterLOC = 1:length(ipointsLOC)
            jcluster = ipointsLOC(jclusterLOC) ;
            DistanceCentroid.Ct_BasisU{icluster,jcluster}  =  CENT(DOFl,jcluster)'*BasisU_cluster{icluster} ;
        end
    end
    
    
end

DistanceCentroid.NEIGHBORING_CLUSTERS = NEIGHBORING_CLUSTERS;
