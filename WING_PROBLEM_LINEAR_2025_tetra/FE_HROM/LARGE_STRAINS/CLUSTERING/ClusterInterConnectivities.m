function NEIGHBORING_CLUSTERS = ClusterInterConnectivities(SNAPSHOTS_PER_CLUSTER_after_overlap,CENT,DATAoffline)
if nargin == 0
    load('tmp_test.mat')
end


disp('DETERMINE POINTS INTERCONNECTIVITIES')


NEW_METHOD =1;

if  NEW_METHOD == 1
%     
%     % STEP 1) COMPUTE DISTANCES BETWEEN CENTROIDS
%     
%     DISTANCES_CENTROIDS = zeros(size(CENT,2)) ;
%     for ipoints = 1:size(CENT,2)
%         d = CENT(:,ipoints) ;
%         DIST = bsxfun(@minus,CENT,d) ;
%         DISTANCES_CENTROIDS(:,ipoints)  = sqrt(sum(DIST.^2,1)) ;
%     end
    
    NEIGHBORING_CLUSTERS = cell(size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;
    
    DATAoffline = DefaultField(DATAoffline,'MaximumNumberNeighborsToSearchInterConnectivities',20) ;  
    
    nbuscar= length(NEIGHBORING_CLUSTERS)-1 ; 
    
    DATAoffline.MaximumNumberNeighborsToSearchInterConnectivities = min(DATAoffline.MaximumNumberNeighborsToSearchInterConnectivities,nbuscar) ; 
    
    for icluster = 1:length(NEIGHBORING_CLUSTERS)
        fprintf(1,'CLUST = %d\n',icluster)
        ipointsLOC = SNAPSHOTS_PER_CLUSTER_after_overlap{icluster} ;
        
        % Potential neighboring clusters  
        d = CENT(:,icluster) ;
        DIST_REST_CENTROIDS = bsxfun(@minus,CENT,d) ;
        DIST_REST_CENTROIDS =sqrt(sum(DIST_REST_CENTROIDS.^2,1)) ;
        [iii,NeighBoringClusters] = sort(DIST_REST_CENTROIDS,'ascend'); 
        NeighBoringClusters = NeighBoringClusters(2:DATAoffline.MaximumNumberNeighborsToSearchInterConnectivities+1) ; 
        
        for jclusterLOC = 1:length(NeighBoringClusters)
            jcluster = NeighBoringClusters(jclusterLOC) ; 
            if jcluster > icluster
            jpoints =  SNAPSHOTS_PER_CLUSTER_after_overlap{jcluster} ;
            INDint = intersect(ipointsLOC,jpoints) ;
            if ~isempty(INDint)
                NEIGHBORING_CLUSTERS{icluster}(end+1) =  jcluster ;
                NEIGHBORING_CLUSTERS{jcluster}(end+1) =  icluster ;
            end
            end
        end
        
    end
    
    
    
else
    
    
    NEIGHBORING_CLUSTERS = cell(size(SNAPSHOTS_PER_CLUSTER_after_overlap)) ;
    
    %     VECTORIAL =0 ;
    %     if  VECTORIAL == 1
    %         error('cellfun is not the solution !!!! ')
    %         for icluster = 1:length(NEIGHBORING_CLUSTERS)
    %             fprintf(1,'CLUST = %d\n',icluster)
    %             ipointsLOC = SNAPSHOTS_PER_CLUSTER_after_overlap{icluster} ;
    %
    %             jcluster = (icluster+1):length(NEIGHBORING_CLUSTERS) ;
    %             jpointsLOC =  SNAPSHOTS_PER_CLUSTER_after_overlap(jcluster) ;
    %
    %             ipointsLOC_cell = cell(size(jpointsLOC)) ;
    %             ipointsLOC_cell(:) = {ipointsLOC} ;
    %
    %             %INDint = intersect(ipointsLOC,jpointsLOC) ;
    %
    %             INDintCELL = cellfun(@intersect,ipointsLOC_cell,jpointsLOC,'UniformOutput',false) ;
    %
    %             INDint = cellfun(@isempty,INDintCELL) ;
    %
    %             IndicesNonEmpty_loc = find(INDint == 0) ;
    %             IndicesNonEmpty = jcluster(IndicesNonEmpty_loc) ;
    %
    %
    %             for jclusterLOC = 1:length(IndicesNonEmpty)
    %                 jclusterGLO = IndicesNonEmpty(jclusterLOC) ;
    %                 NEIGHBORING_CLUSTERS{icluster}(end+1) =  jclusterGLO ;
    %                 NEIGHBORING_CLUSTERS{jclusterGLO}(end+1) =  icluster ;
    %
    %             end
    %
    %
    %             %       %  if ~isempty(INDint)
    %             %       currentLENGTH = length(NEIGHBORING_CLUSTERS{icluster})+1  ;
    %             %       newLENGTH = currentLENGTH + length(IndicesNonEmpty)-1 ;
    %             %             NEIGHBORING_CLUSTERS{icluster}(currentLENGTH:newLENGTH) =  IndicesNonEmpty ;
    %
    %             %
    %             %             NEIGHBORING_CLUSTERS{IndicesNonEmpty}(end+1) =  icluster ;
    %             %       %  end
    %         end
    %
    %
    %
    %     else
    for icluster = 1:length(NEIGHBORING_CLUSTERS)
        fprintf(1,'CLUST = %d\n',icluster)
        ipointsLOC = SNAPSHOTS_PER_CLUSTER_after_overlap{icluster} ;
        
        for jcluster = (icluster+1):length(NEIGHBORING_CLUSTERS)
            jpointsLOC =  SNAPSHOTS_PER_CLUSTER_after_overlap{jcluster} ;
            INDint = intersect(ipointsLOC,jpointsLOC) ;
            if ~isempty(INDint)
                NEIGHBORING_CLUSTERS{icluster}(end+1) =  jcluster ;
                NEIGHBORING_CLUSTERS{jcluster}(end+1) =  icluster ;
            end
        end
        
    end
    
    % end
    
end