function [SNAPSHOTS_PER_CLUSTER_after_overlap,CENT] = ClosePointsClusteringStandard(SNAPdisp,DATAloc)

% Loop over number of points
TOLzero = 1e-14;
normSNAPdisp = sqrt(sum(SNAPdisp.^2,1)) ;
% There may be points which are zero. Such points are not considered
% centroids of any clusters
IndClusters = find(normSNAPdisp>TOLzero) ;  % IndClusters(j) returns the SNAPSHOT INDEX corresponding to CENTROID j
%IndSnapshotToCluster = zeros(npoints,1) ;
CENT = SNAPdisp(:,IndClusters);
SNAPSHOTS_PER_CLUSTER_after_overlap = cell(1,size(CENT,2)) ;
for ipointsLOC = 1:length(IndClusters)
    ipoints = IndClusters(ipointsLOC) ;
    d = SNAPdisp(:,ipoints) ;
    nd = norm(d) ;
    % Compute the distance between this point and the remaining ones
    DIST = bsxfun(@minus,SNAPdisp(:,IndClusters),d) ;
    nDIST  = sqrt(sum(DIST.^2,1)) ; % Euclidean distance
    [minDIST,indDIST] = sort(nDIST)  ;
    indDIST = setdiff(indDIST,ipointsLOC,'stable'); % We exclude the point itself
    
    % We seek to express d as a linear combination of its neighboring
    % points
    jpoints = 1;
    ERROR_approx = 1e20 ;
    while jpoints <=length(indDIST) && ERROR_approx > DATAloc.TolStopSearch
        indCANDIDATES = IndClusters(indDIST(1:jpoints)) ;
        BASIS =  orth(SNAPdisp(:,indCANDIDATES))  ;
        ERROR_approx = norm(d - BASIS*(BASIS'*d))/nd ;
        jpoints = jpoints + 1;
    end
    
%     if ipointsLOC ==688
%         warning('borrar esto')
%     end
%     
    
    if   jpoints == length(indDIST)+1
        disp(['Approximation error =',num2str(ERROR_approx)]) ;
        error(['Not converged for point =',num2str(ipoints)])
        
        
    end
    
    SNAPSHOTS_PER_CLUSTER_after_overlap{ipointsLOC} = [indCANDIDATES,ipoints]' ;
    
end