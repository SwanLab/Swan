function [indMINcl,indMINcl_LOC ]  = ...
    DistanceCheckEuclideanEFF(TransClustDATA,iCL,q,nq2)

if nargin == 0
    load('tmp1.mat')
end

qNEW = [] ;
indMINcl = iCL ;
nclusters = length(TransClustDATA.DistanceCentroid.C2) ;
DIST_cluster  = zeros(nclusters,1) ;
ErrorTransitionALL  = ones(nclusters,1) ; % Transition error

if ~isempty( TransClustDATA.NEIGHBORING_CLUSTERS)
    % Search is restricted to a set of possible neighbors
    
    
    
    NEIGS = [iCL, TransClustDATA.NEIGHBORING_CLUSTERS{iCL}] ; % Indexes neighboring clusters (included iCL itself)
    DIST_cluster = zeros(size(NEIGS)) ;
    for jclusterLOC = 1:length(NEIGS)
        
        jcluster = NEIGS(jclusterLOC) ;
        
      %  DIST_cluster(jclusterLOC) = TransClustDATA.DistanceCentroid.C2(jcluster)...
       %     -2*TransClustDATA.DistanceCentroid.Ct_BasisU{iCL,jcluster}*q +
       %     norm(q)^2 ;   % OLD !!!! 
        
        DIST_cluster(jclusterLOC) = TransClustDATA.DistanceCentroid.C2(jcluster)...
            -2*TransClustDATA.DistanceCentroid.Ct_BasisU{iCL,jclusterLOC}*q + norm(q)^2 ;
        
        
        % dL = BasisUall_cl{iCL}(DistanceCentroid.DOFlFE,DOFl)*VAR.DISP(DOFl)  ;
        %DIST_real(jcluster) = norm(DistanceCentroid.CENTROIDS_cl(DistanceCentroid.DOFlFE,jcluster)-dL)^2 ;
        %         if  jcluster~=iCL
        %             qNEW = TransClustDATA.TransMatrix{jcluster,iCL}*q ;
        %             ErrorTransitionALL(jcluster) = sqrt((nq2 - norm(qNEW)^2)/nq2) ;
        %         end
    end
    ErrorTransitionALL = [] ;
    [DIST_cluster_sorted,III] = sort(DIST_cluster) ;
    minDIST = DIST_cluster_sorted(1) ;
    indMINcl =  NEIGS(III(1));
    
      indMINcl_LOC = III(1)-1; % Local index; if zero, no cluster change 
    
else
    
    % CRITERION BASED ON MINIMUM EUCLIDEAN DISTANCE
    for jcluster = 1:nclusters
        DIST_cluster(jcluster) = TransClustDATA.DistanceCentroid.C2(jcluster)...
            -2*TransClustDATA.DistanceCentroid.Ct_BasisU{jcluster,iCL}*q + norm(q)^2 ;
        % dL = BasisUall_cl{iCL}(DistanceCentroid.DOFlFE,DOFl)*VAR.DISP(DOFl)  ;
        %DIST_real(jcluster) = norm(DistanceCentroid.CENTROIDS_cl(DistanceCentroid.DOFlFE,jcluster)-dL)^2 ;
        if  jcluster~=iCL
            qNEW = TransClustDATA.TransMatrix{jcluster,iCL}*q ;
            ErrorTransitionALL(jcluster) = sqrt((nq2 - norm(qNEW)^2)/nq2) ;
        end
    end
    [DIST_cluster_sorted,III] = sort(DIST_cluster) ;
    indMINcl  = III(1) ;
    
end