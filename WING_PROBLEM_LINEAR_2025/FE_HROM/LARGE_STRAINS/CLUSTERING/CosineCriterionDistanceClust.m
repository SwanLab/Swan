function [indMINcl] = CosineCriterionDistanceClust(nq2,TransClustDATA,iCL,q)


nclusters = length(TransClustDATA.DistanceCentroid.C2) ;
DIST_cluster  = zeros(nclusters,1) ;
ErrorTransitionALL  = ones(nclusters,1) ; % Transition error
indMINcl = iCL ;
qNEW1_all = zeros(size(ErrorTransitionALL));
%   CRITERION  BASED ON "COSINE"
if nq2 >1e-16
    for jcluster = 1:nclusters
        NOMIN =   TransClustDATA.DistanceCentroid.Ct_BasisU{jcluster,iCL}*q  ;
        DENOMIN =   sqrt(nq2*TransClustDATA.DistanceCentroid.C2(jcluster)) ;
        DIST_cluster(jcluster) =  1- (NOMIN/DENOMIN) ;
        if  jcluster~=iCL
            qNEW = TransClustDATA.TransMatrix{jcluster,iCL}*q ;
            %disp(['qNEW(1)_',num2str(jcluster),'=',num2str(qNEW(1)),' qOLD(1) = ',num2str(q(1))]) ;
            ErrorTransitionALL(jcluster) = sqrt((nq2 - norm(qNEW)^2)/nq2) ;
            qNEW1_all(jcluster) = qNEW(1) ;
            % if DATAoffline.EvaluateResidualOtherClusters == 1
            %end
        else
            qNEW1_all(jcluster) = q(1) ;
        end
    end
    [sortDISTANCES,III] = sort(DIST_cluster) ;
    indMINcl  = III(1) ;
    [sortqNEW,III2] = sort(qNEW1_all) ;
    indMINcl_q1  = III2(1) ;
    [minTRANS,indMINcl_transition] = min(ErrorTransitionALL) ;
    disp(['Minimum transition error = ',num2str(minTRANS),' CLUSTER =',num2str(indMINcl_transition)]) ;
    disp(['Maximum q1 = ',num2str(sortqNEW(1)),' CLUSTER =',num2str(indMINcl_q1)]) ;
else
    disp(['...Zero displacements'])
end



if ~isempty(DATAoffline.IdealIndexCluster)  && ~isempty(qNEW1_all)
    % This is just for testing purposes
    if  DATAoffline.IdealIndexCluster==1
        indMINcl = DATAoffline.IdealIndexCluster(istep) ;
    else
        indIDEAL = DATAoffline.IdealIndexCluster(istep) ;
        %      if indIDEAL ~= indMINcl
        disp(['Current cluster = ',num2str(iCL),' DIST centr =',num2str(DIST_cluster(iCL)),' q1 =',num2str(qNEW1_all(iCL))])  ;
        disp(['Next cluster = ',num2str(indMINcl),' DIST centr =',num2str(DIST_cluster(indMINcl)),' q1 =',num2str(qNEW1_all(indMINcl))])  ;
        disp(['Ideal cluster =',num2str(indIDEAL),' DIST centr =',num2str(DIST_cluster(indIDEAL)),' q1 =',num2str(qNEW1_all(indIDEAL))])  ;
        %     end
    end
end