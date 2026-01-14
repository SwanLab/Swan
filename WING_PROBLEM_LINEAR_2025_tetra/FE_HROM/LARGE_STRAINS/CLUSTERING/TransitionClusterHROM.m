function [iCLnew,   VAR,ErrorTransition ]= TransitionClusterHROM(iCL,TransClustDATA,VAR,DATAoffline,DATA_cl,...
    OPERFE_cl,MATPRO_cl,DISP_CONDITIONS_cl,istep,Fbody_cl,Ftrac_cl)
if nargin == 0
    load('tmp1.mat')
end
% Check whether transition between clusters takes place
q = VAR.DISP(DISP_CONDITIONS_cl{iCL}.DOFl) ;  % Current reduced coordinates
nq2  = norm(q)^2 ;  % Norm current reduced coordinates
ErrorTransition = 0 ;  % Error associated to the transition
iCLnew = iCL ;
qNEW = [] ; 
indMINcl = iCL ; 

if   ~isempty(TransClustDATA.SequenceClusters)
    % The sequence of clusters have been computed a priori
    % TransClustDATA.SequenceClusters IS COMPUTED A PRIORI
    if istep < length(TransClustDATA.SequenceClusters)
        indMINcl = TransClustDATA.SequenceClusters(istep + 1) ;
    end
else
    % Here we have to check whether it is necessary to switch from cluster
    % iCL to any other cluster
    
    % This was an attempt  to develop a criterion based on the evaluation
    % of the residual (abandoned )
    DATAoffline = DefaultField(DATAoffline,'EvaluateResidualOtherClusters',0) ;
    if DATAoffline.EvaluateResidualOtherClusters == 1
        ErrorResidualOtherClusters = ResidualErrorOtherClusters(iCL,TransClustDATA,VAR,DATAoffline,DATA_cl,...
            OPERFE_cl,MATPRO_cl,DISP_CONDITIONS_cl,istep,Fbody_cl,Ftrac_cl) ;
    end
    
    if DATAoffline.CriterionSwitchCluster == 0
        % Euclidean distance
        [indMINcl]  =  DistanceCheckEuclidean(TransClustDATA,iCL,q,nq2) ;
    else
        % Cosine distance  (euclidean distance with normalized snapshots)
        [indMINcl]= ...
            CosineCriterionDistanceClust(nq2,TransClustDATA,iCL,q) ;
    end
    
end


%%%  TRANSITION BETWEEN CLUSTERS TAKE PLA

if indMINcl ~= iCL
    
    iCLnew = indMINcl ;
    disp(['Transition from cluster ',num2str(iCL),' to cluster ',num2str(indMINcl)]) ;
    
    qNEW = TransClustDATA.TransMatrix{iCLnew,iCL}*q ;
    
    ErrorTransition  = sqrt(abs(nq2 - norm(qNEW)^2)/nq2) ;
    
    
    
    disp(['Error transition =',num2str(ErrorTransition)])
    
    % Initialization
    VAR.DISP = zeros(size(VAR.DISP )) ;
    VAR.DISP(1:length(qNEW)) = qNEW ;  %
    
end
