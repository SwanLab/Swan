function [DATAoffline,iCL] = ClusterInitialSelection...
    (DATAoffline,VAR,DISP_CONDITIONS_cl,Fbody_cl,Ftrac_cl,istep,TransClustDATA,DATA_cl,OPERFE_cl,MATPRO_cl)


DATAoffline = DefaultField(DATAoffline,'INITIAL_CLUSTER',[]) ;  

if  ~isempty(DATAoffline.INITIAL_CLUSTER) && isempty(TransClustDATA.SequenceClusters)
    if  length(DATAoffline.INITIAL_CLUSTER) >1
        DominantModeCoordinate = zeros(size(DATAoffline.INITIAL_CLUSTER)) ;
        for   icluster =  1:length(DATAoffline.INITIAL_CLUSTER)
            iCL = DATAoffline.INITIAL_CLUSTER(icluster) ;
            disp('*********************************************************************+')
            disp(['Testing initial candidate cluster = ',num2str(iCL)])
            disp('*********************************************************************+')
            DISP_CONDITIONS = DISP_CONDITIONS_cl{iCL} ;
            VAR = BoundaryConditionsImpose(VAR,DISP_CONDITIONS_cl,Fbody_cl,Ftrac_cl,iCL,istep) ;
            
            [VAR,CONVERGED ]= NewtonRapshonStaticLarge(DATA_cl.COMMON,OPERFE_cl{iCL},VAR,MATPRO_cl{iCL},DISP_CONDITIONS_cl{iCL}.DOFl) ;
            DominantModeCoordinate(icluster) = abs(VAR.DISP(1)) ;
        end
        [~,IND_max] = max(DominantModeCoordinate) ;
        iCL = DATAoffline.INITIAL_CLUSTER(IND_max) ;
        disp('*********************************************************')
        disp(['Initial cluster = ',num2str(iCL)]) ;
        disp('*********************************************************')
        if ~isempty(DATAoffline.IdealIndexCluster)
            iCL = DATAoffline.IdealIndexCluster(istep) ;
        end
    else
        iCL = DATAoffline.INITIAL_CLUSTER ; 
    end
else
    iCL = TransClustDATA.SequenceClusters(1) ;
    
end