function NODESV_np1 = UpdateKinematicVariables(OPERfe,PRES_DISP,BND_Dyn,NODESV_np1,istep)

if nargin == 0
    load('tmp1.mat')
end

gBOUND = PRES_DISP.VALUE*PRES_DISP.TIME_FACTOR(istep) ; %  BBnw_s*gBOUND(:,istep)  ;

if  ~isempty(OPERfe.DOFm)
    NODESV_np1.U(OPERfe.DOFs) = gBOUND + OPERfe.Gbound*NODESV_np1.U(OPERfe.DOFm) ;
    if  ~isempty(BND_Dyn.gBOUNDd)
        NODESV_np1.VEL(OPERfe.DOFs) = gBOUNDd(:,istep) + OPERfe.Gbound*NODESV_np1.VEL(OPERfe.DOFm) ;
        NODESV_np1.ACEL(OPERfe.DOFs) = gBOUNDd(:,istep) + OPERfe.Gbound*NODESV_np1.ACEL(OPERfe.DOFm) ;
        
    end
else
    NODESV_np1.U(OPERfe.DOFs) = gBOUND ;
    if  ~isempty(BND_Dyn.gBOUNDd)
        NODESV_np1.VEL(OPERfe.DOFs) = gBOUNDd(:,istep)   ;
        NODESV_np1.ACEL(OPERfe.DOFs) = gBOUNDd(:,istep)   ;
        
    end
    
    
end