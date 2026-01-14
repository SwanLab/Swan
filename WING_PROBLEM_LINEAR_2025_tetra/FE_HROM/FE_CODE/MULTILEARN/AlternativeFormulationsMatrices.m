function DATAwmethod = AlternativeFormulationsMatrices(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
    THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb,BasisUdef,BasisUrb,DATAwmethod,DATAINM,COORref,...
    INDrigR,INDdefR,INDrig,INDdef,uBAR)



if DATAINM.MinimizationBoundaryWork >0
    [ DATAwmethod.Ework]=  AssemblyEwork(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
        THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb,BasisUdef,BasisUrb);
    
end

if DATAINM.ConstrainedWithReactionsAtInterfaces == 1
    % Checking whether it is a 1D -tiling problem
    FACES_CONTACT = sum(THETAfaces,1) ;
    IND_CONTACT= find(FACES_CONTACT) ;
    if sum(IND_CONTACT-[1 3]) ~=0
        error('Option not compatible. ONly 1D-tiling problems are allowed (along x-axis)')
    end
    
    [J]=  AssemblyQint(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
        THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb,BasisUdef,BasisUrb,COORref);
    DATAwmethod.JintRB = J(:,INDrigR) ;
    DATAwmethod.JintDEF = J(:,INDdefR) ;
    Jr_r = J(:,INDrigR)*rRBglo ;
    DATAwmethod.Jr_r =  Jr_r  ;
end

DATAINM =DefaultField(DATAINM,'SimpleProjection',0);
if DATAINM.SimpleProjection == 1
    % Checking whether it is a 1D -tiling problem
    FACES_CONTACT = sum(THETAfaces,1) ;
    IND_CONTACT= find(FACES_CONTACT) ;
    if sum(IND_CONTACT-[1 3]) ~=0
        error('Option not compatible. ONly 1D-tiling problems are allowed (along x-axis)')
    end
    
    [P,b]=  AssemblyPsimple(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
        THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb,BasisUdef,BasisUrb,COORref,uBAR);
    DATAwmethod.Pr = P(:,INDrig) ;
    DATAwmethod.Pd = P(:,INDdef) ;
    DATAwmethod.b  = b  ;
    
     [T]=  AssemblyTsimple(BasisRdef,NODESbound,reactDOMrbGLO,betaBC,COLUMNS_RVE,COLUMNS_RVEloc,...
        THETAfaces,GAMMAfaces,ndim,alphaBC,rRBglo,BasisRrb,BasisUdef,BasisUrb,COORref,uBAR);
     DATAwmethod.Td = T(:,INDdefR) ;
    DATAwmethod.cTr = T(:,INDrigR)*rRBglo ;    
   
end