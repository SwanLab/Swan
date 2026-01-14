function PrintGID_1D_3D_representation(DATAINM,dFACES,DATAONLINE,pINT,alphaBC,BasisINT,COORref,NODESfaces,...
    BasisRdef,BasisRrb,rDEF,rRB,MAXstressVONMISES)

if nargin == 0
    load('tmp.mat')
end

% Interface modes, rigid body
DATAINM.ORTHOGONAL_RIGID_BODY_MODES= 0  ;

[BasisINTrb, xG]=  BasisRigidBodyInterface(COORref,NODESfaces,DATAINM)  ;

ndim  = size(COORref,2)
uBAR = DATAONLINE.uBAR ; 
dFACES_x = dFACES(1) ;
ndom = DATAONLINE.NdomX ;
COORx = linspace(xG(1),dFACES_x*ndom,ndom+1)' ;
COORy = xG(2)*ones(size(COORx)) ;
if ndim ==3 
COORz = xG(3)*ones(size(COORx)) ; 
COOR = [COORx,COORy,COORz] ; 
else
    COOR = [COORx,COORy] ; 
end
 

CN = [(1:ndom)',(2:ndom+1)'] ; 
if DATAINM.MODES_INTERFACE_RIGID_BODY ~=1
    error('Option not compatible')
else
    MODES_RIGID=  DATAINM.ModesINTERFACE_RigidB_select  ;  
    
    % Displacements 
    % ------------
    d = Displac1D_ROM(COORx,uBAR,BasisINTrb,BasisINT,MODES_RIGID,pINT,ndim) ; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Forces and moments 
    % ------------------
    %  
    f1 = NODESfaces{1} ;
    f2 = NODESfaces{3} ;
    f1 = small2large(f1,ndim) ; 
    f2 = small2large(f2,ndim) ; 
    Q_BasisRdef_f1 = BasisINTrb'*BasisRdef{1}(f1,:) ; 
    Q_BasisRdef_f2 = BasisINTrb'*BasisRdef{1}(f2,:) ; 
    Q_BasisRrb_f1 = BasisINTrb'*BasisRrb(f1,:) ; 
    Q_BasisRrb_f2 = BasisINTrb'*BasisRrb(f2,:) ; 
    
    nRB = size(BasisRrb,2) ; 
    nDEFr = size(BasisRdef{1},2) ; 
    GForces = zeros(nRB,ndom+1) ;
    
    ifinRB =0 ;
    ifinDEF = 0 ;
    for idom = 1:ndom
        iiniRB = ifinRB+1 ; ifinRB = iiniRB+nRB-1 ;
        iiniDEF = ifinDEF+1 ; ifinDEF = iiniDEF+nDEFr-1 ;
        Gf1 = Q_BasisRdef_f1*rDEF(iiniDEF:ifinDEF) + Q_BasisRrb_f1*rRB(iiniRB:ifinRB) ;
        Gf2 = Q_BasisRdef_f2*rDEF(iiniDEF:ifinDEF) + Q_BasisRrb_f2*rRB(iiniRB:ifinRB) ;
        if idom == 1
        GForces(:,idom) = -Gf1 ;
        end
        GForces(:,idom+1) = Gf2 ;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DATA = [] ; 
    DATA.Forces1D = [] ;
       DATA.Forces1D.Nx = GForces(1,:)' ;
     DATA.Forces1D.Ty = GForces(2,:)' ; 
     DATA.Forces1D.Mz = GForces(6,:)' ; 
    if ndim ==2  
     DATA.Forces1D.Tz = [] ;
     DATA.Forces1D.My = [] ; 
     DATA.Forces1D.Mx = [] ;
    else          
     DATA.Forces1D.Tz = GForces(3,:)' ; 
     DATA.Forces1D.My = GForces(5,:)' ;  
     DATA.Forces1D.Mx = GForces(4,:)' ;         
    end
    
    strainGLO =[] ; 
    stressGLO = [] ; 
    TypeElement = 'Linear' ;
    React = [] ; 
    posgp = [] ; 
    NAME_INPUT_DATA = 'Beam1D' ; 
    NameFileMesh = '' ; 
    MaterialType = []; 
   DATA.MAXstressVONMISES = MAXstressVONMISES ; 
    
    GidPostProcess(COOR,CN,TypeElement,d,strainGLO, stressGLO, ...
    React,NAME_INPUT_DATA,posgp,NameFileMesh,MaterialType,DATA);
    

end