function strainCOARSE_history = PostProcess_EIFE_vect(COOR_coarse,CN_coarse, PROPMAT,MaterialType_COARSE,TypeElement_COARSE,DATA,TRANSF_COORD,...
    d_COARSE,Vrot,stressesREF,VON_mises_COARSE,InternalVarStrain,OTHERinputs)
% PLOT COARSE-SCALE/FINE-SCALE MESH/POST-ROCESS, EIFE METHOD
% Stresses obtained by the vectorized code
%
% JAHO, 21-march-2023
% ------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
    %DATA.DomainsToShowForCheckingMesh = [1,4,73:90];
end
%-------------------------------------------------------------------------------------------
% DATA = DefaultField(DATA,'ShowResultsGIDperDOMAIN',1) ;
NameFILE = [DATA.NameFile_msh(1:end-4),'EIFE'] ;

% FOLDERprint = [cd,'/GIDPOST/'];
% if ~exist(FOLDERprint)
%     mkdir(FOLDERprint)
% end
NameFileMesh = [NameFILE,'.msh'] ;
NameFileRes  = [NameFILE,'.res'] ;

%OTHERinputs = DefaultField(OTHERinputs,'DATAinpGID',[1]) ; % Just For visualization purposes
%DATA.DATAinpGID = OTHERinputs.DATAinpGID; 

DATA = DefaultField(DATA,'FineScaleDomainsToShow',1:length(MaterialType_COARSE)) ;
if isempty(DATA.FineScaleDomainsToShow)
    DATA.FineScaleDomainsToShow = 1:length(MaterialType_COARSE) ;
else
    DATA.FineScaleDomainsToShow = intersect(1:length(MaterialType_COARSE),DATA.FineScaleDomainsToShow) ; 
    if isempty(DATA.FineScaleDomainsToShow)
        DATA.FineScaleDomainsToShow  = 1; 
    end
end

DATA = DefaultField(DATA,'UNIFORM_SCALING_REFERENCE_ELEMENT',1) ;


nmesh = 1 + length(DATA.FineScaleDomainsToShow) ;
CN =  cell(nmesh,1) ;
TypeElement = cell(nmesh,1) ;
MaterialType = cell(nmesh,1) ;
NAMEMESH = cell(nmesh,1) ;
posgp =  cell(nmesh,1) ;
%VonMisesSTRESS =  cell(nmesh,1) ;
COOR = COOR_coarse ;
DISP = cell(nmesh,1) ; % DISPLACEMENTS

ILOC = 1;
MaterialType{ILOC} = MaterialType_COARSE ;
TypeElement{ILOC} = TypeElement_COARSE ;
nmatCOARSE= length(unique(MaterialType_COARSE)) ;
NAMEMESH{ILOC} ='COARSE';
%    OTHERinputs.COORprintCOARSE = [] ; OTHERinputs.CNprintCOARSE = [] ; OTHERinputs.dPRINTcoarse = [] ; 

if ~isempty(OTHERinputs.CNprintCOARSE)
    % Case 27 Hexahedra elements
    CN{ILOC}  = OTHERinputs.CNprintCOARSE ;
    DISP{ILOC} =  OTHERinputs.dPRINTcoarse  ;
    COOR = OTHERinputs.COORprintCOARSE ;
else
    CN{ILOC} = CN_coarse ;
     DISP{ILOC} = d_COARSE  ;
     
end
posgp{ILOC} = []   ; % THIS SHOULD BE = 1

STRESSES_FINE  =cell(length(DATA.FineScaleDomainsToShow),1 ) ;
INTERNAL_VAR_FINE  =cell(length(DATA.FineScaleDomainsToShow),1 ) ;
VON_MISES_FINE  =cell(length(DATA.FineScaleDomainsToShow),1 ) ;


iacumNODES = size(COOR,1) ;


strainCOARSE_history = cell(size(CN{ILOC},1),1) ; 


disp(['Recovering fine-scale variables domain ='])

for idom = 1:length(DATA.FineScaleDomainsToShow)
    
    e = DATA.FineScaleDomainsToShow(idom) ;
    disp(['e =',num2str(e)])
    
    TRANSF_COOR_loc = TRANSF_COORD{e};  % DATA STRUCTURE CONTAINING INFORMATION ABOUT THE DOMAIN
    IndexDomainLOC = TRANSF_COOR_loc.IndexParentDomain ;  % WHEN THERE ARE SEVERAL CN
    ROTATIONmat = TRANSF_COOR_loc.ROTATION ;
    %    ROTATIONstresses = TRANSF_COOR_loc.ROTATION_STRESSES ;
    TRANSLATION = TRANSF_COOR_loc.TRANSLATION ;
    SCALE_FACTOR_ALLS = TRANSF_COOR_loc.SCALEFACTOR ;
    EIFE_prop = PROPMAT(MaterialType_COARSE(e)).EIFE_prop(IndexDomainLOC) ;
    MESHloc = EIFE_prop.MESH ; % Properties EIF object ("element")
    COORloc = MESHloc.COOR;
    
    C_ref = TRANSF_COOR_loc.CenterRotationREFERENCE ;  % Change 18-March-2023
    COORloc = bsxfun(@minus,COORloc',C_ref)' ;
    
    
    % TRANSFORMATION COORDINATES
    % --------------------------------------------------------------------------
    % COORloc --- Local coordinate fine-scale mesh  (centro = 0,0)
    % First: scaling
    if DATA.UNIFORM_SCALING_REFERENCE_ELEMENT  ==1
        COORloc = SCALE_FACTOR_ALLS*COORloc ;
        % ROTATION
        COORloc = (ROTATIONmat*COORloc')' ;
        % TRANSLATION
        ndim = size(COORloc,2) ;
        for idim = 1:ndim
            COORloc(:,idim) = COORloc(:,idim) + TRANSLATION(idim) ;
        end
    else
        error('Option not implement yet (to implement it, use the polynomial shape functions)')
    end
    % ---------------------------------------------------------------------------
    
    CNloc =  MESHloc.CN ;
    MaterialTypeLOC =MESHloc.MaterialType + nmatCOARSE;
    TypeElementLOC = MESHloc.TypeElement ;
    posgpLOC = MESHloc.posgp ;
    % Coordinate matrix
    COOR = [COOR; COORloc] ;
    % Connectivity matrix
    %  if DATAinputLOC.ShowResultsGIDperDOMAIN == 1
    eFINE = idom + 1;
    CN{eFINE} = CNloc + iacumNODES ;
    TypeElement{eFINE} = TypeElementLOC ;
    posgp{eFINE} = posgpLOC ;
    MaterialType{eFINE} = MaterialTypeLOC ;
    iacumNODES = size(COOR,1) ;
    if e <10
        NAMEMESH{eFINE} = ['SUBDOM_0',num2str(e)] ;
    else
        NAMEMESH{eFINE} = ['SUBDOM_',num2str(e)]  ;
    end
    
    
    % RECOVERING FINE-SCALE DISPLACEMENTS
    % -------------------------------------
    CNloc = CN_coarse(e,:) ;
    DOFSloc = small2large(CNloc,size(COORloc,2));
    dCOARSEloc = d_COARSE(DOFSloc,:) ;
    [DISP{eFINE},strainCOARSE_history{e}] = ReconstDisplacementsEIFE(dCOARSEloc,EIFE_prop,Vrot{e},DATA,SCALE_FACTOR_ALLS,ROTATIONmat,ndim) ;
    
    % RECOVERING FINE-SCALE STRESSES
    % -------------------------------------
    istrain = small2large(e,DATA.MESH.nstrain*DATA.MESH.ngaus_STRESS) ;
    stressesREF_e = stressesREF(istrain,:) ;
     
    
    if DATA.ReconstructionStressesUsingCoarseScaleDISP == 0
        [STRESSES_FINE{idom},VON_MISES_FINE{idom}] = ReconstStressesEIFE_viar(EIFE_prop,stressesREF_e,ROTATIONmat,DATA) ;
    else
        [STRESSES_FINE{idom},INTERNAL_VAR_FINE{idom},VON_MISES_FINE{idom}] = ReconstStressesEIFE_exact(EIFE_prop,dCOARSEloc,DATA,...
            SCALE_FACTOR_ALLS,PROPMAT(MaterialType_COARSE(e)).PROPMAT,Vrot{e},ROTATIONmat) ;
    end
    
    %
    
    
    %
    %     % Recovering stresses
    %     % -----------------------
    %     load(DATAinputLOC.WSstiffnessmatrices{idom},'Bdom','CgloDOM','Wdom') ;
    %     STRESS = CgloDOM*(Bdom*dROM{idom}) ;
    %     % It turns out CgloDOM already includes the weights
    %     nstrain = 6 ;
    %     for istrain =1:nstrain
    %         STRESS(istrain:nstrain:end) =  STRESS(istrain:nstrain:end)./Wdom ;
    %     end
    %     STRESS = reshape(STRESS,nstrain,[]) ;
    %     %   if DATAinputLOC.ShowResultsGIDperDOMAIN == 1
    %     VonMisesSTRESS{idom} =  VonMises_Stress(STRESS)' ;
    %  else
    %     VonMisesSTRESSLOC =  VonMises_Stress(STRESS) ;
    %     VonMisesSTRESS{1} = [ VonMisesSTRESS{1};VonMisesSTRESSLOC'  ] ;
    % end
    
    
end

DATA.ndim  = ndim ;


DATA= DefaultField(DATA,'SeparatedFineScaleDomainsPOST',0);

if DATA.SeparatedFineScaleDomainsPOST == 0
    CNall = CN ;
    MaterialTypeALL = MaterialType;
    TypeElementALL = TypeElement;
    NAMEMESHall = NAMEMESH;
    posgpALL = posgp ;
    CN = cell(2,1) ;
    %  elemINI = size(CNall{1},1) ;
    CN{1} = CNall{1} ;
    CN{2} = cell2mat(CNall(2:end)) ;
    
    MaterialType = cell(2,1) ;
    MaterialType{1} = MaterialTypeALL{1};
    nmat = length(unique(MaterialType{1})) ;
    MaterialType{2} = cell2mat(MaterialTypeALL(2:end)) + nmat ;
    
    TypeElement = cell(2,1) ;
    TypeElement{1} = TypeElementALL{1} ;
    TypeElement{2} = TypeElementALL{2} ;
    
    
    posgp  = cell(2,1) ;
    posgp{1} = posgpALL{1} ;
    posgp{2} = posgpALL{2} ;
    
    NAMEMESH   =cell(2,1) ;
    NAMEMESH{1} = NAMEMESHall{1} ;
    NAMEMESH{2} = 'FINEmesh' ;
    
    STRESSES_FINE = {cell2mat(STRESSES_FINE)} ;   ;
    VON_MISES_FINE = {cell2mat(VON_MISES_FINE)} ;   ;
      INTERNAL_VAR_FINE = {cell2mat(INTERNAL_VAR_FINE)} ;   ;
end



% Writing mesh
% -----------
NAME_INPUT_DATA = 'MeshEIFE' ;
IND_ELEM_MESHES = GidMesh2DFE_multiEIFE(NameFileMesh,COOR,CN,NAME_INPUT_DATA,MaterialType,...
    TypeElement,NAMEMESH);

% strainGLO = [] ;

% React = [] ;
DATAIN = [] ;
DATA.NAMEMESH = NAMEMESH;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NameFileRES ]= GidPostProcessEIFE_vect(COOR,CN,TypeElement,cell2mat(DISP), STRESSES_FINE,...
    posgp,NameFileRes,MaterialType,DATAIN,DATA,VON_mises_COARSE,IND_ELEM_MESHES,InternalVarStrain,INTERNAL_VAR_FINE,VON_MISES_FINE);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
