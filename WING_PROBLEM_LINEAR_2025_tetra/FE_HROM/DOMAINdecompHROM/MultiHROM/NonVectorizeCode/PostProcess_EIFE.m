function PostProcess_EIFE(COOR_coarse,CN_coarse, PROPMAT,MaterialType_COARSE,TypeElement_COARSE,DATA,TRANSF_COORD,...
    d_COARSE,stressGLO_coarse,Vrot,stressesREF);
% PLOT COARSE-SCALE/FINE-SCALE MESH/POST-ROCESS, EIFE METHOD
%
% JAHO, 12-march-2023
% ------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
    %DATA.DomainsToShowForCheckingMesh = [1,4,73:90];
end
%-------------------------------------------------------------------------------------------
% DATA = DefaultField(DATA,'ShowResultsGIDperDOMAIN',1) ;
NameFILE = [DATA.NameFileMesh,'EIFE'] ;

FOLDERprint = [cd,'/GIDPOST/'];
if ~exist(FOLDERprint)
    mkdir(FOLDERprint)
end
NameFileMesh = [FOLDERprint,NameFILE,'.msh'] ;
NameFileRes  = [FOLDERprint,NameFILE,'.res'] ;
%
DATA = DefaultField(DATA,'FineScaleDomainsToShow',1:length(MaterialType_COARSE)) ;
if isempty(DATA.FineScaleDomainsToShow)
    DATA.FineScaleDomainsToShow = 1:length(MaterialType_COARSE) ; 
end


nmesh = 1 + length(DATA.FineScaleDomainsToShow) ;
CN =  cell(nmesh,1) ;
TypeElement = cell(nmesh,1) ;
MaterialType = cell(nmesh,1) ;
NAMEMESH = cell(nmesh,1) ;
posgp =  cell(nmesh,1) ;
%VonMisesSTRESS =  cell(nmesh,1) ;
COOR = COOR_coarse ;
DISP = cell(nmesh) ; % DISPLACEMENTS

ILOC = 1;
MaterialType{ILOC} = MaterialType_COARSE ;
TypeElement{ILOC} = TypeElement_COARSE ;
nmatCOARSE= length(unique(MaterialType_COARSE)) ;
NAMEMESH{ILOC} ='COARSE';
CN{ILOC} = CN_coarse ;
posgp{ILOC} = []   ; % THIS SHOULD BE = 1
DISP{ILOC} =  d_COARSE ;

STRESSES_FINE  =cell(length(DATA.FineScaleDomainsToShow),1 ) ;


iacumNODES = size(COOR_coarse,1) ;

disp(['Recovering fine-scale variables domain ='])
DATA = DefaultField(DATA,'RecoveryStressesFromElasticityMatrix',1) ; 

for idom = 1:length(DATA.FineScaleDomainsToShow)
    
    e = DATA.FineScaleDomainsToShow(idom) ;
    disp(['e =',num2str(e)])
    
    TRANSF_COOR_loc = TRANSF_COORD{e};  % DATA STRUCTURE CONTAINING INFORMATION ABOUT THE DOMAIN
    IndexDomainLOC = TRANSF_COOR_loc.IndexParentDomain ;  % WHEN THERE ARE SEVERAL CN
    ROTATIONmat = TRANSF_COOR_loc.ROTATION ;
    
    if isfield(TRANSF_COOR_loc,'ROTATION_STRESSES')
    ROTATIONstresses = TRANSF_COOR_loc.ROTATION_STRESSES ;
    else
        ROTATIONstresses = [] ; 
    end
    
    
    
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
    dCOARSEloc = d_COARSE(DOFSloc) ;
    DISP{eFINE} = ReconstDisplacementsEIFE(dCOARSEloc,EIFE_prop,Vrot{e},DATA,SCALE_FACTOR_ALLS,ROTATIONmat,ndim) ;
    
    % RECOVERING FINE-SCALE STRESSES
    % -------------------------------------
    if ~isempty(ROTATIONstresses)
    STRESSES_FINE{idom} = ReconstStressesEIFE(dCOARSEloc,EIFE_prop,stressesREF{e},ROTATIONstresses,Vrot{e},DATA,SCALE_FACTOR_ALLS) ;
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
     
     STRESSES_FINE = {cell2mat(STRESSES_FINE')} ;   ; 
     
     
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
[NameFileRES ]= GidPostProcessEIFE(COOR,CN,TypeElement,cell2mat(DISP), STRESSES_FINE,...
    posgp,NameFileRes,MaterialType,DATAIN,DATA,stressGLO_coarse,IND_ELEM_MESHES);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
