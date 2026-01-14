function CheckFineScaleMeshGID(COOR_coarse,CN_coarse, PROPMAT,MaterialType_COARSE,TypeElement_COARSE,DATA,TRANSF_COORD)
% PLOT FINE-SCALE MESH, EIFE METHOD
% NO RESULTS, JUST THE MESH ---FOR CHECKING PURPOSES
% JAHO, 11-march-2023
% ------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
    %DATA.DomainsToShowForCheckingMesh = [1,4,73:90]; 
end
%-------------------------------------------------------------------------------------------
% DATA = DefaultField(DATA,'ShowResultsGIDperDOMAIN',1) ;
NameFILE = [DATA.NameFileMesh,'_check'] ; 

FOLDERprint = [cd,'/GIDPOST/'];
if ~exist(FOLDERprint)
    mkdir(FOLDERprint)
end
NameFileMesh = [FOLDERprint,NameFILE,'.msh'] ;
NameFileRes  = [FOLDERprint,NameFILE,'.res'] ;
%
DATA = DefaultField(DATA,'DomainsToShowForCheckingMesh',1:length(MaterialType_COARSE)) ; 


nmesh = 1 + length(DATA.DomainsToShowForCheckingMesh) ; 
CN =  cell(nmesh,1) ;
TypeElement = cell(nmesh,1) ;
MaterialType = cell(nmesh,1) ;
NAMEMESH = cell(nmesh,1) ;
posgp =  cell(nmesh,1) ;
%VonMisesSTRESS =  cell(nmesh,1) ;
COOR = COOR_coarse ; 

ILOC = 1; 
MaterialType{ILOC} = MaterialType_COARSE ; 
TypeElement{ILOC} = TypeElement_COARSE ; 
nmatCOARSE= length(unique(MaterialType_COARSE)) ; 
NAMEMESH{ILOC} ='COARSE'; 
CN{ILOC} = CN_coarse ; 
posgp{ILOC} = []   ; % THIS SHOULD BE = 1 
 
iacumNODES = size(COOR_coarse,1) ;
     
    
for idom = 1:length(DATA.DomainsToShowForCheckingMesh)     
    
    e = DATA.DomainsToShowForCheckingMesh(idom) ; 
    
     TRANSF_COOR_loc = TRANSF_COORD{e};  % DATA STRUCTURE CONTAINING INFORMATION ABOUT THE DOMAIN
     IndexDomainLOC = TRANSF_COOR_loc.IndexParentDomain ;  % WHEN THERE ARE SEVERAL CN
     ROTATIONmat = TRANSF_COOR_loc.ROTATION ; 
     TRANSLATION = TRANSF_COOR_loc.TRANSLATION ; 
     SCALE_FACTOR_ALLS = TRANSF_COOR_loc.SCALEFACTOR ; 
     MESHloc = PROPMAT(MaterialType_COARSE(e)).EIFE_prop(IndexDomainLOC).MESH ; % Properties EIF object ("element")
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
     for idim = 1:size(COORloc,2)
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

% if DATA.ShowResultsGIDperDOMAIN == 0
%     CN = {cell2mat(CN)} ;
%     MaterialType = {cell2mat(MaterialType)} ;
%     TypeElement = {TypeElement{1}} ; 
%     NAMEMESH = {'ALL'} ; 
%     posgp = {posgp{1}}; 
%   %  VonMisesSTRESS = {cell2mat(VonMisesSTRESS)} ; 
% end


% Writing mesh
% -----------
NAME_INPUT_DATA = 'MeshROM' ;
IND_ELEM_MESHES = GidMesh2DFE_multiEIFE(NameFileMesh,COOR,CN,NAME_INPUT_DATA,MaterialType,...
    TypeElement,NAMEMESH);

disp(['Open GID POST FILE'])
disp(NameFileMesh)

% % -------------------
% % WRITING .RES FILE
% % -------------------
% %STRESS3D = STRESS3D(:) ;
% strainGLO = [] ;
% React = [] ;
% RESIDUAL = [] ;
% d = cell2mat(dROM) ;
% 
% strainGLO = [] ;
% STRESS3D = [] ;
% React = [] ;
% NODES = [] ;
% DATAIN = [] ;
% DATA.NAMEMESH = NAMEMESH;
% STRESSDATA.VONMISES =VonMisesSTRESS ;
% FORCES_2_PRINT = [] ;
% DATA_REFMESH = [] ; DATARUN = [] ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [NameFileRES ]= GidPostProcessMULTI_coar(COOR,CN,TypeElement,d,strainGLO, STRESS3D, ...
%     React,posgp,NameFileRes,MaterialType,NODES,RESIDUAL,DATAIN,DATA,STRESSDATA,IND_ELEM_MESHES,...
%     FORCES_2_PRINT,DATA_REFMESH,DATARUN);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
