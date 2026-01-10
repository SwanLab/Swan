function strainCOARSE_history = PostProcess_EIFE_vectBUBcorot(COOR_coarse,CN_coarse, PROPMAT,MaterialType_COARSE,TypeElement_COARSE,DATA,TRANSF_COORD,...
    d_COARSE,stressesREF,VON_mises_COARSE,InternalVarStrain,OTHERinputs,...
    QrotTIME,QrotINI,LboolCall)
%PostProcess_EIFE_vectBUBcorot Co-rotational postprocessing for EIFEM simulations with bubble enrichment
%
%   This function performs advanced postprocessing of results from the Empirical
%   Interscale Finite Element Method (EIFEM), incorporating bubble modes and
%   co-rotational corrections. It reconstructs fine-scale meshes and variables
%   (displacements, stresses, internal strains) based on coarse-scale data and
%   geometric transformation maps. Results are written to GiD-compatible files.
%
%   INPUT:
%     COOR_coarse       - Nodal coordinates of the coarse-scale mesh
%     CN_coarse         - Element connectivity of the coarse-scale mesh
%     PROPMAT           - Material properties, including EIFE domains
%     MaterialType_COARSE - Material type for each coarse element
%     TypeElement_COARSE - Element type for the coarse-scale mesh
%     DATA              - Structure with postprocessing configuration and filenames
%     TRANSF_COORD      - Geometric maps (scaling, rotation, centering) per domain
%     d_COARSE          - Coarse-scale displacements (with bubbles)
%     stressesREF       - (Optional) Reference stress field (e.g., PK2)
%     VON_mises_COARSE  - (Optional) Precomputed von Mises values for coarse elements
%     InternalVarStrain - (Optional) Internal variable values per element (e.g. strain energy)
%     OTHERinputs       - Struct with auxiliary data (interpolated displacements, node mappings, etc.)
%     QrotTIME          - Rotation matrices per time step for each domain (co-rotational)
%     QrotINI           - Initial rotation matrices per domain
%     LboolCall         - Logical matrix selecting which DOFs are active (e.g., after SVD)
%
%   OUTPUT:
%     strainCOARSE_history - Cell array with time histories of strain fields (coarse elements)
%
%   FUNCTIONALITY:
%   ------------------------------------------------------------------------------
%   1. Builds GiD-compatible mesh structure, including fine-scale subdomains
%   2. Transforms initial reference coordinates to global frame using scaling and rotations
%   3. Recovers fine-scale displacements using `ReconstDisplacementsEIFEbubCOROT`
%   4. Reconstructs stress and strain fields using `ReconstStressesEIFE_exactBUBcorot`
%   5. Optionally merges fine-scale subdomains into a single mesh for visualization
%   6. Calls GiD postprocessing routines to output *.msh and *.res files
%
%   ASSUMPTIONS AND SPECIAL CASES:
%   - Coarse elements with 27 → 26 nodes (hexes) can be corrected via `OTHERinputs`
%   - Supports anisotropic scaling and arbitrary local orientations (via Qrot)
%   - Allows per-domain postprocessing or aggregation depending on DATA.SeparatedFineScaleDomainsPOST
%
%   SEE ALSO:
%     ReconstDisplacementsEIFEbubCOROT, ReconstStressesEIFE_exactBUBcorot,
%     GidPostProcessEIFE_vect, GidMesh2DFE_multiEIFE
%
%   AUTHOR:
%     Joaquín A. Hernández Ortega (JAHO), UPC BarcelonaTech, 6-Nov-2024
%     Adapted from PostProcess_EIFE_vectBUB.m to include large rotations and co-rotational corrections.
%     COmments by ChatGPT-4
 % /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROT.tex
% ------------------------------ ------------------------------------------------------------
% JAHO, 6-Nov-2024, Balmes 185, Barcelona.
% ------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
%-------------------------------------------------------------------------------------------
% DATA = DefaultField(DATA,'ShowResultsGIDperDOMAIN',1) ;
NameFILE = [DATA.NameFile_msh(1:end-4),'EIFE'] ;
MESHextended = DATA.MESHextended ;

% FOLDERprint = [cd,'/GIDPOST/'];
% if ~exist(FOLDERprint)
%     mkdir(FOLDERprint)
% end
NameFileMesh = [NameFILE,'.msh'] ;
NameFileRes  = [NameFILE,'.res'] ;
%
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
OTHERinputs = DefaultField(OTHERinputs,'MaterialType_coarse_GID',MaterialType_COARSE) ; 

MaterialType{ILOC} = OTHERinputs.MaterialType_coarse_GID ;  

TypeElement{ILOC} = TypeElement_COARSE ;
nmatCOARSE= length(unique(MaterialType_COARSE)) ;
NAMEMESH{ILOC} ='COARSE';
%    OTHERinputs.COORprintCOARSE = [] ; OTHERinputs.CNprintCOARSE = [] ; OTHERinputs.dPRINTcoarse = [] ;

IndicesBoundaryNODES = DATA.MESHextended.DOFS_bLOC ;
ndimFINE = size(PROPMAT(MaterialType_COARSE(1)).EIFE_prop(1).MESH.COOR,2) ;
ndimCOAR =  size(COOR_coarse,2) ;

if ~isempty(OTHERinputs.CNprintCOARSE)
    % Case 27 Hexahedra elements
    CN{ILOC}  = OTHERinputs.CNprintCOARSE ;
    DISP{ILOC} =  OTHERinputs.dPRINTcoarse  ;
    COOR = OTHERinputs.COORprintCOARSE ;
else
    CN{ILOC} = CN_coarse ;
    
    if  ndimFINE ~= size(COOR_coarse,2)
        
        
        IndicesBoundaryNODES_all = small2large(IndicesBoundaryNODES,ndimFINE) ;
        DISP{ILOC}  = zeros(length(IndicesBoundaryNODES_all),size(d_COARSE,2)) ;
        DISP{ILOC}(1:ndimFINE:end) =  d_COARSE(IndicesBoundaryNODES,:)  ;
        
    else
        DISP{ILOC} = d_COARSE(IndicesBoundaryNODES,:)  ;
        
    end
    
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
    
    %% INITIAL COORDINATES, FINE SCALE
    % m The\textbf{ initial fine-scale coordinates} of element $\domCe{e}$ are given by
    %  \XfNe{I}{e}  =  \CentCgloE{e} + \QrotINIe{e}  \lambdaLENe{e} (\XfREFnE{I}{e} -\CentClocREFe{e}), \hspace{0.5cm}  I =1,2 \ldots \nnodeF
    eIND = ((e-1)*ndimFINE+1):e*ndimFINE ;
    QrotINIe = QrotINI(eIND,:) ;
    
    lambdaLENe = TRANSF_COORD{e}.lambdaLEN ;
    CentCgloE = TRANSF_COORD{e}.CentCglo ;
    CentClocREFe = TRANSF_COORD{e}.CentClocREF ;
    IndexDomainLOC = TRANSF_COORD{e}.IndexParentDomain ;  % Family of EIF elements to which Omega_e belongs
    EIFE_prop = PROPMAT(MaterialType_COARSE(e)).EIFE_prop(IndexDomainLOC) ;
    
    
    %     TRANSF_COOR_loc = TRANSF_COORD{e};  % DATA STRUCTURE CONTAINING INFORMATION ABOUT THE DOMAIN
    %     ROTATIONmat = TRANSF_COOR_loc.ROTATION ;
    %     %    ROTATIONstresses = TRANSF_COOR_loc.ROTATION_STRESSES ;
    %     TRANSLATION = TRANSF_COOR_loc.TRANSLATION ;
    %     SCALE_FACTOR_ALLS = TRANSF_COOR_loc.SCALEFACTOR ;
    MESHloc = EIFE_prop.MESH ; % Properties EIF object ("element")
    XfREF = MESHloc.COOR;  % Mesh parent domain
    
    %  \lambdaLENe{e} (\XfREFnE{I}{e} -\CentClocREFe{e})
    Xf   =   lambdaLENe*bsxfun(@minus,XfREF',CentClocREFe)' ;
    %  \QrotINIe{e}  \lambdaLENe{e} (\XfREFnE{I}{e} -\CentClocREFe{e})
    Xf = (QrotINIe*Xf')' ;
    % \XfNe{I}{e}  =  \CentCgloE{e} + \QrotINIe{e}  \lambdaLENe{e} (\XfREFnE{I}{e} -\CentClocREFe{e})
    Xf   =    bsxfun(@plus,Xf',CentCgloE)' ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CNloc =  MESHloc.CN ;
    MaterialTypeLOC =MESHloc.MaterialType + nmatCOARSE;
    TypeElementLOC = MESHloc.TypeElement ;
    posgpLOC = MESHloc.posgp ;
    % Coordinate matrix
    if size(COOR,2) ~= size(Xf,2)
        % Case fine/coarse meshes different dimensions
        COORamp = zeros(size(COOR,1),size(Xf,2)) ;
        COORamp(:,1:size(COOR,2))  = COOR;
        COOR = [COORamp; Xf] ;
    else
        COOR = [COOR; Xf] ;
    end
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
    %     CNloc = MESHextended.CN(e,:) ; %  CN_coarse(e,:) ;
    %     DOFSloc_all = small2large(CNloc,MESHextended.NDOFS_pernode);
    %     [AAA,INDLOC] = setdiff(DOFSloc_all,MESHextended.DOFS_ghost,'stable') ;
    %     [DOFSloc,INDLOCa,INDLOC ]= intersect(AAA,MESHextended.DOFS_TO_KEEP,'stable') ;
    
    dCe = LboolCall(DATA.MESH.IndexDOFS_per_element{e},:)*d_COARSE;
    
    [DISP{eFINE},strainCOARSE_history{e},dClocINCREe_time] =...
        ReconstDisplacementsEIFEbubCOROT(XfREF,dCe,EIFE_prop,QrotINIe,DATA,lambdaLENe,QrotTIME,size(COOR,2),eIND) ;
    
    % RECOVERING FINE-SCALE STRESSES
    % -------------------------------------
    %istrain = small2large(e,DATA.MESH.nstrain*DATA.MESH.ngaus_STRESS) ;
    
    %if ~isempty(stressesREF)
    %    stressesREF_e = stressesREF(istrain,:) ;
    %end
    
    
    %if DATA.ReconstructionStressesUsingCoarseScaleDISP == 0
    %    [STRESSES_FINE{idom},VON_MISES_FINE{idom}] = ReconstStressesEIFE_viar(EIFE_prop,stressesREF_e,ROTATIONmat,DATA) ;
    %else
    [STRESSES_FINE{idom},INTERNAL_VAR_FINE{idom},VON_MISES_FINE{idom}] = ...
        ReconstStressesEIFE_exactBUBcorot(EIFE_prop,dClocINCREe_time,DATA,...
        lambdaLENe,PROPMAT(MaterialType_COARSE(e)).PROPMAT) ;
    %end
    
    %BUB
    
    
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

DATA.ndim  = DATA.MESH.ndim ;


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
DATA = DefaultField(DATA,'OPEN_GID_POST',0) ; 
DATAIN.OPEN_GID_POST = DATA.OPEN_GID_POST ; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NameFileRES ]= GidPostProcessEIFE_vect(COOR,CN,TypeElement,cell2mat(DISP), STRESSES_FINE,...
    posgp,NameFileRes,MaterialType,DATAIN,DATA,VON_mises_COARSE,IND_ELEM_MESHES,InternalVarStrain,INTERNAL_VAR_FINE,VON_MISES_FINE);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
