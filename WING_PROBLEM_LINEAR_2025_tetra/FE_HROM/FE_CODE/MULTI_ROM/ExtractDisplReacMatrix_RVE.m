function [dDOM,reactDOM,stressDOM,DATA_REFMESH,nDOMx,nDOMy,timestepINCLUDE,REMOVE_FIRST_DOMAIN ]= ...
    ExtractDisplReacMatrix_RVE(INPUT_DATAFILE,DATAIN,iproject)

% Output: Matrix of displacements (dDOM), and matrix of reactions (reactDOM)
% DATAOUT: Structure containing mesh information of each domain
if nargin == 0
    load('tmp.mat')
end
NameFileMeshLOC_coarse = [] ;
% -------------------------
% (EXTRACTION OF PROPERTIES FROM FE FILES)
% -------------------------
load(INPUT_DATAFILE,'DATA_INPUT_FE') ;  % DATA_INPUT_FE --> Input data file
if exist(INPUT_DATAFILE,'file') %
    LIST_VAR = who('-file',INPUT_DATAFILE) ;
    if ~ismember(LIST_VAR,'Bst' ) & ~isempty(DATA_INPUT_FE.nameWORKSPACE_Kstiff)
        % This means that the Bst matrix is not stored in this mat
        % file.
        % Then the required information is in
        load(DATA_INPUT_FE.nameWORKSPACE_Kstiff,'MaterialType','COOR','CN','TypeElement','posgp',...
            'CNb','TypeElementB','DOMAINVAR','CONNECTb','IndicesRenumberingElements',...
            'Nst',...
            'wSTs_RHS','posgp_RHS','Fpnt','Tnod','CNb','Cglo','wSTs','Bst','NameFileMesh','density','AREA',...
            'NameFileMeshLOC_coarse') ;
        load(INPUT_DATAFILE,'d','fNOD','Fpnt','Tnod','CNb','stressGLO') ; % 11-Jun-09, Nonlinear
        % Displecement field (and project-specific variablres) are
        % retrieved from the actual mat-file ssociated to the project the project
    else
        
        load(INPUT_DATAFILE,'MaterialType','COOR','CN','d','TypeElement','posgp',...
            'CNb','TypeElementB','DOMAINVAR','CONNECTb','IndicesRenumberingElements',...
            'wSTs_RHS','fNOD','posgp_RHS','Fpnt','Tnod','wSTs','NameFileMesh',...
            'AREA','density','NameFileMeshLOC_coarse') ;
        
        if ~ismember(LIST_VAR,'Bst' )
            % Bst, Nst, and Cglo were not stored in memory ---
            % (because  DATA.DO_NOT_STORE_STIFFNESS_AND_B_MATRIX  has been
            % set to data in the pertinent FE analyis
            Bst = [] ; Nst = [] ; Cglo = [] ;
            
        else
            load(INPUT_DATAFILE,'Bst','Nst','Cglo','stressGLO') ; % Nonlinear, 11-Jun-09
        end
    end
else
    error(['Non existing FE data !!! '])
end


%DATAIN = DefaultField(DATAIN,'ORTHOGONALITY_RIGIDB_MASS_MATRIX',1) ;


%%%%%%%%%%
% IndicesRenumberingElements is created when changing the numbering of
%  % elements
% ELASTOSTATIC_GEN/FE_CODE/SolveElastFE.m
% IndicesRenumberingElements
[dummy IndicesRenumberingElements_INV] = sort(IndicesRenumberingElements) ;
DOMAINVAR.ListElements = DOMAINVAR.ListElements(:) ; % We convert the ndomX x ndomY cell array into a ndomX*ndomY x 1

DOMAINVAR = DefaultField(DOMAINVAR,'ROTATIONS',[]) ;
if ~isempty(DOMAINVAR.ROTATIONS)
    DOMAINVAR.ROTATIONS = DOMAINVAR.ROTATIONS(:) ; % Rotation for each domain
else
    DOMAINVAR.ROTATIONS  = cell(size(DOMAINVAR.ListElements)) ;
end
% cell array
for idom = 1:length(DOMAINVAR.ListElements)
    elemLOC = DOMAINVAR.ListElements{idom} ; % List of elements, old numbering
    DOMAINVAR.ListElements{idom} = IndicesRenumberingElements_INV(elemLOC) ;
end




% After this renumbering, we are sure that
% DOMAINVAR.ListElements{idom}(jelem) is paired with
% DOMAINVAR.ListElements{jdom}(jelem) for all jelem
% --------------------------------------------------------------------
NODESrve = DOMAINVAR.ListNodesDom(:) ;   % List of nodes of each domain. They were constructed so that
% NODESrve{idom}(inode) is paired with NODESrve{jdom}(inode)
%
% DATAOUT.TypeElement = TypeElement ;
% DATAOUT.TypeElementB = TypeElementB ;
% DATAOUT.CONNECTb = CONNECTb ;
% DATAOUT.posgp = posgp ;
% DATAOUT.COOR = COOR ;
% DATAOUT.CN  = CN ;
% DATAOUT.DOMAINVAR = DOMAINVAR ;
nDOM = length(DOMAINVAR.ListElements);  % Number of domains
ndim = size(COOR,2) ;
nDOMtotal  =nDOM ;

domainsINCLUDE = 1:nDOM ;  % 11-Jun-2019---Avoid operating on non-included domains
[nDOMx,nDOMy] = size(DOMAINVAR.ListNodesDom) ;
% Domains to be included.
if ~isempty(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})   &&  length(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})<nDOM
    domainsINCLUDE =DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject} ;
    selected_columnsLOC = SelectedComlumnsFun(DATAIN,iproject,nDOMx,nDOMy) ;
    
    %domainsINCLUDE =DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject} ;
    domainsINCLUDE = intersect(1:nDOM,selected_columnsLOC)' ;
    % THE FIRST DOMAIN IS ALWAYS INCLUDED ---AND THEN REMOVED
    domainsINCLUDE = [1,domainsINCLUDE] ;
    REMOVE_FIRST_DOMAIN =1 ;
    nDOM = length(domainsINCLUDE) ;
else
    REMOVE_FIRST_DOMAIN = 0 ;
end





% --------------------------------------------
% Loop over number of DOMAINS
% ----------------------------------------------
%dbstop('45')
if ndim == 3
    nstrain = 6 ;
else
    nstrain = 4;
end


% ------------------------------------------
% TIME STEPS TO INCLUDE  IN THE ANALYSIS , 11-June-2019
% ----------------------------------------
[~,nstep_total ]= size(d) ; % Non-linear , nodal variables .
[~,nstep_totalG ]= size(stressGLO) ; % Stress variables ... The number of stored time steps might be different
nstep = min(nstep_total,nstep_totalG) ;
timestepINCLUDE = 1:nstep;
if ~isempty(DATAIN.TIME_STEPS_TO_INCLUDE_TRAINING{iproject})
    % Select time steps
    nstepINCLUDE= length(DATAIN.TIME_STEPS_TO_INCLUDE_TRAINING{iproject}) ;
    [ nstep,iii ]= min([nstepINCLUDE,nstep]) ;
    if iii==1
        timestepINCLUDE = DATAIN.TIME_STEPS_TO_INCLUDE_TRAINING{iproject} ;
    end
    
end


nnodeDOM = length(NODESrve{1}) ; %number of nodes per domain
nelemDOM = length(DOMAINVAR.ListElements{1}) ; %number of elements per domain
ngaus = size(posgp,2) ;  % Number of Gauss points per element
% dDOM = zeros(nnodeDOM*ndim,nDOM) ; % Matrix to store displacements
% reactDOM =  zeros(nnodeDOM*ndim,nDOM) ; % Matrix to store reactions
% stressDOM =  zeros(nelemDOM*ngaus*nstrain,nDOM) ; % Matrix to store stresses
% Non-linear regime
dDOM = zeros(nnodeDOM*ndim,nDOM*nstep) ; % Matrix to store displacements
reactDOM =  zeros(size(dDOM)) ; % Matrix to store reactions
ndofG = nelemDOM*ngaus*nstrain;
stressDOM =  zeros(ndofG,nDOM*nstep) ; % Matrix to store stresses


DATALOC.ORTHOGONAL_RIGID_BODY_MODES = 0 ; % Normalization of rigid body modes

%% DOFs boundary, for domain 1 (reference)
% ----------------------------------------------------------
bndDOFS = FaceNodesRefDom(DOMAINVAR,ndim) ;
% ----------------------------------------------------------


DATA_REFMESH = [] ;



nDOMx = size(CONNECTb,1) ;
nDOMy = size(CONNECTb,2) ;

for idomLOC  = 1:nDOM
    idom = domainsINCLUDE(idomLOC) ; % 11-June-2019
    disp(['Domain = ',num2str(idom), ' of ',num2str(nDOMtotal)])
    [idomX,idomY] = ind2sub([nDOMx,nDOMy],idom) ;
    % ------------------------------------------------
    % IDENTIFICATION elements and nodes of domain "i"
    % ------------------------------------------------
    COORabs = COOR(NODESrve{idom},:) ;  % Coordinates of DOmain i
    DOFs = small2large(NODESrve{idom},ndim) ;  % Degrees of Freedom
    dLOC = d(DOFs,timestepINCLUDE) ;  % Displacements of domain i
    
    
    
    
    % FpntLOC = Fpnt(DOFs) ;
    % ------------------------------------------------------------------------------
    % Purging rotations  and translation from displacement vector (main ouput dDOM)
    %-------------------------------------------------------------------------------
    % Construct basis matrix from centroid  (columns are normalized, )
    % ------------------------------------
    GEOMETRIC_PROPERTIES_VOLUME = [] ;
    % STEP 2: List of Gauss Points of domain idom
    eDOM = DOMAINVAR.ListElements{idom}  ; % List of elements domain idom
    ListGauss =  small2large(eDOM,ngaus) ;
    idomTOT= small2large(idomLOC,nstep) ; % INDEXES snapshots where to store the solution
    
    %     if DATAIN.ORTHOGONALITY_RIGIDB_MASS_MATRIX == 0
    %         BasisUrb = ConstructBasisRigidBody_centroid(COORabs,DATALOC) ; % Basis matrix for rigid body motions, relative to centroid
    %         % Purging rotations
    %         dDOM(:,idom) = dLOC -BasisUrb*(BasisUrb\dLOC) ;
    %
    %     else
    % Mass-matrix orthogonal (geometric mass matrix)
    ListGaussDOF =  small2large(ListGauss,ndim) ;
    %if idom == 1
    if ~isempty(Nst)
        Ndom = Nst(ListGaussDOF,DOFs) ;
        
    else
        load(DATAIN.NAME_WS_MODES,'Ndom') ;
        % In some cases, Nst is not stored in memory, and the global
        % matrix is assumed to be the same as the reference one
    end
    WdomN = wSTs_RHS(ListGauss)  ;
    %   else
    %      Ndom = Nst ;
    
    % end
    [BasisUrb,VOLUME,Rbar ]= ConstructBasisRigidBody_MASSM(COORabs,Ndom,WdomN,DATALOC) ; %
    % Purging rotations
    GEOMETRIC_PROPERTIES_VOLUME = (BasisUrb'*Rbar) ;  % Volume and moment of inertias
    
    if all(GEOMETRIC_PROPERTIES_VOLUME == 0)
        disp('When using triangles, you have to emply the same number of GAUSS points for  M and K')
        error('Something is wrong with Nst')
        
    end
    
    dLOCnew = dLOC -BasisUrb*(GEOMETRIC_PROPERTIES_VOLUME\(Rbar'*dLOC)) ;
    %   end
    %% ----------------------------------------------------------------------
    
    %%%%%
    ROTATION =   DOMAINVAR.ROTATIONS{idom} ;  % Maps local to global
    
    if ~isempty(ROTATION)  %  Rotation, several time steps 11-June-2019
        for istep = 1:nstep
            dLOCnew_istep = reshape(dLOCnew(:,istep),ndim,[]) ;
            dLOCnew_istep = (ROTATION'*dLOCnew_istep);
            dDOM(:,idomTOT(istep))  = dLOCnew_istep(:) ;
        end
    else
        dDOM(:,idomTOT)  = dLOCnew ;
    end
    
    %
    %     if ~isempty(ROTATION)
    %         dLOCnew = reshape(dLOCnew,ndim,[]) ;
    %         dLOCnew = (ROTATION'*dLOCnew);
    %         dDOM(:,idom)  = dLOCnew(:) ;
    %     else
    %         dDOM(:,idom)  = dLOCnew ;
    %     end
    
    
    
    %-----------------------------------------
    % Determining reactions at each domain
    % ----------------------------------------
    
    
    
    
    % Determining stresses and reactions
    CONNECTb_local = CONNECTb(idomX,idomY,:) ;
    CONNECTb_local = CONNECTb_local(:) ;
    DATAIN.ROTATION_LOCAL = ROTATION;
    [stressDOM(:,idomTOT),reactDOM(:,idomTOT),Ndom,WdiagRHS,Bdom,Wdom,K,CgloDOM] = ...
        ReactionsRVE(Nst,fNOD,ListGauss,ndim,DOFs,Fpnt,CNb,...
        NODESrve{idom} ,Tnod,TypeElementB,CONNECTb_local,Bst,dLOC,wSTs_RHS,wSTs,Cglo,COOR,...
        nnodeDOM,bndDOFS,BasisUrb,nstrain,idom,iproject,DATAIN,stressGLO(:,timestepINCLUDE),...
        DATA_INPUT_FE,timestepINCLUDE)  ;
    
    if iproject == 1 && idom == 1
        DATA_REFMESH.BasisUrb = BasisUrb ; % Rigid body modes
        DATA_REFMESH.Nst = Ndom ; % Matrix of shape functions
        DATA_REFMESH.WdiagRHS = WdiagRHS ; % Matrix of weights right-hand side
        % DATA_REFMESH.Bst = Bdom ;  % Strain-displacmeent matrix
        % DATA_REFMESH.Wdom = Wdom ;  % Weights left-hand side
        %  DATA_REFMESH.K = K ;  % Stiffness matrix
        %  DATA_REFMESH.CgloDOM = CgloDOM ;  % Stiffness matrix
        DATA_REFMESH.posgp_RHS = posgp_RHS ;
        DATA_REFMESH.GEOMETRIC_PROPERTIES_VOLUME =GEOMETRIC_PROPERTIES_VOLUME ;
        disp('Saving project variables....')
        tic
        %         Bst_old = Bst;
        %        Bst = Bdom ;
        
        save(DATAIN.NAME_WS_MODES,'Ndom','WdiagRHS','Bdom','Wdom','K','CgloDOM')
        toc
        disp('Done')
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA REFERENCE DOMAIN --> First domain, first project
% Variables to be stored in memory
% Name mesh file, nodes faces 1 and 2, length, coordinates, connectivities,
% type of elements...
%% ..................................................
if iproject == 1
    DATA_REFMESH.AREA = AREA;
    DATA_REFMESH =  PropertiesReferenceRVE(DATA_REFMESH,NameFileMesh,DOMAINVAR,COOR,...
        CN,CONNECTb,TypeElement,TypeElementB,MaterialType,posgp,density,DATAIN,NameFileMeshLOC_coarse)  ;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

