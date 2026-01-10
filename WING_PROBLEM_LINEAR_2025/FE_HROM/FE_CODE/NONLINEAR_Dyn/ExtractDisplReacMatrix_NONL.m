function [dDOM,reactDOM,stressDOM,DATA_REFMESH,nstep,nstepG ]= ...
    ExtractDisplReacMatrix_NONL(INPUT_DATAFILE,DATAIN,iproject)

% Output: Matrix of displacements (dDOM), and matrix of reactions (reactDOM)
% DATAOUT: Structure containing mesh information of each domain
if nargin == 0
    load('tmp1.mat')
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
            'wSTs_RHS','posgp_RHS','Fpnt','Tnod','CNb','wSTs','Bst','NameFileMesh','density','AREA',...
            'NameFileMeshLOC_coarse') ;
        load(INPUT_DATAFILE,'d','fNOD','Fpnt','Tnod','CNb','stressGLO','Cglo') ; %
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
            load(INPUT_DATAFILE,'Bst','Nst','stressGLO','Cglo') ;
        end
    end
else
    error(['Non existing FE data !!! '])
end


DATAIN = DefaultField(DATAIN,'ORTHOGONALITY_RIGIDB_MASS_MATRIX',1) ;


%%%%%%%%%%
% IndicesRenumberingElements is created when changing the numbering of
%  % elements
% ELASTOSTATIC_GEN/FE_CODE/SolveElastFE.m
% IndicesRenumberingElements
[dummy IndicesRenumberingElements_INV] = sort(IndicesRenumberingElements) ;
for idom = 1:length(DOMAINVAR.ListElements)
    elemLOC = DOMAINVAR.ListElements{idom} ; % List of elements, old numbering
    DOMAINVAR.ListElements{idom} = IndicesRenumberingElements_INV(elemLOC) ;
end




% After this renumbering, we are sure that
% DOMAINVAR.ListElements{idom}(jelem) is paired with
% DOMAINVAR.ListElements{jdom}(jelem) for all jelem
% --------------------------------------------------------------------
NODESrve = DOMAINVAR.ListNodesDom ;   % List of nodes of each domain. They were constructed so that
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
domainsINCLUDE = 1:nDOM ;

if ~isempty(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})   &&  length(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})<nDOM
    domainsINCLUDE =DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject} ;
    nDOM = length(domainsINCLUDE) ;
end
ndim = size(COOR,2) ;
% --------------------------------------------
% Loop over number of DOMAINS
% ----------------------------------------------
%dbstop('45')
if ndim==3
    nstrain = 6 ;
else
    nstrain = 4;
end
nnodeDOM = length(NODESrve{1}) ; %number of nodes per domain
nelemDOM = length(DOMAINVAR.ListElements{1}) ; %number of elements per domain
ngaus = size(posgp,2) ;  % Number of Gauss points per element
[~,nstep ]= size(d) ;
ndof = nnodeDOM*ndim ;
dDOM = zeros(ndof,nDOM*nstep) ; % Matrix to store displacements
reactDOM =  zeros(size(dDOM)) ; % Matrix to store reactions
[~,nstepG ]= size(stressGLO) ;
ndofG = nelemDOM*ngaus*nstrain;
stressDOM =  zeros(ndofG,nDOM*nstepG) ; % Matrix to store stresses
DATALOC.ORTHOGONAL_RIGID_BODY_MODES = 0 ; % Normalization of rigid body modes

 


%% DOFs boundary, for domain 1 (reference)
idom = 1;
faceNODES_1  = DOMAINVAR.NODES_faces12{idom,1} ;
faceNODES_2  = DOMAINVAR.NODES_faces12{idom,2} ;
faceDOFS_1 = small2large(faceNODES_1,ndim) ;
faceDOFS_2 = small2large(faceNODES_2,ndim) ;


DATA_REFMESH = [] ;



% Therefore, the total list of interface DOFs of domain idom is given
% by
faceDOFS = [faceDOFS_1; faceDOFS_2] ;
for idomLOC  = 1:nDOM
    
    idom = domainsINCLUDE(idomLOC) ;
    
    disp(['Domain = ',num2str(idom)])
    % ------------------------------------------------
    % IDENTIFICATION elements and nodes of domain "i"
    % ------------------------------------------------
    COORabs = COOR(NODESrve{idom},:) ;  % Coordinates of DOmain i
    DOFs = small2large(NODESrve{idom},ndim) ;  % Degrees of Freedom
    dLOC = d(DOFs,:) ;  % Displacements of domain i
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
    %     if DATAIN.ORTHOGONALITY_RIGIDB_MASS_MATRIX == 0
    %         BasisUrb = ConstructBasisRigidBody_centroid(COORabs,DATALOC) ; % Basis matrix for rigid body motions, relative to centroid
    %         % Purging rotations
    %         dDOM(:,idomLOC) = dLOC -BasisUrb*(BasisUrb\dLOC) ;
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
    
    idomTOT = small2large(idomLOC,nstep) ;
    
    
    dDOM(:,idomTOT) = dLOC -BasisUrb*(GEOMETRIC_PROPERTIES_VOLUME\(Rbar'*dLOC)) ;
    %  end
    %% ----------------------------------------------------------------------
    
    %-----------------------------------------
    % Determining reactions at each domain
    % ----------------------------------------
    
    
    
    
    % Determining stresses and reactions
    idomREACT = small2large(idomLOC,nstep) ;
    idomSTRESS = small2large(idomLOC,nstepG) ;
    [stressDOM(:,idomSTRESS),reactDOM(:,idomREACT),Ndom,WdiagRHS,Bdom,Wdom,K,CgloDOM] = ...
        ReactionsDOMnonl(Nst,fNOD,ListGauss,ndim,DOFs,Fpnt,CNb,...
        NODESrve{idom} ,Tnod,TypeElementB,CONNECTb,Bst,wSTs_RHS,wSTs,COOR,...
        nnodeDOM,faceDOFS,BasisUrb,nstrain,idom,iproject,DATAIN,stressGLO,DATA_INPUT_FE,Cglo)  ;
    
    if iproject == 1 && idomLOC == 1
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
    DATA_REFMESH =  PropertiesReferenceDomain(DATA_REFMESH,NameFileMesh,DOMAINVAR,COOR,faceNODES_1,faceNODES_2,...
        CN,CONNECTb,TypeElement,TypeElementB,MaterialType,posgp,density,DATAIN,NameFileMeshLOC_coarse)  ;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

