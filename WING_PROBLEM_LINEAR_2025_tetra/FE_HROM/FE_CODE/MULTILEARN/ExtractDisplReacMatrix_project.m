function [dDOM,reactDOM,stressDOM,DATA_REFMESH,REMOVE_FIRST_DOMAIN,timestepINCLUDE ]= ...
    ExtractDisplReacMatrix_project(INPUT_DATAFILE,DATAIN,iproject)

% Output: Matrix of displacements (dDOM), and matrix of reactions (reactDOM)
% DATAOUT: Structure containing mesh information of each domain
if nargin == 0
    load('tmp4.mat')
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
        % load(INPUT_DATAFILE,'d','fNOD','Fpnt','Tnod','CNb') ; %
        load(INPUT_DATAFILE,'d','fNOD','Fpnt','Tnod','CNb','stressGLO') ; % 27-May-2019,nonlinear
        
        
        
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
[dummy, IndicesRenumberingElements_INV] = sort(IndicesRenumberingElements) ;
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
nDOMtotal  =nDOM ;
%%%%%%%%%%%%%%%%%%%%


domainsINCLUDE = 1:nDOM ;  % 27-May-2019---Avoid operating on non-included domains
% Domains to be included.
if ~isempty(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})   &&  length(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})<nDOM
    domainsINCLUDE =DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject} ;
    domainsINCLUDE = intersect(1:nDOM,domainsINCLUDE) ;
    % THE FIRST DOMAIN IS ALWAYS INCLUDED ---AND THEN REMOVED
    domainsINCLUDE = [1,domainsINCLUDE] ;
    REMOVE_FIRST_DOMAIN =1 ;
    nDOM = length(domainsINCLUDE) ;
else
    REMOVE_FIRST_DOMAIN = 0 ;
end

% ------------------------------------------
% TIME STEPS TO INCLUDE  IN THE ANALYSIS , 27-mAY-2019
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



%%%%%%%%%%%%%%%%%%%%%%


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
ngaus_RHS = size(posgp_RHS,2) ;  % 22-Apr-2021 


dDOM = zeros(nnodeDOM*ndim,nDOM*nstep) ; % Matrix to store displacements
reactDOM =  zeros(size(dDOM)) ; % Matrix to store reactions
ndofG = nelemDOM*ngaus*nstrain;
stressDOM =  zeros(ndofG,nDOM*nstep) ; % Matrix to store stresses
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

DOMAINVAR = DefaultField(DOMAINVAR,'rotGLO',cell(1,nDOM)) ;

for idomLOC  = 1:nDOM
    
    idom = domainsINCLUDE(idomLOC) ;
    
    disp(['Domain = ',num2str(idom),' of ',num2str(nDOMtotal)]) ;
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
        ListGauss_RHS =  small2large(eDOM,ngaus_RHS) ; % 22-Apr-2021

    idomTOT= small2large(idomLOC,nstep) ; % INDEXES snapshots where to store the solution
    
    %     if DATAIN.ORTHOGONALITY_RIGIDB_MASS_MATRIX == 0
    %         error('Old option. Set DATAIN.ORTHOGONALITY_RIGIDB_MASS_MATRIX  ==1')
    %         BasisUrb = ConstructBasisRigidBody_centroid(COORabs,DATALOC) ; % Basis matrix for rigid body motions, relative to centroid
    %         % Purging rotations
    %         dDOM(:,idomTOT) = dLOC -BasisUrb*(BasisUrb\dLOC) ;
    %
    %     else
    %
    % Mass-matrix orthogonal (geometric mass matrix)
    ListGaussDOF =  small2large(ListGauss,ndim) ;
    ListGaussDOF_RHS =  small2large(ListGauss_RHS,ndim) ; % 22-Apr-2021
    %if idom == 1
    if ~isempty(Nst)
        Ndom = Nst(ListGaussDOF_RHS,DOFs) ;  % 22-Apr-2021
        
    else
        load(DATAIN.NAME_WS_MODES,'Ndom') ;
        % In some cases, Nst is not stored in memory, and the global
        % matrix is assumed to be the same as the reference one
    end
    WdomN = wSTs_RHS(ListGauss_RHS)  ;
    %   else
    %      Ndom = Nst ;
    
    % end
    [BasisUrb,VOLUME,Rbar ]= ConstructBasisRigidBody_MASSM(COORabs,Ndom,WdomN,DATALOC) ; %
    % Purging rotations
    GEOMETRIC_PROPERTIES_VOLUME = (BasisUrb'*Rbar) ;  % Volume and moment of inertias
    dLOCnew = dLOC -BasisUrb*(GEOMETRIC_PROPERTIES_VOLUME\(Rbar'*dLOC)) ;
    
    %%%%%
    if iscell(DOMAINVAR.rotGLO)
        ROTATION =   DOMAINVAR.rotGLO{idom} ;  % Maps local to global
    else
        ROTATION = [] ;
    end
    
    % Rotation must be performed time-step-wise
    
    
    if ~isempty(ROTATION)
        for istep = 1:nstep
            dLOCnew_istep = reshape(dLOCnew(:,istep),ndim,[]) ;
            dLOCnew_istep = (ROTATION'*dLOCnew_istep);
            dDOM(:,idomTOT(istep))  = dLOCnew_istep(:) ;
        end
    else
        dDOM(:,idomTOT)  = dLOCnew ;
    end
    
    
    
    
    
    DATAIN.ROTATION = ROTATION ;
    %   end
    %% ----------------------------------------------------------------------
    
    %-----------------------------------------
    % Determining reactions at each domain
    % ----------------------------------------
    
    
    
    
    % Determining stresses and reactions
    % idomREACT = small2large(idomLOC,nstep) ;
    
    % Change 27-May-2019
    [stressDOM(:,idomTOT),reactDOM(:,idomTOT),Ndom,WdiagRHS,Bdom,Wdom,K,CgloDOM] = ...
        ReactionsDOMAIN(Nst,fNOD,ListGauss,ndim,DOFs,Fpnt,CNb,...
        NODESrve{idom} ,Tnod,TypeElementB,CONNECTb,Bst,dLOC,wSTs_RHS,wSTs,Cglo,COOR,...
        nnodeDOM,faceDOFS,BasisUrb,nstrain,idom,iproject,DATAIN,stressGLO(:,timestepINCLUDE),...
        DATA_INPUT_FE,timestepINCLUDE,ListGauss_RHS)  ;
    
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
    DATA_REFMESH =  PropertiesReferenceDomain(DATA_REFMESH,NameFileMesh,DOMAINVAR,COOR,faceNODES_1,faceNODES_2,...
        CN,CONNECTb,TypeElement,TypeElementB,MaterialType,posgp,density,DATAIN,NameFileMeshLOC_coarse)  ;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

