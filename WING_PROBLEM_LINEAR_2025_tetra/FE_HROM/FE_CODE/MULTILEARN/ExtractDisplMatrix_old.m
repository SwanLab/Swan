function [dRVE,NODESfaces,IDXglo,BasisRrb,NODESfaces_orig,DATAOUT,MaterialTypeLOC,ndim,DOFS_reference ]= ...
    ExtractDisplMatrix_old(INPUT_DATAFILE,DATA,DATAIN)

if nargin == 0
    load('tmp2.mat')
end
%
% INCLUDE_MODES_RVE = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COORrve = {} ;
CNrve = {}  ;
dRVE = {}  ;
NODESrve = {}  ;
NODESfaces = {}  ;
RIGID_BODY_MOTIONglo = {} ;
CoordinatesChangeGLO = {} ;
%reactDOM = {} ; % REactions of each domain
IDXelemGLO = {} ;
IDXgaussGLO = {} ;
IndElements = {} ;
% -------------------------
% (EXTRACTION OF PROPERTIES FROM FE FILES)
% -------------------------
 load(INPUT_DATAFILE,'DATA_INPUT_FE') ;  % DATA_INPUT_FE --> Input data file 
 if exist(INPUT_DATAFILE,'file')
     LIST_VAR = who('-file',INPUT_DATAFILE) ;
     if ~ismember(LIST_VAR,'COOR' )
         % This means that the COORDINATE matrix is not stored in this mat
         % file.
         % Then the required information is in
         load(DATA_INPUT_FE.nameWORKSPACE_Kstiff,'MaterialType','COOR','CN','TypeElement','posgp',...
             'CNb','TypeElementB','DOMAINVAR','CONNECTb','IndicesRenumberingElements') ;
         load(INPUT_DATAFILE,'d') ; % Then the displecement field is loaded from the actual file of the project       
     else
         load(INPUT_DATAFILE,'MaterialType','COOR','CN','d','TypeElement','posgp',...
             'CNb','TypeElementB','DOMAINVAR','CONNECTb','IndicesRenumberingElements') ;
     end
 else
     error(['Non existing FE data !!! '])
 end
 
 
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
NODESrve = DOMAINVAR.ListNodesDom ; 


DATAOUT.TypeElement = TypeElement ;
DATAOUT.TypeElementB = TypeElementB ;
DATAOUT.CONNECTb = CONNECTb ;
DATAOUT.posgp = posgp ;
DATAOUT.COOR = COOR ;
nRVE = length(DOMAINVAR.ListElements); 
% Store the coordinates, connectivitites and displacememts of each cell in
% a cell array
ndim = size(COOR,2) ;
%d = reshape(d,ndim,[])' ; % Displacement vector


% --------------------------------------------
% Loop over number of DOMAINS
% ----------------------------------------------
%dbstop('45')
for i  = 1:nRVE
    disp(['Domain = ',num2str(i), ' of ',num2str(nRVE)])
    % ------------------------------------------------
    % IDENTIFICATION elements and nodes of domain "i"
    % ------------------------------------------------
    COORabs = COOR(NODESrve{i},:) ;  % Coordinates of DOmain i
    DOFs = small2large(NODESrve{i},ndim) ; 
    dRVE{i} = d(DOFs) ;  % Displacements of domain i
    % ------------------------------------------------------------------------------
    % Purging rotations  and translation from displacement vector (main ouput dRVE)
    %-------------------------------------------------------------------------------
    % Construct basis matrix from centroid 
    % ------------------------------------
    R = ConstructBasisRigidBody_centroid(COORabs) ; % Basis matrix for rigid body motions, relative to centroid
    dRVE{i} = dRVE{i} - R*(R\dRVE{i}) ; 
    
    [NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,COORrve,DOFs_ref] ...
        = PurgingRotations(i,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
        COORrve,DATAIN);
    if i==1
        DOFS_reference= DOFs_ref ;
    end
    % ---------------------------------------------------------------
    % Computing self-equilibrated reactions for each domain
    % -------------------------------------------------------
    %     [reactDOMloc] = ReactionsRVE(i,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
    %         COORrve,Nst,wSTs_RHS,fNOD,NODESrve,CNrve,posgp,IndElements{i},posgp_RHS,...
    %         CNb,TypeElementB,Fpnt,Tnod,CONNECTb,COOR,stressGLO,wSTs,Bst,TypeElement);
    %     reactDOM{i} =  reactDOMloc  ;
end

NODESfaces_orig = NODESfaces ;

NCELLS = length(dRVE) ;

%%% REFERENCE MESH.
refMESH = 1 ;
%% COORDINATES OF REFERENCE MESH
COORref = COORrve{refMESH} ;
DATAOUT.COORref = COORref ;
DATAOUT.CNref = CNrve{refMESH} ;
DATAOUT.NODESref = NODESrve{refMESH} ;

COORrefORD = reshape(COORref',prod(size(COORref)),1) ;
% Now we re-order
COORnew = {} ;
ElementsREF = IndElements{refMESH} ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identification of similar nodes between domains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(COORrve)
    if i ~= refMESH
        IDX = knnsearch(COORrve{i},COORref) ;   % So that COORrve{i}(IDX(j),:) = COORref(j,:)
        % See Implementation.pdf. It turns out there is no need to obtain
        % IDXelem, and IDXgauss
        %          if  DATA.CUBATURE.ACTIVE == 1
        %              % Mapping between meshes
        %              dbstop('102')
        %              TypeIntegrand = 'K' ;
        %              IDX_or = NODESrve{i}(IDX) ;
        %              [IDXelem,IDXgauss ]= EquivalenceBetweenMeshes(COORref,IDX_or,CNrve{i},CNrve{refMESH},TypeIntegrand,...
        %                  TypeElement) ;
        %              IDXelemGLO{i} = IDXelem ;
        %              IDXgaussGLO{i} = IDXgauss ;
        %          end
        
        IDXglo{i} = IDX ;
        COORrve{i} =  COORrve{i}(IDX,:) ;
        CoordinatesChange=  CoordinatesChangeGLO{i}(IDX,:) ;
        CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChange)),1 );
        COORnew{i} =COORrefORD + CoordinatesChange ;
        
        %    RIGID_BODY_MOTIONglo{i} =  RIGID_BODY_MOTIONglo{i}(IDX,:) ;
        %    RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
        dRVE{i} =  dRVE{i}(IDX,:)  ;
        dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
        IDXdof = small2large(IDX,ndim) ;
        %      reactDOM{i} = reactDOM{i}(IDXdof) ;
        for j = 1:length(NODESfaces{i})
            NODESfacesALL = zeros(size(NODESrve{i})) ;
            NODESfacesALL(NODESfaces{i}{j}) = 1 ;
            NODESfacesALL = NODESfacesALL(IDX) ;
            NODESfaces{i}{j}  =find(NODESfacesALL == 1) ;
        end
        NODESrve{i} =  NODESrve{i}(IDX) ;
    else
        dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
        %   RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
        CoordinatesChange=  CoordinatesChangeGLO{i}  ;
        CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChangeGLO{i})),1 );
        COORnew{i} =COORrefORD + CoordinatesChange ;
    end
end

%%%%%%%%%%%%%%%%%
%%% RIGID BODY MODES reference mesh
% -----------------
dRVE = cell2mat(dRVE) ;
COORrel = COORref ;
BasisRrb = ConstructBasisRigidBody(COORrel,  DATAIN) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
