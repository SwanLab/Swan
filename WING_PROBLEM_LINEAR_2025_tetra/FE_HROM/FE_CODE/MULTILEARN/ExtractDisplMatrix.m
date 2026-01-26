function [dDOM,DATAOUT ]= ExtractDisplMatrix(INPUT_DATAFILE,DATAIN)

% Output: Matrix of displacements (dDOM), 
% DATAOUT: Structure containing mesh information of each domain
if nargin == 0
    load('tmp2.mat')
end 
 
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
DATAOUT.CN  = CN 
DATAOUT.DOMAINVAR = DOMAINVAR ;
nDOM = length(DOMAINVAR.ListElements);
% Store the coordinates, connectivitites and displacememts of each cell in
% a cell array
ndim = size(COOR,2) ;
%d = reshape(d,ndim,[])' ; % Displacement vector


% --------------------------------------------
% Loop over number of DOMAINS
% ----------------------------------------------
%dbstop('45')
nnodeDOM = length(NODESrve{1}) ;
dDOM = zeros(nnodeDOM*ndim,nDOM) ;
 DATALOC.ORTHOGONAL_RIGID_BODY_MODES = 2 ; % Normalization of rigid body modes 
for i  = 1:nDOM
    disp(['Domain = ',num2str(i), ' of ',num2str(nDOM)])
    % ------------------------------------------------
    % IDENTIFICATION elements and nodes of domain "i"
    % ------------------------------------------------
    COORabs = COOR(NODESrve{i},:) ;  % Coordinates of DOmain i
    DOFs = small2large(NODESrve{i},ndim) ;  % Degrees of Freedom
    dLOC = d(DOFs) ;  % Displacements of domain i
    % ------------------------------------------------------------------------------
    % Purging rotations  and translation from displacement vector (main ouput dDOM)
    %-------------------------------------------------------------------------------
    % Construct basis matrix from centroid  (columns are normalized, )
    % ------------------------------------
    R = ConstructBasisRigidBody_centroid(COORabs,DATALOC) ; % Basis matrix for rigid body motions, relative to centroid
    % Purging rotations
    dDOM(:,i) = dLOC - R*(R\dLOC) ;
    %
    
end


