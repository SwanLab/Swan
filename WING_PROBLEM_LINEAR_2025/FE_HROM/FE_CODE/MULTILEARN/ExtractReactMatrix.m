function [reactDOM stressDOM BdomRED Wdom]= ...
    ExtractReactMatrix(INPUT_DATAFILE,DATA,dRVE_strain,BasisUdef,IDXglo,...
    BasisUrb,NODESfaces,DATACUBATURE,COLUMNS_RVE,COLUMNS_RVEloc,DATAIN)
% See Implementation.pdf
% INCLUDE_MODES_RVE = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    load('tmp1.mat')
end
COORrve = {} ;
CNrve = {}  ;
%dRVE = {}  ;
NODESrve = {}  ;
RIGID_BODY_MOTIONglo = {} ;
CoordinatesChangeGLO = {} ;
nnode  = size(BasisUrb,1);
reactDOM = zeros(nnode,size(dRVE_strain,2))  ; % REactions of each domain
stressDOM = [] ;
IndElements = {} ;
% -------------------------
% (EXTRACTION OF PROPERTIES FROM FE FILES)
% -------------------------
NAME_MODES =[INPUT_DATAFILE,'.msh'];
outputFE = ['DATAWS/',INPUT_DATAFILE,'_WS.mat'] ;
eval(INPUT_DATAFILE);%
load(outputFE,'MaterialType','COOR','CN','TypeElement','posgp','Nst','wSTs_RHS','fNOD','posgp_RHS',...
    'CNb','TypeElementB','Fpnt','Tnod','CONNECTb','Cglo','wSTs','Bst','d') ;

% Number of materials/bodies
IND_MAT =unique(MaterialType) ;
% Number of original materiasl
nMAT = length(MATERIAL.PLY) ;
% Number of RVEs
nRVE = size(dRVE_strain,2)        ;
% Store the coordinates, connectivitites and displacememts of each cell in
% a cell array
ndim = size(COOR,2) ;
% Tolerance for matching nodes %
%TOL = ChooseTolerance(CN,COOR) ;
% --------------------------------------------
% Loop over number of DOMAINS
% ----------------------------------------------
ngausDOF = size(Bst,1)/nRVE ;
stressDOM = zeros(ngausDOF,size(dRVE_strain,2))  ; % REactions of each domain
nTYPErve = length(COLUMNS_RVE) ;
DATAIN = DefaultField(DATAIN,'REACTION_MODES_FROM_FE_SIMULATIONS',0)  ;


for i  = 1:nRVE
    disp(['Domain = ',num2str(i), ' of ',num2str(nRVE)])
    % ------------------------------------------------
    % IDENTIFICATION elements and nodes of domain "i"
    % ------------------------------------------------
    IndElementsLOC = [] ;
    for imat = 1:nMAT
        IMAT = (i-1)*nMAT + imat ;
        IndElementsLOC =[IndElementsLOC;  find(IND_MAT(IMAT) == MaterialType)] ;   % Elements of domain i
    end
    IndElementsLOC = sort(IndElementsLOC) ;
    IndElements{i} = IndElementsLOC ;
    CNrve{i} =  CN(IndElements{i},:) ; % Connectivities of Domain i
    NODESrve{i} =unique( CNrve{i}(:)) ;  % Nodes of Domain i
    COORabs = COOR(NODESrve{i},:) ;  % Coordinates of DOmain i
    [nnode ndim] = size(COORabs) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % With no truncation
    % -----------------------
    % dRVE{i} = d( NODESrve{i},:) ;  % Displacements of domain i
    % With truncation
    %
    % dbstop('71')
    if DATAIN.REACTION_MODES_FROM_FE_SIMULATIONS == 0
        %  dRVEloc = BasisUdef{COLUMNS_RVEloc(i)}*(BasisUdef{COLUMNS_RVEloc(i)}\dRVE_strain(:,i)) ; % Projection
        dRVEloc = BasisUdef{1}*(BasisUdef{1}\dRVE_strain(:,i)) ; % Revise this !!!
    else
        dRVEloc = dRVE_strain(:,i) ;
    end
    %  dRVEloc =  reshape(dRVEloc,size(COORabs,2),[])' ;% REshaping
    
    if   ~isempty(IDXglo) && ~isempty(IDXglo{i} )
        IDX = IDXglo{i} ;
        [dummy, IDXinv ]= sort(IDX) ;
        IDXinvDOFs = Nod2DOF(IDXinv,size(COORabs,2)) ;
        %   dbstop('71')
        BasisUrbLOC = BasisUrb(IDXinvDOFs,:) ; % We transform it back to the original configuration
        dRVEloc =dRVEloc(IDXinvDOFs);
        IDXDOFs = Nod2DOF(IDX,size(COORabs,2)) ;
    else
        BasisUrbLOC = BasisUrb ;  IDX = [] ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     % -----------------------------------------------------------
    %     %  NODESfaces, relative coordintates...
    %     %------------------------------------------------------------
    %       [NODESfaces,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,COORrve] ...
    %           = RotationModes(i,DATA,COORabs,TOL,NODESfaces,dRVE,RIGID_BODY_MOTIONglo,CoordinatesChangeGLO,...
    %           COORrve);
    % ---------------------------------------------------------------
    % Computing self-equilibrated reactions for each domain
    % -------------------------------------------------------
    DATA.MakeMeshByRepetition = DefaultField(DATA.MakeMeshByRepetition,'DIRECTION',[]);
    %   dbstop('104')
    [reactDOMloc, stressDOMloc,BdomLOC,WdomLOC] = ReactionsRVEnew(DATA,NODESfaces{i},dRVEloc,...
        nnode,ndim,Nst,wSTs_RHS,fNOD,NODESrve{i},CNrve,posgp,IndElements{i},posgp_RHS,...
        CNb,TypeElementB,Fpnt,Tnod,CONNECTb,COOR,wSTs,Bst,TypeElement,BasisUrbLOC,Cglo,d,...
        DATACUBATURE,DATAIN);
    %   dbstop('91')
    if i==1 % This is the reference mesh; accordingly, we take the strain-displacement basis matrix from this one
        BdomRED ={} ;
        for itype = 1:length(BasisUdef)
            BdomRED{itype} = BdomLOC*BasisUdef{itype};
        end
        Wdom = WdomLOC ;
    end
    
    if ~isempty(IDX)
        reactDOM(:,i) =  reactDOMloc(IDXDOFs)  ;
    else
        reactDOM(:,i) =  reactDOMloc  ;
    end
    if ~isempty(stressDOMloc)
        stressDOM(:,i) =  stressDOMloc  ;
    end
end

%
% NCELLS = length(dRVE) ;
%
% %%% REFERENCE MESH.
% refMESH = 1 ;
% %% COORDINATES OF REFERENCE MESH
% COORref = COORrve{refMESH} ;
% COORrefORD = reshape(COORref',prod(size(COORref)),1) ;
% % Now we re-order
% COORnew = {} ;
% ElementsREF = IndElements{refMESH} ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Identification of similar nodes between domains
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:length(COORrve)
%     if i ~= refMESH
%         IDX = knnsearch(COORrve{i},COORref) ;
%         COORrve{i} =  COORrve{i}(IDX,:) ;
%         CoordinatesChange=  CoordinatesChangeGLO{i}(IDX,:) ;
%         CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChange)),1 );
%         COORnew{i} =COORrefORD + CoordinatesChange ;
%
%         %    RIGID_BODY_MOTIONglo{i} =  RIGID_BODY_MOTIONglo{i}(IDX,:) ;
%         %    RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
%         dRVE{i} =  dRVE{i}(IDX,:)  ;
%         dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
%         IDXdof = small2large(IDX,ndim) ;
%         %      reactDOM{i} = reactDOM{i}(IDXdof) ;
%         for j = 1:length(NODESfaces{i})
%             NODESfacesALL = zeros(size(NODESrve{i})) ;
%             NODESfacesALL(NODESfaces{i}{j}) = 1 ;
%             NODESfacesALL = NODESfacesALL(IDX) ;
%             NODESfaces{i}{j}  =find(NODESfacesALL == 1) ;
%         end
%         NODESrve{i} =  NODESrve{i}(IDX) ;
%     else
%         dRVE{i} =  reshape(dRVE{i}',prod(size(dRVE{i})),1) ;
%         %   RIGID_BODY_MOTIONglo{i} =  reshape(RIGID_BODY_MOTIONglo{i}',prod(size(RIGID_BODY_MOTIONglo{i})),1) ;
%         CoordinatesChange=  CoordinatesChangeGLO{i}  ;
%         CoordinatesChange =  reshape(CoordinatesChange',prod(size(CoordinatesChangeGLO{i})),1 );
%         COORnew{i} =COORrefORD + CoordinatesChange ;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
