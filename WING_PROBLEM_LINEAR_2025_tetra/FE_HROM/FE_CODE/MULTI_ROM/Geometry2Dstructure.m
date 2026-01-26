function MESH2D = Geometry2Dstructure(MESH2Dinput,MESH3D,DATAIN)

if nargin == 0
    load('tmp.mat')
end
% -----------------------------------------
% Reading coordinate structure (2D skeleton, linear)
% -----------------------------------------
[MESH2D]= ReadMeshFileStr(MESH2Dinput.NAME,'READ_MATERIAL_COLUMN',1)  ;
[ELEMENTS_SURFACES,~] = DAT_file_2dskeletonGID(MESH2Dinput.NAME);  % 1-July-2019

DATAIN = DefaultField(DATAIN,'RenumberingElements2D',0)  ;  % To diminish the band of the K stiffness matrix

if DATAIN.RenumberingElements2D == 1
    % Change, 24-June-2019
    % Renumbering
    %DATA.RENUMBERED_OUTSIDE = 1 ; % To avoid performing this operation later, when computing stiffness matrix
    [~,IndicesRenumberingElements]  = sort(MESH2D.CN(:,1)) ;
    MESH2D.CN = MESH2D.CN(IndicesRenumberingElements,:) ;
    if ~isempty(MESH2D.MaterialType)
        MESH2D.MaterialType = MESH2D.MaterialType(IndicesRenumberingElements) ;
    end
    
    [~,IndicesRenumberingElements_old ]= sort(IndicesRenumberingElements) ;
    
    if ~isempty(ELEMENTS_SURFACES)   % 1-July-2019
        for isurf = 1:length(ELEMENTS_SURFACES)
            ELEMENTS_SURFACES{isurf} = IndicesRenumberingElements_old(ELEMENTS_SURFACES{isurf}) ; % Surface elements selected by the user (old mesh)
        end
    end
    
end





MESH2D.NAME = MESH2Dinput.NAME ;
MESH2D.PROP = MESH2Dinput.PROP ;
ndim = size(MESH3D.RVES.DATA3D.COOR,2) ;
if size(MESH2D.COOR,2) ==2 & ndim == 3 % 3D
    MESH2D.COOR = [MESH2D.COOR,zeros(size(MESH2D.COOR,1),1)] ;
end

if ~isempty(DATAIN.angDOM) &&  DATAIN.angDOM ~= 0
    % Center of the cylinder, height
    % ------------------------------
    % Points plane y = ymin
    TOL = 1e-10 ;
    y = MESH2D.COOR(:,2)' ;
    ymin = min(y) ;
    INDICESmin = find(abs(y-ymin)<=TOL) ;
    % Pick up 3 points
    COORcirc = MESH2D.COOR(INDICESmin(1:3),[1,3]) ;
    % Fit a circle
    %DATAIN = DefaultField(DATAIN,'CIRCLE') ;
    % DATAIN.CIRCLE = DefaultField(DATAIN.CIRCLE,'RADIUS',[]) ;
    %  DATAIN.CIRCLE = DefaultField(DATAIN.CIRCLE,'CENTER',[]) ;
    
    %   if isempty(DATAIN.CIRCLE.CENTER)
    [CIRCLE.RADIUS,CIRCLE.CENTER] = fit_circle_through_3_points(COORcirc) ;
    disp(['RADIUS CYLINDER 2D MESH = ',num2str(CIRCLE.RADIUS)])
    %     else
    %         CIRCLE  =  DATAIN.CIRCLE  ;
    %      end
    
    
end

%% ELEMENTS SHOULD BE QUADRATIC QUADRILATERALS
% ---------------------------------------------
switch MESH2D.TypeElement
    case 'Quadrilateral'
        nnodes= size(MESH2D.CN,2) ;
        if nnodes ~= 8
            error('2D elements must be Quadratic Quadrilaterals')
        end
    otherwise
        error('2D elements must be Quadratic Quadrilaterals')
end

%% Now we split the matrix of coordinates and connectivities into ``vertices" nodes and midside nodes
% -----------------------------------------------------------------------------------------------------
MESH2D.CNall = MESH2D.CN ; % Quadratic mesh- All 8-nodes  per element
MESH2D.COORall = MESH2D.COOR ;
MESH2D.CNbALL = MESH2D.CNb ;
% Midside nodes
% -------------
COLUMNS = [5:size(MESH2D.CN,2)] ; % Gid's criterion. Columns assigned to midside nodes
[COORmidside,CNmidside,NODESmid] = ExtractCoorSeparate(MESH2D.CN,COLUMNS,MESH2D.COOR) ;
MESH2D.NODESmid = NODESmid ;
MESH2D.COLUMNSmid = COLUMNS ;

PRINT_MESH_MIDSIDES = 1;
if PRINT_MESH_MIDSIDES == 1
    [FFF,NNN,TTT] = fileparts(MESH2Dinput.NAME) ;
    NameFile_msh = [FFF,filesep,'MIDSIDE_MESH.msh'] ;
    IND_ELEM_MESHES = GidMesh2DFE_multi(NameFile_msh,COORmidside,{CNmidside},'',{ones(size(CNmidside,1),1)},{MESH2D.TypeElement},{'Prueba'});
end

 
% -------------------------------------------------------------------------------
%  The order of the connectivities should be such that  the following
%  correspondence between local and global labels should hold
%  FACE xMiN = 1, FACE yMIN =2 . FACE xMAX = 3, FACE yMAX = 4
if isempty(DATAIN.angDOM) || DATAIN.angDOM == 0
    [FACEmidside_CN,xMIN,xMAX ]= ClassificationNodesMidSide(CNmidside,COORmidside) ;
    MESH2D.xMIN = xMIN ;
    MESH2D.xMAX  = xMAX ;
    MESH2D.rotDOM = [] ;
    MESH2D.angDOMall = [] ;  % angle formed by the global x-axis and the local x-axis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECKING THAT LOCAL NUMBERING CORRESPOND TO GLOBAL NUMBERING
    % --------------------------------------------------------------
    [ FACESloc,INFOREP1,INFOREP2 ]= unique(FACEmidside_CN,'rows') ;
    if size(FACESloc,1) > 1
        % Why is this an error ? Because you should do it one by one. Forget
        % vectorization if an error pops up (4-Apr-2019)
        warning('Local numbering (3D) and global numbering (2D) does not coincide')
        FACESorder = zeros(size(FACESloc)) ; % 31-Jan-2020
        
        for irwos = 1:size(FACESloc,1)
            [~, FACESorder(irwos,:)]= sort(FACESloc(irwos,:)) ;
        end
    else
        [~, FACESorder]= sort(FACESloc);  % 5-Apr-2019. Error amended
        
    end
    
    
else
    [FACEmidside_CN, MESH2D.rotDOM,MESH2D.angDOMall]= ClassificationNodesMidSideCylin(CNmidside,COORmidside,CIRCLE,DATAIN) ;
    FACESloc = unique(FACEmidside_CN,'rows') ;
    
    MESH2D.xMIN = [] ;
    MESH2D.xMAX  = [] ;
    FACESorder = FACESloc;
    warning('Check this part of the code !!! ')  % JAHO, 9-April-2019. It seems that
    % there is some issue with local and global numbering
    INFOREP2 = ones(size(CNmidside,1),1) ; 
    
    %     if size(FACESloc,1) ~=1
    %         error('CHECK THE ORDER OF NUMBERING OF THE MIDSIDES')
    %     end
end

CNmidsideNEW = CNmidside ;

for irows = 1:size(FACESorder,1)   % Modification 31-Jan-2020
    FFF = find(INFOREP2==irows) ;
    CNmidsideNEW(FFF,:) = CNmidside(FFF,FACESorder(irows,:)) ;
end

CNmidside = CNmidsideNEW ;

% We have to modify also MESH2D.CNall   % 31-Jan-2020
CNallNEW = MESH2D.CNall ;
for irows = 1:size(FACESloc,1)
    newCOLUMNS = COLUMNS(FACESloc(irows,:)) ;
    FFF = find(INFOREP2==irows) ;
    CNallNEW(FFF,newCOLUMNS) =  MESH2D.CNall(FFF,COLUMNS) ;
end
MESH2D.CNall = CNallNEW ;

% RIGHT_ORDER_COlUMS = zeros(1,4) ;
% for inode = 1:size(FACEmidside_CN,2)
%     DIFF_meth =  FACEmidside_CN(:,inode)-inode ;
%     if any(DIFF_meth)
%         %warning('Local numbering (3D) and global numbering (2D) does not coincide')
%         %warning('Changing local connectivitities')
%         FACE = unique(FACEmidside_CN(:,inode)) ;
%         if length(FACE) ~=1
%             error('Local numbering (3D) and global numbering (2D) does not coincide')
%         end
%         RIGHT_ORDER_COlUMS(FACE) = inode ;
%
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VERTICES nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COLUMNS = [1:4] ;
[COORvert,CNvert,NODESvert] = ExtractCoorSeparate(MESH2D.CN,COLUMNS,MESH2D.COOR) ;
MESH2D.NODESvert = NODESvert ;
MESH2D.COLUMNSvert = COLUMNS ;
%%% Boundary connectivities
% ---------------------------




% ---------------------------------------------
% FACES 4-1: 1, 1-2: 2, 2-3:3 , 3:4, 4
%b --------------------------------------
if  isempty(DATAIN.angDOM) || DATAIN.angDOM == 0
    % this is for classifying the nodes located at the vertices 
    [FACEvert_CN ]= ClassificationNodesVert(CNvert,COORvert,xMIN,xMAX) ;
    [ FACESloc,INFOREP1,INFOREP2 ] = unique(FACEvert_CN ,'rows') ;
    %[~, FACESorder]= sort(FACESloc);
      FACESorder = zeros(size(FACESloc)) ; % 31-Jan-2020
        
        for irwos = 1:size(FACESloc,1)
            [~, FACESorder(irwos,:)]= sort(FACESloc(irwos,:));
        end
    
%          [ FACESloc,INFOREP1,INFOREP2 ]= unique(FACEmidside_CN,'rows') ;
%     if size(FACESloc,1) > 1
%         % Why is this an error ? Because you should do it one by one. Forget
%         % vectorization if an error pops up (4-Apr-2019)
%         warning('Local numbering (3D) and global numbering (2D) does not coincide')
%         FACESorder = zeros(size(FACESloc)) ; % 31-Jan-2020
%         
%         for irwos = 1:size(FACESloc,1)
%             [~, FACESorder(irwos,:)]= sort(FACESloc(irwos,:))
%         end
%     else
%         [~, FACESorder]= sort(FACESloc);  % 5-Apr-2019. Error amended
%         
%     end
    
    
    
    
%     if size(FACESloc,1) > 1
%         error('Local numbering (3D) and global numbering (2D) does not coincide')
%         % 4-Apr-2019: What to do in this case ? Why should them be equal ?
%     end
else
    [FACESloc ]= ClassificationNodesVertCylin(CNvert,COORvert,CIRCLE) ;
    
end

% 
% CNmidsideNEW = CNmidside ;
% 
% for irows = 1:size(FACESorder,1)   % Modification 31-Jan-2020
%     FFF = find(INFOREP2==irows) ;
%     CNmidsideNEW(FFF,:) = CNmidside(FFF,FACESorder(irows,:)) ;
% end
% 
% CNmidside = CNmidsideNEW ;
% 
% % We have to modify also MESH2D.CNall   % 31-Jan-2020
% CNallNEW = MESH2D.CNall ;
% for irows = 1:size(FACESloc,1)
%     newCOLUMNS = COLUMNS(FACESloc(irows,:)) ;
%     FFF = find(INFOREP2==irows) ;
%     CNallNEW(FFF,newCOLUMNS) =  MESH2D.CNall(FFF,COLUMNS) ;
% end
% MESH2D.CNall = CNallNEW ;


%CNvert = CNvert(:,FACESorder) ;  % Before Jan-31-2020

CNvertNEW = CNvert ;

for irows = 1:size(FACESorder,1)   % Modification 31-Jan-2020
    FFF = find(INFOREP2==irows) ;
    CNvertNEW(FFF,:) = CNvert(FFF,FACESorder(irows,:)) ;
end

CNvert = CNvertNEW ;


% We have to modify also MESH2D.CNall  - Before 31th-Jan-2020
% newCOLUMNS = COLUMNS(FACESloc) ;
% MESH2D.CNall(:,newCOLUMNS) =  MESH2D.CNall(:,COLUMNS) ;

% % We have to modify also MESH2D.CNall   % 31-Jan-2020
 CNallNEW = MESH2D.CNall ;
 for irows = 1:size(FACESloc,1)
     newCOLUMNS = COLUMNS(FACESloc(irows,:)) ;
     FFF = find(INFOREP2==irows) ;
     CNallNEW(FFF,newCOLUMNS) =  MESH2D.CNall(FFF,COLUMNS) ;
 end
 MESH2D.CNall = CNallNEW ;




% Vertices Nodes.
MESH2D.CNvert = CNvert ;
MESH2D.COORvert = COORvert ;

% Midside nodes represent faces.
MESH2D.CNmidside = CNmidside ;
MESH2D.COORmidside = COORmidside ;

% DATAIN.ContinuumStructuresWithoutCorners

if isempty(MESH3D(1).RVES.DATA3D.NODES_CORNERS) || DATAIN.ContinuumStructuresWithoutCorners == 1 % IMP-MAY19-A
    % Open cells (or corners removed by option ContinuumStructuresWithoutCorners)
    % ----------
    MESH2D.CN = CNmidside ;
    MESH2D.COOR = COORmidside ;
else
    % Plate elements
    MESH2D.CN = CNvert ;
    MESH2D.COOR = COORvert ;
end

%
% Identification of lines and surfaces
[~,NODES_LINES] = DAT_file_2dskeletonGID(MESH2D.NAME);
MESH2D.NODES_LINESall = NODES_LINES ;
% We may be only interested in midside nodes:

NODES_LINES=  DetectMidSideNodes(NODES_LINES,NODESmid) ;

MESH2D.NODES_LINES = NODES_LINES ;
MESH2D.ELEMENTS_SURFACES = ELEMENTS_SURFACES ;




end

function [COORmidside,CNmidsideNEW,MDNODES] = ExtractCoorSeparate(CN,COLUMNS,COOR)

CNmidside =  CN(:,COLUMNS) ;
% List of midside nodes
MDNODES = unique(CNmidside(:)) ;
% COORDINATES
COORmidside =  COOR(MDNODES,:) ;
% Renumber connectivities
CNmidsideNEW = RenumberConnectivities(CNmidside,1:length(COORmidside)) ;



end

function  NODES_LINES=  DetectMidSideNodes(NODES_LINES,NODESmid)
for iline = 1:length(NODES_LINES)
    SETNODES = NODES_LINES{iline} ;
    NODE_LINE_MIDSIDE = [] ;
    for innode = 1:length(SETNODES)
        INDD =  find(SETNODES(innode) == NODESmid) ;
        if ~isempty(INDD)
            NODE_LINE_MIDSIDE(end+1) = INDD ;
        end
    end
    NODES_LINES{iline} = NODE_LINE_MIDSIDE ;
    
end

end


function   [FACES_CN,xMIN,xMAX ]= ClassificationNodesMidSide(CNmidside,COORmidside)

% Classification of midside nodes according to the faces they belong to
% ----------------------------------------------------------------------
xMIN = 1e20*ones(size(CNmidside,1),2) ;
xMAX =-1e20*ones(size(CNmidside,1),2) ;
for inodeE = 1:size(CNmidside,2)
    COLUMN_loc = CNmidside(:,inodeE) ;
    COOR_loc = COORmidside(COLUMN_loc,:) ;
    for idim=1:2
        xMIN(:,idim) = min([xMIN(:,idim), COOR_loc(:,idim)],[],2) ;
        xMAX(:,idim) = max([xMAX(:,idim), COOR_loc(:,idim)],[],2) ;
    end
end

FACES_CN = zeros(size(CNmidside)) ;

for inodeE = 1:size(CNmidside,2)  % Loop over the 4 sides of the quadrilateral
    COLUMN_loc = CNmidside(:,inodeE) ;
    COOR_loc = COORmidside(COLUMN_loc,:) ;
    IND =   COOR_loc(:,1) == xMIN(:,1)  ;
    FACES_CN(IND,inodeE) =  1 ;
    IND =   COOR_loc(:,1) == xMAX(:,1)  ;
    FACES_CN(IND,inodeE) =  3 ;
    IND =   COOR_loc(:,2) == xMIN(:,2)  ;
    FACES_CN(IND,inodeE) =  2 ;
    IND =   COOR_loc(:,2) == xMAX(:,2)  ;
    FACES_CN(IND,inodeE) =  4 ;
end

end





function   [FACES_CN  ] = ClassificationNodesVert(CNvert,COORvert,xMIN,xMAX) ;

% Classification of vertices nodes
% AC = 1, AD = 2, BC=3 , BC=4


FACES_CN = zeros(size(CNvert)) ;
TOL = 1e-5 ;
for inodeE = 1:size(CNvert,2)  % loops over vertices
    COLUMN_loc = CNvert(:,inodeE) ;  % Column corresponding to node inodeE
    COOR_loc = COORvert(COLUMN_loc,:) ; % Coordinates of this column
%     IND_A =  COOR_loc(:,1) == xMIN(:,1) ;
%     IND_B =  COOR_loc(:,1) == xMAX(:,1) ;
%     IND_C =  COOR_loc(:,2) == xMIN(:,2) ;
%     IND_D =  COOR_loc(:,2) == xMAX(:,2) ;
    
    % AMENDED ON MAY 8th 2024
    IND_A =  abs(COOR_loc(:,1) - xMIN(:,1)) < TOL; %  COOR_loc(:,1) - xMIN(:,1) ;
   
      IND_B =  abs(COOR_loc(:,1) - xMAX(:,1)) < TOL;
          IND_C = abs(COOR_loc(:,2) - xMIN(:,2)) < TOL;  % COOR_loc(:,2) == xMIN(:,2) ;
       IND_D =  abs(COOR_loc(:,2) - xMAX(:,2)) < TOL;  %COOR_loc(:,2) == xMAX(:,2) ;
   
   if all(IND_A)+all(IND_B)+all(IND_C)+all(IND_D) ~=2 
       error('Correspondence not found...Set a lower tolerance ')
   end
    
    % FACES_CN(IND,inodeE) =  4 ;
    %%%%% NODE 4-1
    IND = IND_A & IND_D ;
    FACES_CN(IND,inodeE) = 1 ;
    %%%%% NODE AC
    IND = IND_A & IND_C ;
    FACES_CN(IND,inodeE) = 2 ;
    %%%%% NODE BC
    IND = IND_B & IND_C ;
    FACES_CN(IND,inodeE) = 3 ;
    %%%%% NODE BD
    IND = IND_B & IND_D ;
    FACES_CN(IND,inodeE) = 4 ;
    
end

end

