function [COOR,CN,TypeElement,CONNECTb,TypeElementB,Materials]= MeshGenerationRepeat(NAMEMSH,nDOM,DATA)
% Automatic generation of meshes by tiling copies (translation in one direction)
% ----------------------------------------------------------------
%addpath('DATA_input')
% INPUTS
% ---------


if nargin ==0
    load('tmp2.mat')
    DATA.TypeUnitCell = 'HEXAG_·D_SQUARE' ; 
end
%%%%

disp('REpeating cells ...')

% 1. Reading mesh
% ----------------
%dbstop('20')
DATA = DefaultField(DATA,'TypeUnitCell','HEXAG_2D_SQUARE') ;

[COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
    ReadMeshFile(NAMEMSH,'READ_MATERIAL_COLUMN',1)  ;
MaterialType = ones(size(CN,1),1) ;

COOR = COOR(:,1:2);


nmat = length(unique(MaterialType)) ;
% FAces F1 and F2
[NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COOR,CN,DATA) ;

% Remove f2NOD from this set
% NODESb = setdiff(unique(CONNECTb(:)),f2NOD) ;
% [CNb setBelem]= ElemBnd(CONNECTb,NODESb);
% CONNETbGLO = [CNb] ;



CONNECTbGLO = [CONNECTb] ;


% CONSTRUCTING THE MATRIX OF COORDINATES, CONNECTIVITIES AND MATERIAL
% LIBRARY
% -------------------------------------------------------------------
COORglo = [COOR];
CNglo = [CN]  ;
Materials = ones(size(CN,1),1) ;
translationVECTOR = COOR(f2NOD(1),:)-COOR(f1NOD(1),:);
TRANSLATION =0 ;


for e = 2:nDOM
    disp(['Cell =',num2str(e)])
    TRANSLATION = TRANSLATION +  translationVECTOR ;
    % New coordinates
    COORnew=  COOR + repmat(TRANSLATION,size(COOR,1),1) ;
    % New connectivity matrix.
    %-------------------------
    CNnew = CN +(e-1)*size(COOR,1);
    CONNECTbNEW = CONNECTb +(e-1)*size(COOR,1);
    f1NEW = f1NOD + (e-1)*size(COOR,1) ;
    f2OLD = f2NOD + (e-2)*size(COOR,1) ;
    % Now we have to replace nodes (f1NEW) by f2NOD
    % in CNnew and CONNECTbNEW
    CNnewREN  = CNnew ;
    CONNECTbNEWren = CONNECTbNEW ;
    for ifacen = 1:length(f1NOD)
        nodeLOC = f1NEW(ifacen);
        INDnodes = find(CNnew==nodeLOC) ;
        INDnodesb = find(CONNECTbNEW==nodeLOC) ;
        
        CNnewREN(INDnodes) = f2OLD(ifacen) ;
        CONNECTbNEWren(INDnodesb) =f2OLD(ifacen) ;
    end
    CNglo = [CNglo ; CNnewREN] ;
    CONNECTbGLO = [CONNECTbGLO; CONNECTbNEWren] ;
    NewMaterial = MaterialType +  (e-1)*nmat ;
    Materials =[Materials; NewMaterial]  ;
    COORglo = [COORglo; COORnew] ;
end

% Boundary nodes
% --------------

disp('Boundary NODES...')
[NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COORglo,CNglo,DATA) ;

NODESB = [] ;
for i = 1:length(NODESfaces)
    NODESB = [NODESB;NODESfaces{i}(:)] ;
end

CONNECTbGLO = ElemBnd(CONNECTbGLO,NODESB) ;




 METHOD_RENUMBER = 1;

disp('REnumbering...')


if METHOD_RENUMBER == 1
    NODES = unique(CNglo(:)) ;
    COOR = COORglo(NODES,:) ;
    NODES_new = 1:length(NODES) ;
    % Interior CN
    CN = RenumberConnectivities(CNglo,NODES_new) ;
    % Boundary CN 
    NODESbnd = unique(CONNECTbGLO(:)) ;
    [~,NODES_bnd,~] = intersect(NODES,NODESbnd) ; 
     CONNECTb = RenumberConnectivities(CONNECTbGLO,NODES_bnd) ;
    
    
else
    
    
    %% NOW WE REMOVE THOSE NODES THAT ARE NOT EMPLOYED
    % -----------------------------------------------
    
    
    
    nnode = length(NODES);
    CN = zeros(size(CNglo)) ;
    CONNECTb = zeros(size(CONNECTbGLO)) ;
    disp('Renumbering connectivities')
    for inode = 1:length(NODES)
        nodeLOC = NODES(inode) ;
        INDnodes = find(CNglo==nodeLOC) ;
        CN(INDnodes) = inode ;
        INDnodes = find(CONNECTbGLO==nodeLOC) ;
        if ~isempty(INDnodes)
            CONNECTb(INDnodes) = inode ;
            
        end
    end
    
end




