function [COOR,CN,TypeElement,CONNECTb,TypeElementB,Materials]= MeshGenerationRepeat3D(NAMEMSH,nDOM,DATA)
% Automatic generation of meshes by tiling copies (translation in one
% direction). See DomainDecom_SVD.m
% ----------------------------------------------------------------
%addpath('DATA_input')
% INPUTS
% ---------

%dbstop('10')
if nargin ==0
    load('tmp2.mat')
    DATA.TypeUnitCell = 'HEXAHEDRA' ;
end
%%%%

disp('REpeating cells ...')

% 1. Reading mesh
% ----------------
%dbstop('17')

[COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
    ReadMeshFile(NAMEMSH,'READ_MATERIAL_COLUMN',1)  ;
MaterialType = ones(size(CN,1),1) ;




nmat = length(unique(MaterialType)) ;
% FAces F1 and F2
DATA.CalculateMasterSlaves  = 1 ;
[NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COOR,CN,DATA) ;

ELIMINATE_INTERFACE_CONNECTIVITIES = 0; 

if ELIMINATE_INTERFACE_CONNECTIVITIES == 1
% ---------------------------------------
% We distinguish between connectivities corresponding to face f1NOD, to
% face f2NOD, and the remaining connectivities
 [CNbF1NOD setBelemF1NOD]= ElemBnd(CONNECTb,f1NOD);
 [CNbF2NOD setBelemF2NOD]= ElemBnd(CONNECTb,f2NOD);
 CONNECTbF1 = CONNECTb(setBelemF1NOD,:) ; 
 CONNECTbF2 = CONNECTb(setBelemF2NOD,:) ;
 setBelemREST = setdiff(1:size(CONNECTb,1),[setBelemF1NOD setBelemF2NOD]) ; 
 CONNECTbREST = CONNECTb(setBelemREST,:) ;


end 



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
    % New coordinates (of domain e)
    COORnew=  COOR + repmat(TRANSLATION,size(COOR,1),1) ;
    % New connectivity matrix. (of domain e)
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
    CONNECTbGLO = [CONNECTbGLO; CONNECTbNEWren] ; % Notice that boundary connectivities at the interfaces are repeated
    NewMaterial = MaterialType +  (e-1)*nmat ;
    Materials =[Materials; NewMaterial]  ;
    COORglo = [COORglo; COORnew] ;
end


% Change, 19-th-July-2017. Connectivities should be computed after
% renumbering
%
% Boundary nodes
% --------------
% 
% BNODES_AFTER = 0;
% 
% if  BNODES_AFTER == 0
    
    disp('Boundary NODES...')
    DATA.CalculateMasterSlaves  = 0 ;
    [NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COORglo,CNglo,DATA) ;
    
    NODESB = [] ;
    for i = 1:length(NODESfaces)
        NODESB = [NODESB;NODESfaces{i}(:)] ;
    end
    
    CONNECTbGLO = ElemBnd(CONNECTbGLO,NODESB) ;
    
% end


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
  %  if  BNODES_AFTER == 0
        CONNECTb= RenumberConnectivities(CONNECTbGLO,NODES_bnd) ;
        
   % else
    %    CONNECTbGLO= RenumberConnectivities(CONNECTbGLO,NODES_bnd) ;
   % end
    
    
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
%
%
%
% 
% if  BNODES_AFTER == 1
%     disp('Boundary NODES...')
%     DATA.CalculateMasterSlaves  = 0 ;
%     [NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COOR,CN,DATA) ;
%     
%     NODESB = [] ;
%     for i = 1:length(NODESfaces)
%         NODESB = [NODESB;NODESfaces{i}(:)] ;
%     end
%     
%     CONNECTb = ElemBnd(CONNECTbGLO,NODESB) ;
% end
% %
% 
% 
