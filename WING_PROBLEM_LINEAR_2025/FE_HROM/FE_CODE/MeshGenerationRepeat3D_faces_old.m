function [COOR,CN,TypeElement,CONNECTb,TypeElementB,Materials,MATERIALNEW,NODESfaces]= ...
    MeshGenerationRepeat3D_faces_old(NAMEMSH,DATA,MATERIAL)
% Automatic generation of meshes by tiling copies (translation in one
% direction). See DomainDecom_SVD.m
% Copy of MeshGenerationRepeat3D, but with enhanced properties. (April 2018)
% JaHO
% ----------------------------------------------------------------
%addpath('DATA_input')
% INPUTS
% ---------

%dbstop('11')
if nargin ==0
    load('tmp4.mat')
    % DATA.TypeUnitCell = 'HEXAHEDRA' ;
end
%%%%

nDOMglo = DATA.MakeMeshByRepetition.nDOM ; % number of domains per direction
nDIR = length(nDOMglo) ; % Number of directions
if  ~iscell(DATA.MakeMeshByRepetition.DIRECTION)
    DATA.MakeMeshByRepetition.DIRECTION = { DATA.MakeMeshByRepetition.DIRECTION} ;
end


disp('REpeating cells ...')

%dbstop('27')





nDOM = nDOMglo;

% Rerefence mesh= NAMEMSH

% 1. Reading mesh
% ----------------
%
%          [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
%              ReadMeshFile(NAMEMSH,'READ_MATERIAL_COLUMN',1)  ;

SLICE.NAME =NAMEMSH ;
DATA3D = GeometrySlice(SLICE) ; % New function, for reading mesh and node faces
%
COOR = DATA3D.COOR ;
CN = DATA3D.CN ;
TypeElement = DATA3D.TypeElement   ;
TypeElementB = DATA3D.TypeElementB ;
MaterialType = DATA3D.MaterialType ;
CONNECTb = DATA3D.CNb ;

% else
%     MaterialType =   Materials ;
%     MATERIAL = MATERIALNEW ;
% end




MATERIALNEW.PLY = MATERIAL.PLY ;


nmat = length(unique(MaterialType)) ;

if nmat ~=length(MATERIAL.PLY)
    error('The number of materials should coincide with the number of materials defined in the input file')
end
NODESfaces = DATA3D.NODES_FACES ;  % All face nodes (labels from GID)

% New version. We only take into account connectivities of the
% faces specified by NODESfaces (by the user, in GID)
CONNECTb_faces = {} ;
for iface = 1:length(NODESfaces)
    [dummy setBelemLOC]= ElemBnd(CONNECTb,NODESfaces{iface}); % elements face "iface"
    % Connectivities faces f1 and f2
    CONNECTb_faces{iface} = CONNECTb(setBelemLOC,:) ;
end


if nDOM== 1
    MATERIALNEW =    MATERIAL   ;
    Materials =  MaterialType   ;
    CONNECTb = CONNECTb_faces ;
    for iface = 1:length(NODESfaces)
        CONNECTb{iface} = {CONNECTb_faces{iface}} ;
    end
else
    % FAces F1 and F2
    % ----------------
    % 
    % NEW VERSION
    % -----------
    f1NOD = DATA3D.NODES_FACES{1} ;  % Nodes face 1, domain 1
    f2NOD = DATA3D.NODES_FACES{2} ;  % Nodes face 2, domain 1
    
    
    % ---------------------------------------
    % We distinguish between connectivities corresponding to face f1
    % and
    % face f2, and the remaining connectivities
    
    
    CONNECTbF1 = CONNECTb_faces{1}  ; % Connectivity FACE 1
    CONNECTbF2 = CONNECTb_faces{2}  ; % Connectivity FACE 2
    
    CONNECTbREST = CONNECTb_faces(3:end) ;  % Remaining faces 
    
    %%%%5
    % CONSTRUCTING THE MATRIX OF COORDINATES, CONNECTIVITIES AND MATERIAL
    % LIBRARY
    % -------------------------------------------------------------------
    COORglo = [COOR];  % This is 
    CNglo = [CN]  ;
    %  CNdom = cell(1,nDOM) ;  % Domain-wise connectivities
    % CNdom{1} = CN ;
    CONNECTbRESTglo =  cell(nDOM,1) ;
    CONNECTbRESTglo{1} =  CONNECTbREST;
    
    Materials =MaterialType;
    translationVECTOR = COOR(f2NOD(1),:)-COOR(f1NOD(1),:);  % This should be changed for curved elements
    TRANSLATION =0 ;
    for e = 2:nDOM
        disp(['Cell =',num2str(e)])
        TRANSLATION = TRANSLATION +  translationVECTOR ;
        % New coordinates (of domain e)
        COORnew=  COOR + repmat(TRANSLATION,size(COOR,1),1) ;
        % New connectivity matrix. (of domain e)
        %-------------------------
        CNnew = CN +(e-1)*size(COOR,1);
        
        % New boundary elements ---------------------
        % Labeled surfaces
        for iface =1:length(CONNECTbREST)
            CONNECTbRESTglo{e}{iface} = CONNECTbREST{iface} +(e-1)*size(COOR,1);
        end
        % End surface (only at the end of the beam-like structure)
        if e==nDOM
            CONNECTbF2 = CONNECTbF2 +(e-1)*size(COOR,1);
        end
        % -------------------------------------------
        f1NEW = f1NOD + (e-1)*size(COOR,1) ;  % Numbering of new face 1
        f2OLD = f2NOD + (e-2)*size(COOR,1) ;  % Numbering of old face 2
        % Renumbering
        % We have to replace nodes (f1NEW) by f2NOD
        % in CNnew and CONNECTbRESTnew
        CNnewREN  = CNnew ;
        CONNECTbRESTnewREN = CONNECTbRESTglo{e} ;
        for ifacen = 1:length(f1NOD)
            nodeLOC = f1NEW(ifacen);
            % Replacing it in the connectivity matrix
            INDnodes = find(CNnew==nodeLOC) ;
            CNnewREN(INDnodes) = f2OLD(ifacen) ;
            % Replacing it in the face connectivity matrix
            for iface = 1:length(CONNECTbRESTnewREN)
                INDnodesb = find(CONNECTbRESTglo{e}{iface}==nodeLOC) ;
                CONNECTbRESTnewREN{iface}(INDnodesb) =f2OLD(ifacen) ;
            end
        end
        CNglo = [CNglo ; CNnewREN] ;
        %  CNdom{e} = CNnewREN;
        CONNECTbRESTglo{e} = [CONNECTbRESTnewREN] ; %
        NewMaterial = MaterialType +  (e-1)*nmat ;
        Materials =[Materials; NewMaterial]  ;
        
        for imat = 1:length(MATERIAL.PLY)
            matREF = (e-1)*nmat ;
            MATERIALNEW.PLY(matREF+imat) =  MATERIAL.PLY(imat) ;
        end
        
        COORglo = [COORglo; COORnew] ;
    end
    
    % Boundary CONNECTIVITIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Face 1, Face 2, and outer and inner surfaces
    CONNECTbGLO = cell(length(NODESfaces),1);
    CONNECTbGLO{1} = {CONNECTbF1};
    CONNECTbGLO{2} = {CONNECTbF2};
    for iface = 3:length(NODESfaces)
        CONNECTbGLO{iface} = cell(nDOM,1) ;
        for idom = 1:nDOM
            CONNECTbGLO{iface}{idom} =  CONNECTbRESTglo{idom}{iface-2} ;
        end
    end
    CONNECTb =  CONNECTbGLO ;
    
    disp('REnumbering...')
    
    NODES = unique(CNglo(:)) ;
    COOR = COORglo(NODES,:) ;
    NODES_new = 1:length(NODES) ;
    % Interior CN
    CN = RenumberConnectivities(CNglo,NODES_new) ;
    % Boundary CN
    for iface = 1:length(CONNECTbGLO)
        for idom = 1:length(CONNECTbGLO{iface})
            NODESbnd = unique(CONNECTbGLO{iface}{idom}) ;
            [~,NODES_bnd,~] = intersect(NODES,NODESbnd) ;
            
            CONNECTb{iface}{idom}= RenumberConnectivities(CONNECTbGLO{iface}{idom},NODES_bnd) ;
        end
        
    end
    
end


end
