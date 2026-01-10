function [COOR,CN,TypeElement,CONNECTb,TypeElementB,Materials,MATERIALNEW]= ...
    MeshGenerationRepeat3D(NAMEMSH,DATA,MATERIAL)
% Automatic generation of meshes by tiling copies (translation in one
% direction). See DomainDecom_SVD.m
% ----------------------------------------------------------------
%addpath('DATA_input')
% INPUTS
% ---------

%dbstop('11')
if nargin ==0
    load('tmp5.mat')
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
for iDIR = 1:nDIR  % Loop over number of directions
    
    DATA.MakeMeshByRepetition.DIRECTIONloc = DATA.MakeMeshByRepetition.DIRECTION{iDIR} ;
    nDOM = nDOMglo(iDIR) ;
    if iDIR ==1
        
        % Rerefence mesh= NAMEMSH
        
        % 1. Reading mesh
        % ----------------
        
         [COOR,CN,TypeElement,CONNECTb,TypeElementB,MaterialType]=...
             ReadMeshFile(NAMEMSH,'READ_MATERIAL_COLUMN',1)  ;
        
%         SLICE.NAME =NAMEMSH ;
%         DATA3D = GeometrySlice(SLICE) ; % New function, for reading mesh and node faces 
%         % 
%         COOR = DATA3D.COOR ;
%         CN = DATA3D.CN ; 
%         TypeElement = DATA3D.TypeElement   ; 
%         TypeElementB = DATA3D.TypeElementB ; 
%         MaterialType = DATA3D.MaterialType ;  
        
        
    else
        MaterialType =   Materials ;
        MATERIAL = MATERIALNEW ;
    end
    
    
    
    
    MATERIALNEW.PLY = MATERIAL.PLY ;
    
    
    nmat = length(unique(MaterialType)) ;
    
    if nmat ~=length(MATERIAL.PLY)
        error('The number of materials should coincide with the number of materials defined in the input file')
    end
    
    if nDOM== 1
        MATERIALNEW =    MATERIAL   ;
        Materials =  MaterialType   ;
    else
        % FAces F1 and F2
        DATA.CalculateMasterSlaves  = 1 ;
        [NODESfaces,NODEREF,f1NOD, f2NOD] =  PointPlanesRBODY(COOR,CN,DATA) ;
        % ---------------------------------------
        % We distinguish between connectivities corresponding to face f1NOD, to
        % face f2NOD, and the remaining connectivities
        [CNbF1NOD setBelemF1NOD]= ElemBnd(CONNECTb,f1NOD);
        [CNbF2NOD setBelemF2NOD]= ElemBnd(CONNECTb,f2NOD);
        CONNECTbF1 = CONNECTb(setBelemF1NOD,:) ;
        CONNECTbF2 = CONNECTb(setBelemF2NOD,:) ;
        setBelemREST = setdiff(1:size(CONNECTb,1),[setBelemF1NOD setBelemF2NOD]) ;
        CONNECTbREST = CONNECTb(setBelemREST,:) ;
        % CONSTRUCTING THE MATRIX OF COORDINATES, CONNECTIVITIES AND MATERIAL
        % LIBRARY
        % -------------------------------------------------------------------
        COORglo = [COOR];
        CNglo = [CN]  ;
        CONNECTbRESTglo = CONNECTbREST ; % Boundary surfaces excluding face 1 and face 2
        Materials =MaterialType;
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
            % New boundary elements ---------------------
            % Outer and inner surfaces
            CONNECTbRESTnew = CONNECTbREST +(e-1)*size(COOR,1);
            % End surface (only at the end of the beam-like structure)
            if e==nDOM
                CONNECTbF2 = CONNECTbF2 +(e-1)*size(COOR,1);
            end
            % -------------------------------------------
            f1NEW = f1NOD + (e-1)*size(COOR,1) ;
            f2OLD = f2NOD + (e-2)*size(COOR,1) ;
            % Renumbering
            % We have to replace nodes (f1NEW) by f2NOD
            % in CNnew and CONNECTbRESTnew
            CNnewREN  = CNnew ;
            CONNECTbRESTnewREN = CONNECTbRESTnew ;
            for ifacen = 1:length(f1NOD)
                nodeLOC = f1NEW(ifacen);
                INDnodes = find(CNnew==nodeLOC) ;
                INDnodesb = find(CONNECTbRESTnew==nodeLOC) ;
                CNnewREN(INDnodes) = f2OLD(ifacen) ;
                CONNECTbRESTnewREN(INDnodesb) =f2OLD(ifacen) ;
            end
            CNglo = [CNglo ; CNnewREN] ;
            CONNECTbRESTglo = [CONNECTbRESTglo; CONNECTbRESTnewREN] ; % Notice that boundary connectivities at the interfaces are repeated
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
        CONNECTbGLO = [CONNECTbF1;CONNECTbF2;CONNECTbRESTglo] ;
        
        disp('REnumbering...')
        
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
        
    end
    
end

% else
%    CONNECTbGLO= RenumberConnectivities(CONNECTbGLO,NODES_bnd) ;
% end

