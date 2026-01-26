function  CreateNewJoint_old(NAMEMESH_SLICE,CNint,MaterialType1D,itype,DATAIN,SLICEglo,rotMATglo,coorINT)
% This function ...
% JAHO, 29-Dec-2017
if nargin == 0
    load('tmp.mat')
end

DATAIN = DefaultField(DATAIN,'WriteBatch',1) ;

% Detect 1D boundary nodes:  OUT --> BND.nodes, and BND.slices
[BND ] = DetectBoundaryNodesJoints(CNint,MaterialType1D,itype) ;
indLINESgid = 0 ;  TEXT = {} ;
nvolACUM = 0 ;
for inode = 1:length(BND.nodes)
    node = BND.nodes(inode) ;
    % Coordinate centroid
    coorCENT = coorINT(node,:);
    % Number of slice (1D element)
    slice = BND.slices(inode) ;
    % Type slice
    typeSlice = MaterialType1D(slice) ;
    % Rotation matrix
    % Numbering slice within its type
    INDX = find(MaterialType1D == typeSlice) ;
    sliceLOC = find(INDX == slice) ;
    dof = small2large(sliceLOC,length(coorCENT)) ;
    R = rotMATglo{typeSlice}(:,dof ) ;
    % Determining which face is node "node", as well as the associated 3D
    % nodes, and the centroid in local coordinates
    NODES = CNint(slice,:) ;
    if node == NODES(1) ;
        %
        NODES_face = SLICEglo{typeSlice}.NODES_face1 ;
        CENTR = SLICEglo{typeSlice}.CENTRf1 ;
        CENTR_opp = SLICEglo{typeSlice}.CENTRf2 ;
        
    else
        NODES_face = SLICEglo{typeSlice}.NODES_face2 ;
        CENTR = SLICEglo{typeSlice}.CENTRf2 ;
        CENTR_opp = SLICEglo{typeSlice}.CENTRf1 ;
    end
    % Coordinates of 3D nodes
    COOR3D = SLICEglo{typeSlice}.COOR ;
    CNb3D = SLICEglo{typeSlice}.CNb ;
    % Associated boundary elements
    [CNbF1 setBelem] =  ElemBnd(CNb3D,NODES_face) ;
    % Writting GID's batch
    % Transformation of coordinates
    %
    TranslationVector = coorCENT-CENTR_opp' ;  % Translation vector (for GID's batch file)
    COORtransf = CoordinatesChange(COOR3D',CENTR,coorCENT',R) ;
    % Batch GID for writing the corresponding faces
    
    %     CREATE_SURFACES = 0 ;
    %     if CREATE_SURFACES == 1
    %         [ indLINESgid,TEXT] = BatchGID_FaceNodes(COORtransf',CNbF1,TEXT,indLINESgid) ;
    %     else
    % Write face points ---just to make sure conformity is met.
    if DATAIN.WriteBatch == 1
        [dummy1 NAMEMESH_JOINT    dummy2]= fileparts(NAMEMESH_SLICE{itype}) ;
        [TEXT] = BatchGID_FaceNodes_points(COORtransf',CNbF1,TEXT) ;
        
        % Insert GID geometries
        % ---------------------
        TEXT{end+1} = 'Mescape Files InsertGeom' ;
        TEXT{end+1} = [cd,'/',NAMEMESH_SLICE{typeSlice}(1:end-4),'.gid'] ;
        % Hom many volumes  conform this slice ?
        nmatLOC = length(unique( SLICEglo{typeSlice}.MaterialType)) ;
        TEXT{end+1} = ['Mescape Utilities Move Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin ',num2str(TranslationVector(1)),',',num2str(TranslationVector(2)),',',num2str(TranslationVector(3))] ;
        VOLUMES = nvolACUM + (1:nmatLOC) ;
        TEXT{end+1} = [num2str(VOLUMES)] ;
        nvolACUM = nvolACUM +nmatLOC ;
        % 1
        TEXT{end+1} =  'escape' ;
        TEXT{end+1} = 'Mescape' ;
        
    end
    
    
    
end

if DATAIN.WriteBatch == 1
    FILE_BATCH = ['MESH_slices/',NAMEMESH_JOINT,'_START','.bch'] ;
    %
    %     Mescape Files SaveAs
    % -alsoresults:1 --
    % geoversion:current
    % /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/aaaaa.gid
    % *****SAVE MARK /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/aaaaa
    % *****SAVE STATE
    TEXT{end+1} =   ['Mescape Files SaveAs '];
    TEXT{end+1} =   ['-alsoresults:1 --'];
    TEXT{end+1} =   ['geoversion:current'];
    TEXT{end+1} =   [cd,'/',NAMEMESH_SLICE{itype}(1:end-4)];
    TEXT{end+1} =   ['*****SAVE MARK ',cd,'/',NAMEMESH_SLICE{itype}(1:end-4)];
    
    fid =fopen(FILE_BATCH,'w');
    for i = 1:length(TEXT)
        fprintf(fid,[TEXT{i},'\n']);
    end
    fod =fclose(fid);
    
    % Run gid
    %
    % Mescape Files SaveAs
    % -alsoresults:1 --
    % geoversion:current
    % /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/aaaaa.gid
    % *****SAVE MARK /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/MESH_slices/aaaaa
    % *****SAVE STATE
    
    comm  = ['gid    -b ',FILE_BATCH,' '];
    %iux
    %comm{iux}
    unix(comm );
    
    
    error('Finish to draw the geometry in GID and then proceed')
    
end


end



