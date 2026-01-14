function  CONNECTslices = CreateNewJoint(MESH1D,MESH3D,itypeJOINT,WRITE_BATCH_FILE) ;
% See BeamROM;.pdf
% JAHO, 29-Dec-2017
%
if nargin == 0
    load('tmp1.mat')
end

% ---------------------------------------------------------------------
% Skeleton 1D mesh file may be different from the one defined for the
% problem under consideration
% Skeleton 1D mesh is the 1D mesh defined by the user when the joint was
% created. The centroids of the 3D mesh should coincide with the nodes of
% the 1D mesh
% ---------------------------------------------------------------------
JOINT = MESH3D.JOINTS(itypeJOINT)  ;
JOINT = DefaultField(JOINT,'MESH1D_SKELETON',[]) ;
JOINT.MESH1D_SKELETON = DefaultField(JOINT.MESH1D_SKELETON,'STRUCTURE',MESH1D) ;

if strcmp(JOINT.MESH1D_SKELETON.STRUCTURE.NAME,MESH1D.NAME)
    %  disp('No need to recalculate')
else
    disp('Re-computing structure variable MESH1D')
    MESH3D_old = MESH3D ;
%   run(JOINT.MESH1D_SKELETON.STRUCTURE.NAME) ;
    MESH1D_new = MESH1D ; 
    MESH1D_new.NAME = JOINT.MESH1D_SKELETON.STRUCTURE.NAME ; 
    %DATAIN.PlotNormalsX = 0 ;
    DATAIN = [] ; 
    MESH1D = Geometry1Dstructure(MESH1D_new,MESH3D_old,DATAIN) ;
    MESH3D = MESH3D_old ;
end


% First we check which portions of the 1D mesh has assigned this type of
% joint -----------------------------------------------------------------
ientity = 1 ;
END_NODES = [] ;
SELECTED_JOINT = 1 ;  % We assume that the reference joint will be the first one
while ientity <= length(MESH1D.PROP)
    switch MESH1D.PROP(ientity).TYPE
        case 'JOINT'
            if MESH1D.PROP(ientity).INDEX_ENTITY ==itypeJOINT
                
                END_NODES = MESH1D.INFOLINES.END_NODES{ientity}
                if ~isempty(END_NODES)
                    END_NODES = END_NODES{SELECTED_JOINT} ;
                end
                END_ELEMENTS = MESH1D.INFOLINES.END_ELEMENTS{ientity}
                if ~isempty(END_ELEMENTS)
                    END_ELEMENTS = END_ELEMENTS{SELECTED_JOINT} ;
                end
                if isempty(END_NODES)
                    % Two-faces joint
                    % ---------------
                    END_NODES = MESH1D.INFOLINES.NODES{ientity}{SELECTED_JOINT}([1,end]) ;
                    END_ELEMENTS = MESH1D.INFOLINES.ELEMENTS{ientity}{SELECTED_JOINT}([1,end]) ;
                end
                break
            end
    end
    ientity = ientity + 1;
end
if isempty(END_NODES)
    error('MESH1D.PROP(ientity).INDEX_ENTITY does not coincide with any existing joint')
end

% Once we know the "end" nodes of the studied joint, we have to check which
% are the neighboring slices, as well as the area of the contacting surface
TEXT = {} ;
nvolACUM = 0 ;


% Coordinates of the 3D nodes of the slices connected with joints whose end
% nodes and end elements are END_NODES and END_ELEMENTS
% OUTPUT
CONNECTslices = ConnectingSlices_JOINT(END_NODES,END_ELEMENTS,MESH1D,MESH3D) ;

FOLDER_B = fileparts(MESH3D.JOINTS(itypeJOINT).NAME) ;



if WRITE_BATCH_FILE == 1
    
    
    for  inodeLOC = 1:length(END_NODES)
        
        if ~isempty(CONNECTslices.COOR{inodeLOC})
            
            COORtransf = CONNECTslices.COOR{inodeLOC} ;
            TranslationVector = CONNECTslices.TRANSLATION{inodeLOC} ;
            [dummy1 NAMEMESH_JOINT    dummy2]= fileparts(MESH3D.JOINTS(itypeJOINT).NAME) ;
            for inode = 1:size(COORtransf,2)
                TEXT{end+1} = 'Mescape Geometry Create Point' ;
                COORini = COORtransf(:,inode) ;
                TEXT{end+1} = [num2str(COORini(1)),',',num2str(COORini(2)),',',num2str(COORini(3))] ;
            end
            TEXT{end+1} = 'escape' ;
            % Insert GID geometries
            % ---------------------
            TEXT{end+1} = 'Mescape Files InsertGeom' ;
            indexSLICE = CONNECTslices.indexSLICE_GLO(inodeLOC) ;
            
            MESH3D.SLICES(indexSLICE) = DefaultField(MESH3D.SLICES(indexSLICE),'NAME_GIDPROJECT',MESH3D.SLICES(indexSLICE).NAME) ;
            TEXT{end+1} = [cd,filesep,MESH3D.SLICES(indexSLICE).NAME_GIDPROJECT,'.gid'] ;
            % Hom many volumes  conform this slice ?
            nmatLOC = length(unique( MESH3D.SLICES(indexSLICE).DATA3D.MaterialType)) ;
            TEXT{end+1} = ['Mescape Utilities Move Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin ',num2str(TranslationVector(1)),',',num2str(TranslationVector(2)),',',num2str(TranslationVector(3))] ;
            VOLUMES = nvolACUM + (1:nmatLOC) ;
            TEXT{end+1} = [num2str(VOLUMES)] ;
            nvolACUM = nvolACUM +nmatLOC ;
            % 1
            TEXT{end+1} =  'escape' ;
            TEXT{end+1} = 'Mescape' ;
        end
        
    end
    
    
    
    FILE_BATCH = [FOLDER_B,filesep,NAMEMESH_JOINT,'_START','.bch'] ;
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
    TEXT{end+1} =   [cd,filesep,MESH3D.JOINTS(itypeJOINT).NAME];
    TEXT{end+1} =   ['*****SAVE MARK ',cd,filesep,MESH3D.JOINTS(itypeJOINT).NAME];
    
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



