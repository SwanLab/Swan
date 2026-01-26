function [DISP3D,DISP3D_lateral,stressGLO,STRESS_DATA,DATAIN,DATA_REFMESH,REACTIONS3D] ...
    = Displacement_stress_3D_JOINTslice(DATAROM,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH,DATARUN,rDEF,rRB,DATAadd)


if nargin ==0
    load('tmp2.mat')
elseif nargin == 9 
    DATAadd = [] ; 
end
DATAIN = DefaultField(DATAIN,'DOMAINS_POSTPROCESS_SELECT',[]) ;
DATAIN.DOMAINS_POSTPROCESS_SELECT = DefaultField(DATAIN.DOMAINS_POSTPROCESS_SELECT,'NUMBER',[]) ;
DATAIN.DOMAINS_POSTPROCESS_SELECT = DefaultField(DATAIN.DOMAINS_POSTPROCESS_SELECT,'VARIABLE','VONMISES') ;
DATAIN = DefaultField(DATAIN,'POST_PROCESS_LATERAL_SURFACES',[]) ;
DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'ACTIVE',0) ;
DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'COARSE_MESH',1) ;
%DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'NAME_MESH_COARSE',[]) ;

DATAIN = DefaultField(DATAIN,'PRINT_DISTRIBUTED_FORCES',0) ; % Print in GID external traction forces (over external surfaces)
% If this option is enabled, then no inner elements are plotted (only lateral surfaces)
% Likewise, the plot is made using the finer mesh
nDOM = size(qDEF,2) ;
if DATAIN.PRINT_DISTRIBUTED_FORCES == 1
    DATAIN.DOMAINS_POSTPROCESS = [] ;
    DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER = [] ;
    NO_STRESSES = 1;
else
    DATAIN = DefaultField(DATAIN,'DOMAINS_POSTPROCESS',[]) ;
    if isempty(DATAIN.DOMAINS_POSTPROCESS)  && isempty(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER)
        DATAIN.DOMAINS_POSTPROCESS = 1:nDOM ;
    end
    NO_STRESSES = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%555
DISP3D = [] ;
DISP3D_lateral = [] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DISP3D_all = DISP3D ;
%--------------
%% STRESSES
% ---------
[stressGLO,STRESS_DATA,SELECTED_DOMAINS] =...
    StressDom_SlicesJoints(DATAROM,MESH1D,DATAIN,qDEF,DATA_REFMESH,NO_STRESSES,DATARUN,DATAadd) ;

%% DISPLACEMENTS
% ------------------
DISP3D = ComputeDisplacements(DATAIN,DATAROM,MESH1D,qRB,DATA_REFMESH,qDEF,SELECTED_DOMAINS) ;



% TO BE IMPLEMENTED
% ---------------------
DATAIN = DefaultField(DATAIN,'PRINT_REACTIONS_INTERFACE_SELF_EQUILIBRATED',0) ;

if DATAIN.PRINT_REACTIONS_INTERFACE_SELF_EQUILIBRATED == 1
    REACTIONS3D = ComputeReactions(DATAIN,DATAROM,MESH1D,rRB,DATA_REFMESH,rDEF,SELECTED_DOMAINS) ;
else
    REACTIONS3D = [] ;
end

% DISPLACEMENT LATERAL SURFACES
% --------------------------------
[DISP3D,DATAIN ]= DisplacementLateralSurfacesJ(DATAIN,qDEF,DATA_REFMESH,DISP3D,...
    SELECTED_DOMAINS,qRB,DATAROM,MESH1D) ;
%b



%----------------------------------------------------------------------------------

DATAIN.SELECTED_DOMAINS = SELECTED_DOMAINS ;

end

% --------------------------------



function  DISP3D = ComputeDisplacements(DATAIN,DATAROM,MESH1D,qRB,DATA_REFMESH,qDEF,SELECTED_DOMAINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATAIN = DefaultField(DATAIN,'POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL',0) ;
nentities =length(DATAROM) ;
DISP3D = cell(1,nentities) ;
for ientity = 1:length(DATAROM)
    ELEMS = find(MESH1D.MaterialType == ientity) ;
    ELEMS = ELEMS(SELECTED_DOMAINS{ientity}) ;
    
    % ELEMS =
    BasisUdef = DATAROM{ientity}.BasisUdef ;
    BasisUrb = DATA_REFMESH{ientity}.BasisUrb;
    nDOM = length(ELEMS) ;
    DISP3D{ientity} =zeros(size(BasisUdef,1),nDOM) ;
    if DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 0
        DISP3D{ientity} = BasisUdef*qDEF(:,ELEMS) + BasisUrb*qRB(:,ELEMS) ;
    elseif DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 1
        DISP3D{ientity}=   BasisUrb*qRB(:,ELEMS);
    elseif DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 2
        DISP3D{ientity} =   BasisUdef*qDEF(:,ELEMS)  ;
    else
        error('OPtion not implemented')
    end
    
    %%% ROTATION OF DISPLACEMENTS 
    % ---------------------------
    DATA_REFMESH{ientity} = DefaultField( DATA_REFMESH{ientity},'RotationMatrixFace',[]) ;
    if ~isempty(DATA_REFMESH{ientity}.RotationMatrixFace)
        Rotation = DATA_REFMESH{ientity}.RotationMatrixFace ;
        if ~isempty(Rotation{2})
            %%% Rotation of displacements (from domain coordinates to global coordinates )
            % -----------------------------------------------------------------------------
            % We begin by implementing the non-vectorized version
            ndim = size(MESH1D.ROTATIONS,1) ;
            for ielemLOC = 1:length(ELEMS)
                ielem = ELEMS(ielemLOC) ;
                finIND =  ndim*ielem ;
                iniIND = ndim*ielem - ndim+1 ;
                R = MESH1D.ROTATIONS(:,iniIND:finIND) ;
                dispLOC =  reshape(DISP3D{ientity}(:,ielemLOC),ndim,[]) ;
                dispLOC = R*dispLOC ;
                DISP3D{ientity}(:,ielemLOC)  =dispLOC(:) ;
            end
        end
    end
    
    
     
    
    
end


end




function  REACTIONS3D = ComputeReactions(DATAIN,DATAROM,MESH1D,rRB,DATA_REFMESH,rDEF,SELECTED_DOMAINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nentities =length(DATAROM) ;
REACTIONS3D = cell(1,nentities) ;
for ientity = 1:length(DATAROM)
    ELEMS = find(MESH1D.MaterialType == ientity) ;
    ELEMS = ELEMS(SELECTED_DOMAINS{ientity}) ;
    
    % ELEMS =
    BasisRdef = DATAROM{ientity}.BasisRdef ;
    %  BasisRrb = DATA_REFMESH{ientity}.BasisRrb;
    nDOM = length(ELEMS) ;
    REACTIONS3D{ientity} =zeros(size(BasisRdef,1),nDOM) ;
    %  if DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 0
    REACTIONS3D{ientity} = BasisRdef*rDEF(:,ELEMS) ;% + BasisRrb*rRB(:,ELEMS) ;
    REACTIONS3D{ientity} =   REACTIONS3D{ientity}(:) ;
    % elseif DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 1
    %   DISP3D{ientity}=   BasisUrb*qRB(:,ELEMS);
    %elseif DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 2
    %  DISP3D{ientity} =   BasisUdef*qDEF(:,ELEMS)  ;
    %else
    %   error('OPtion not implemented')
    %end
end


end

