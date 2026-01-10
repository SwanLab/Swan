function [DISP3D,DATAIN ]= DisplacementLateralSurfacesRVE(DATAIN,qDEF,DATA_REFMESH,DISP3D,SELECTED_DOMAINS,...
    qRB,DATAROM,MESH1D)

if nargin == 0
    load('tmp.mat')
end






% DISPLACEMENTS LATERAL SURFACES
% -------------------------------
ndim = 3;
CONNECTb = [] ;
DISP3D_lateral = 0 ;


DATAIN = DefaultField(DATAIN,'PRINT_DISTRIBUTED_FORCES',0) ;

if DATAIN.PRINT_DISTRIBUTED_FORCES == 1
    DATAIN.POST_PROCESS_LATERAL_SURFACES.COARSE_MESH  = 0 ;
end

if  length(DATAIN.DOMAINS_POSTPROCESS)  ~=size(qDEF,2) && DATAIN.POST_PROCESS_LATERAL_SURFACES.ACTIVE == 1
    
    % Connectivities of COARSE MESH
    % -------------------------------
    if DATAIN.POST_PROCESS_LATERAL_SURFACES.COARSE_MESH == 1
        for ientity = 1:length(DATA_REFMESH)
            if isempty(DATA_REFMESH{ientity}.CONNECTb_coarse)
            else
                DATA_REFMESH{ientity}.CONNECTb =  DATA_REFMESH{ientity}.CONNECTb_coarse  ;
                DATA_REFMESH{ientity}.TypeElementB = DATA_REFMESH{ientity}.TypeElementB_coarse  ;
            end
        end
    else
        % COONECTb is determined from DATA_REFMESH
    end
    % -------------------------------------------------
    DISP3D_lateral = 1;
    nfacesINTER = length(DATA_REFMESH{1}.CENTRf);
    
    for ientity = 1:length(DATA_REFMESH)
        %DISP3D{ientity} = BasisUdef*qDEF(:,ELEMS) + BasisUrb*qRB(:,ELEMS) ;
        ELEMS = find(MESH1D.MaterialType == ientity) ;
        BasisUdef = DATAROM{ientity}.BasisUdef ;
        BasisUrb = DATA_REFMESH{ientity}.BasisUrb;
        % Structural entity "ientity"
        ilateralSURF = (nfacesINTER+1):length( DATA_REFMESH{ientity}.CONNECTb) ;
        CONNECTb = DATA_REFMESH{ientity}.CONNECTb(ilateralSURF) ; % Connectiviites lateral surfaces
        NODES = cell2mat(CONNECTb') ;  % Nodes lateral surfaces
        NODES = unique(NODES(:)) ;
        DOFS = small2large(NODES,ndim) ; % DOFs lateral surface
        % New array
        DISP3D_new{ientity} =cell(1,length(ELEMS)) ;
        SELECTED = zeros(1,length(ELEMS)) ;
        SELECTED(SELECTED_DOMAINS{ientity}) = 1; % Selected domains
        icum = 0 ;
        for idom = 1:length(ELEMS)
            ISEL = ELEMS(idom) ; 
            if SELECTED(idom) == 1
                % Nothing is done
                icum = icum  + 1;
                DISP3D_new{ientity}{idom} = DISP3D{ientity}(:,icum) ;
            else
                % Only some DOFs are selected (those corresponding to lateral surfaces )
                DISPLOC =  BasisUdef(DOFS,:)*qDEF(:,ISEL) + BasisUrb(DOFS,:)*qRB(:,ISEL) ;
                DISP3D_new{ientity}{idom} = DISPLOC ;
            end
        end
        % Then everything is converted into a MATRIX
        DISP3D{ientity} = cell2mat(DISP3D_new{ientity}(:)) ;
    end
    
else
    
    DISP3D_lateral =0;
    for ientity = 1:length(DISP3D)
        %      DISP3D{ientity} = DISP3D{ientity}(:,SELECTED_DOMAINS{ientity}) ;
        DISP3D{ientity} = DISP3D{ientity}(:) ;
    end
    
end

DATAIN.DISP3D_lateral   = DISP3D_lateral ;

end