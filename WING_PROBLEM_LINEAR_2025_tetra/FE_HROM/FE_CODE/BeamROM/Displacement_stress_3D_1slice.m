function [DISP3D,DISP3D_lateral,stressGLO,STRESS_DATA,DATAIN,DATA_REFMESH] ...
    = Displacement_stress_3D_1slice(DATAROM,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH)


if nargin ==0
    load('tmp.mat')
end


DATAIN = DefaultField(DATAIN,'DOMAINS_POSTPROCESS_SELECT',[]) ;
DATAIN.DOMAINS_POSTPROCESS_SELECT = DefaultField(DATAIN.DOMAINS_POSTPROCESS_SELECT,'NUMBER',[]) ;
DATAIN.DOMAINS_POSTPROCESS_SELECT = DefaultField(DATAIN.DOMAINS_POSTPROCESS_SELECT,'VARIABLE','VONMISES') ;
DATAIN = DefaultField(DATAIN,'POST_PROCESS_LATERAL_SURFACES',[]) ;
DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'ACTIVE',1) ;
DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'COARSE_MESH',1) ;
DATAIN.POST_PROCESS_LATERAL_SURFACES = DefaultField(DATAIN.POST_PROCESS_LATERAL_SURFACES,'NAME_MESH_COARSE',[]) ;

DATAIN = DefaultField(DATAIN,'PRINT_DISTRIBUTED_FORCES',0) ; % Print in GID external traction forces (over external surfaces)
% If this option is enabled, then no inner elements are plotted (only lateral surfaces)
% Likewise, the plot is made using the finer mesh
nDOM = size(qDEF,2) ;
if DATAIN.PRINT_DISTRIBUTED_FORCES == 1
    DATAIN.DOMAINS_POSTPROCESS = [] ;
    DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER = [] ;
    NO_STRESSES = 1;
else
    DATAIN = DefaultField(DATAIN,'DOMAINS_POSTPROCESS',1:nDOM) ;
    if isempty(DATAIN.DOMAINS_POSTPROCESS)  && isempty(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER)
        DATAIN.DOMAINS_POSTPROCESS = 1:nDOM ;
    end
    NO_STRESSES = 0;
end









DISP3D = [] ;
DISP3D_lateral = [] ;


% DOMAINS_POSTPROCESS = DATAIN.DOMAINS_POSTPROCESS ;
%
nDOM = size(qDEF,2) ;
% qDEF =  (qDEF(:,DOMAINS_POSTPROCESS)) ;
% qRB =  (qRB(:,DOMAINS_POSTPROCESS)) ;

BasisUdef = DATAROM.BasisUdef ;
BasisUrb = DATA_REFMESH.BasisUrb;
%DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL  = 1 ; % Print rigid body. == 2 Print deformational

DATAIN = DefaultField(DATAIN,'POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL',0) ;

if DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 0
    DISP3D = BasisUdef*qDEF + BasisUrb*qRB ;
elseif DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 1
    DISP3D =   BasisUrb*qRB ;
elseif DATAIN.POSTPROCESS_SEPARATE_DISPLACEMENTS_RB_DEFORMATIONAL == 2
    DISP3D =   BasisUdef*qDEF  ;
else
    error('OPtion not implemented')
end


DISP3D_all = DISP3D ;

% DISPLACEMENTS LATERAL SURFACES
% -------------------------------


%--------------
%% STRESSES
% ---------
disp('Retrieving OFFLINE data')
tic
load(DATAIN.NAME_WS_MODES,'CgloDOM','Wdom','Bdom')
toc
disp('DONE')
%CdomW = DATA_REFMESH.CgloDOM ;  %
%Wdom = DATA_REFMESH.Wdom ;
%Bdom = DATA_REFMESH.Bst ;

Bdom = Bdom*BasisUdef ;
CgloDOM = CgloDOM*Bdom ;

% STRAIN = Bdom*DISP3D ;
% STRESS = CdomW*STRAIN ;
nstrain = 6 ;
% ncomp = size(STRESS,1) ;

DATAIN = DefaultField(DATAIN,'PRINT_AVERAGE_STRESSES_ON_ELEMENTS',1) ;

if DATAIN.PRINT_AVERAGE_STRESSES_ON_ELEMENTS == 0
    error(['Option not implemented'])
end

ngaus = size(DATA_REFMESH.posgp,2) ;

Wst = repmat(Wdom',nstrain,1) ;
Wst = Wst(:) ;
nelem =  length(Wdom)/ngaus ;

stressGLO = zeros(nelem*nstrain,size(DISP3D,2)) ;
stressVONMISES = zeros(nelem,size(DISP3D,2)) ;
MAXstressVONMISES = zeros(1,size(DISP3D,2))  ;

DATAIN = DefaultField(DATAIN,'DO_NOT_COMPUTE_STRESSES',0) ;


if DATAIN.DO_NOT_COMPUTE_STRESSES == 0
    tic
    disp('Computing average stress on each FE element')
    
    for idom = 1:size(DISP3D,2)
        stressDOM_e = CgloDOM*qDEF(:,idom) ;  %
        for istrain =1:nstrain
            stressDOM_e(istrain:nstrain:end) = stressDOM_e(istrain:nstrain:end)./Wdom ;   % Cglo already includes WEIGHTS
        end
        [stressDOM_e] = AverageStressOnElements(stressDOM_e,Wst,nelem,nstrain,ngaus) ;
        
        % Von Mises
        stressDOM_ELEM = reshape(stressDOM_e,nstrain,[]) ;
        [ stressVONMISES_e ] =  VonMises_Stress(stressDOM_ELEM) ;
        MAXstressVONMISES(idom) = max(stressVONMISES_e) ;
        stressGLO(:,idom) = stressDOM_e ;
        stressVONMISES(:,idom) = stressVONMISES_e' ;
    end
    disp('...Done')
    toc
else
    stressVONMISES=[] ;
    stressGLO = [] ;
    MAXstressVONMISES = [] ;
    
end




%
% DATAIN.DOMAINS_POSTPROCESS_SELECT = DefaultField(DATAIN.DOMAINS_POSTPROCESS_SELECT,'NUMBER',[]) ;
% DATAIN.DOMAINS_POSTPROCESS_SELECT = DefaultField(DATAIN.DOMAINS_POSTPROCESS_SELECT,'VARIABLE','VONMISES') ;

if ~isempty(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER) & DATAIN.DO_NOT_COMPUTE_STRESSES ==0
    NUMBERprint = min(DATAIN.DOMAINS_POSTPROCESS_SELECT.NUMBER,size(stressGLO,2)) ;
    switch DATAIN.DOMAINS_POSTPROCESS_SELECT.VARIABLE
        case 'VONMISES'
            [~,INDEXES] =  sort(MAXstressVONMISES,'descend') ;
            SELECTED_DOMAINS = INDEXES(1:NUMBERprint) ;
            
            SELECTED_DOMAINS = [1,SELECTED_DOMAINS,length(MAXstressVONMISES)] ; % First and last domains are included
            SELECTED_DOMAINS = unique(SELECTED_DOMAINS) ;
            stressVONMISES =  stressVONMISES(:,SELECTED_DOMAINS) ;
            stressGLO =  stressGLO(:,SELECTED_DOMAINS) ;
            %      DISP3D = DISP3D(:,SELECTED_DOMAINS) ;
            DATAIN.DOMAINS_POSTPROCESS = SELECTED_DOMAINS ;
             
        otherwise
            error('Option not implemented')
    end
    
else
    if ~isempty(DATAIN.DOMAINS_POSTPROCESS) && NO_STRESSES == 0 && DATAIN.DO_NOT_COMPUTE_STRESSES ==0
        stressVONMISES =  stressVONMISES(:,DATAIN.DOMAINS_POSTPROCESS) ;
        stressGLO =  stressGLO(:,DATAIN.DOMAINS_POSTPROCESS) ;
        SELECTED_DOMAINS = DATAIN.DOMAINS_POSTPROCESS ;
    elseif  isempty(DATAIN.DOMAINS_POSTPROCESS) && NO_STRESSES == 1
        stressVONMISES =  [] ;
        stressGLO =  [] ;
        SELECTED_DOMAINS = [];
    end
    
    
end




%%% GID PRINTING
if ~isempty(stressGLO)
    
  %  indGID = [1 2 3 6 4 5] ;
    
     if size(MESH2D.COOR,2) == 3
            indGID = [1 2 3 6 4 5] ;
        else
            indGID = [1 2 3 4] ;
        end
    
    stressGLO_old = stressGLO;
    for iglo = 1:length(indGID)
        stressGLO(iglo:nstrain:end)  =  stressGLO_old(indGID(iglo):nstrain:end) ;
    end
    
    ngaus = 1;
    stressGLO =  reshape(stressGLO,ngaus*nstrain,[]) ;
    STRESS_DATA.VONMISES =  reshape(stressVONMISES,ngaus,[]) ;
    STRESS_DATA.MAX_VONMISES =   MAXstressVONMISES' ;
    
else
    STRESS_DATA.VONMISES  = [] ;
    STRESS_DATA.MAX_VONMISES  = [] ;
end



% DISPLACEMENTS LATERAL SURFACES
% -------------------------------
ndim = 3;
CONNECTb = [] ;

DATAIN = DefaultField(DATAIN,'PRINT_DISTRIBUTED_FORCES',0) ;

if DATAIN.PRINT_DISTRIBUTED_FORCES == 1
    DATAIN.POST_PROCESS_LATERAL_SURFACES.COARSE_MESH  = 0 ;
end

if  length(DATAIN.DOMAINS_POSTPROCESS)  ~=nDOM && DATAIN.POST_PROCESS_LATERAL_SURFACES.ACTIVE == 1
    % We plot the displacements of the lateral surfaces on GID
    % Nodes lateral surfaces
    % ----------------------
    if DATAIN.POST_PROCESS_LATERAL_SURFACES.COARSE_MESH == 1
        
        if isempty(DATAIN.POST_PROCESS_LATERAL_SURFACES.NAME_MESH_COARSE)
            NAME_MESH_COARSE = DATA_REFMESH.NameFileMeshLOC_coarse ;
            if isempty(DATA_REFMESH.CONNECTb_coarse)
                % Nothing is done
            else
                DATA_REFMESH.CONNECTb =  DATA_REFMESH.CONNECTb_coarse  ;
                DATA_REFMESH.TypeElementB = DATA_REFMESH.TypeElementB_coarse  ;
            end
        else
            NAME_MESH_COARSE = DATAIN.POST_PROCESS_LATERAL_SURFACES.NAME_MESH_COARSE ;
            error('Implement this option')
        end
        
        
    else
        % COONECTb is determined from DATA_REFMESH
        
    end
    
    
    
    %DISP3D_all
    DISP3D_lateral = 1;
    nfacesINTER = 2;
    
    %     if isempty(CONNECTb)
    %         CONNECTb = DATA_REFMESH.CONNECTb ;
    %     end
    
    ilateralSURF = (nfacesINTER+1):length( DATA_REFMESH.CONNECTb) ;
    CONNECTb = DATA_REFMESH.CONNECTb(ilateralSURF) ;
    NODES = cell2mat(CONNECTb') ;
    NODES = unique(NODES(:)) ;
    DOFS = small2large(NODES,ndim) ;
    DISP3D_new =cell(1,size(DISP3D,2)) ;
    SELECTED = zeros(1,size(DISP3D,2)) ;
    SELECTED(SELECTED_DOMAINS) = 1;
    for idom = 1:size(DISP3D,2)
        if SELECTED(idom) == 1
            DISP3D_new{idom} = DISP3D(:,idom) ;
        else
            DISP3D_new{idom} = DISP3D(DOFS,idom) ;
            
        end
    end
    DISP3D = cell2mat(DISP3D_new(:)) ;
    
    
    %     NODES_1domain = unique(cell2mat(CONNECTb')) ;
    %     DOFS = small2large(NODES_1domain,ndim) ;
    %     DISP3D_lateral.VALUE = DISP3D_all(DOFS,:) ;
    %  %   DISP3D_lateral.NODES_1domain =  NODES_1domain ;
    %     DISP3D_lateral.CONNECTb =  CONNECTb ;
else
    
    DISP3D_lateral =0;
    DISP3D = DISP3D(:,DATAIN.DOMAINS_POSTPROCESS) ;
    DISP3D = DISP3D(:) ;
    
end
end

% --------------------------------
