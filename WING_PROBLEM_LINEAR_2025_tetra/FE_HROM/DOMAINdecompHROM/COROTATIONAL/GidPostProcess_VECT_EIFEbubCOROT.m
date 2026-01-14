function  strainCOARSE_history = GidPostProcess_VECT_EIFEbubCOROT(MESH,DATA,icluster,DATAinpGID,PROPMAT,DATALOC,...
    QrotTIME,QrotINI,LboolCall)
%
%   This function processes the results of a coarse-scale simulation using the
%   Empirical Interscale Finite Element (EIFE) method with bubble enrichment under 
%   a co-rotational framework. It prepares mesh and result files for GiD visualization, 
%   interpolates displacements (including midpoints if needed), and reconstructs 
%   strain and stress tensors.
%
%   INPUT:
%     MESH         - Mesh structure with geometry, connectivity, and optional auxiliary mesh
%     DATA         - Structure containing global printing options, file names, steps, etc.
%     icluster     - Index of the current cluster being postprocessed
%     DATAinpGID   - Optional input data specific to GiD; can be empty
%     PROPMAT      - Material properties
%     DATALOC      - Local configuration options for postprocessing
%     QrotTIME     - Rotation matrices at each time step (for co-rotational correction)
%     QrotINI      - Initial rotation matrices
%     LboolCall    - Logical flag to determine whether to call GiD postprocessing or not
%
%   OUTPUT:
%     strainCOARSE_history - Strain field history across all time steps, ready for GiD
%
%   MAIN TASKS PERFORMED:
%   ---------------------------------------------------------------------------
%   1. Load mesh data, rotation matrices, and stored snapshots (with or without SVD compression)
%   2. Decode and reconstruct displacement fields
%   3. (If needed) Restore midpoint data for 26-node hexahedral elements to comply with GiD
%   4. Extract internal variable fields (e.g., PK2 stress, von Mises, strain history)
%   5. Optionally smooth displacements using mapping fields from uncoupled to support mesh
%   6. Pass all data to the final postprocessing driver (PostProcess_EIFE_vectBUBcorot)
%
%   NOTES:
%   - The function assumes a co-rotational formulation and handles rotation recovery explicitly.
%   - GiD mesh (`*.msh`) and result (`*.res`) files are written as defined in DATA.PRINT.
%   - For hex27-to-hex26 mesh simplification, interior nodes are recomputed and added.
%
%   AUTHOR:
%     Joaquín A. Hernández Ortega (JAHO), UPC BarcelonaTech, 5-Nov-2024
%     Adapted from GidPostProcess_VECT_EIFEbub.m for co-rotational corrections
%     Comments by CHatGPT-4
% -----------------------------------------------------------------------------


if nargin==0
    load('tmp.mat')
elseif nargin == 3
    DATAinpGID = [] ;
    
end
COOR = MESH.COOR; % COARSE-SCALE MESH
CN = MESH.CN ; % CONNECTIVITIES

MESH = DefaultField(MESH,'AUXILIAR',[]);
if isempty(MESH.AUXILIAR)
    CN = MESH.CN ;
    TypeElement = MESH.TypeElement;
    MaterialType_coarse_GID = MESH.MaterialType;
    
else
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
    CN = MESH.AUXILIAR.CN  ;
    TypeElement = MESH.AUXILIAR.TypeElement  ;
    MaterialType_coarse_GID = ones(size(CN,1),1) ;
end

MaterialType = MESH.MaterialType;  % Type of EIF element


NameFile_msh = DATA.PRINT.NAME_FILE_MSH{icluster} ;
NameFile_res= DATA.PRINT.NAME_FILE_RES{icluster} ;

if exist(DATA.STORE.NAME_MATFILE_STORE{icluster})
    if DATA.STORE.COMPRESS_WITH_SVD == 1
        load( DATA.STORE.NAME_MATFILE_STORE{icluster},'SNAP_cluster','STEP_LOC') ;
        % Typically
        % SNAP_cluster =
        %
        %   struct with fields:
        %
        %                       DISP: [1×1 struct]
        %                      RESID: [1×1 struct]
        %              STRAIN_ENERGY: [1×1 struct]
        %     VONMISES_CAUCHY_STRESS: [1×1 struct]
    else
        load( DATA.STORE.NAME_MATFILE_STORE{icluster},'SNAP','STEP_LOC') ;
    end
else
    return
end


[STEPS_PRINT,iaaa,ibbb ]= intersect(STEP_LOC,DATA.PRINT.NSTEPS) ; % Steps to print

if length(STEPS_PRINT) ~= length(STEP_LOC)
    error('7-Nov-2024:There is an unresolved issue when not all time steps are printed... Set up for the time being FREQ = 1')
end

if DATA.STORE.COMPRESS_WITH_SVD == 1
    %Uncompressing information.  SNAP = U*S*V^T
    % Be wary of almost zero matrices....
    fff= fieldnames(SNAP_cluster) ;
    SNAP = [];
    for  iii = 1:length(fff)
        if ~isempty( SNAP_cluster.(fff{iii}).V)
            V = SNAP_cluster.(fff{iii}).V(iaaa,:) ;
            S = SNAP_cluster.(fff{iii}).S  ;
            SV = bsxfun(@times,V',S)';
            SNAP.(fff{iii}) = SNAP_cluster.(fff{iii}).U*SV' ;
            
        else
            SNAP.(fff{iii}) = [];
        end
        
    end
else
    fff= fieldnames(SNAP) ;
    for  iii = 1:length(fff)
        SNAP.(fff{iii}) =   SNAP.(fff{iii})(:,iaaa) ;
    end
    
end
TIME_PRINT = DATA.STEPS(STEPS_PRINT) ;
DATA.TIME_PRINT = TIME_PRINT ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = SNAP.DISP;%(:,STEPS_PRINT) ; Coarse-scale DOFs, including bubble DOFs

OTHERinputs.COORprintCOARSE = [] ; OTHERinputs.CNprintCOARSE = [] ; OTHERinputs.dPRINTcoarse = [] ;

if isfield(MESH,'SMOOTH_from_UNCOUP_TO_SUPPORT')
    OTHERinputs.dPRINTcoarse = MESH.SMOOTH_from_UNCOUP_TO_SUPPORT*d ;
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/06_UNCOUPLED.mlx
    OTHERinputs.COORprintCOARSE =  MESH.COOR_SUPPORT ;
    OTHERinputs.CNprintCOARSE = MESH.CN_SUPPORT ;
    
end



if size(CN,2) == 26
    % Hexahedra element in which we have eliminated the interior point
    % We need now to restore such a point in order to be processesd by GID
    COORintPOINT = zeros(size(CN,1),size(COOR,2)) ;
    nnodeE = 26 ;
    ndim = 3;
    for  inodeE = 1:nnodeE
        COORintPOINT  = COORintPOINT + COOR(CN(:,inodeE),:) ;
    end
    COORintPOINT = COORintPOINT/nnodeE ;
    % DISPLACEMENTS
    dMIDPOINT = zeros(size(COORintPOINT,1)*ndim,size(d,2)) ;
    for itime = 1:size(d,2)
        dLOC = reshape(d(:,itime),ndim,[] )' ;
        dMIDPOINT_loc = zeros(size(CN,1),ndim) ;
        for  inodeE = 1:nnodeE
            dMIDPOINT_loc   = dMIDPOINT_loc + dLOC(CN(:,inodeE),:) ;
        end
        dMIDPOINT_loc = dMIDPOINT_loc'/nnodeE;
        dMIDPOINT(:,itime) =  dMIDPOINT_loc(:) ;  ;
    end
    
    nnodeOLD  = size(COOR,1) ;
    OTHERinputs.COORprintCOARSE = [COOR;COORintPOINT] ;
    CNintpoint = nnodeOLD + (1:size(CN,1))';
    OTHERinputs.CNprintCOARSE = [CN,CNintpoint] ;
    OTHERinputs.dPRINTcoarse = [d; dMIDPOINT] ;
    
end

if isfield(SNAP,'PK2STRESS')
    stressGLO_REF = SNAP.PK2STRESS; %(:,STEPS_PRINT);  % These are the stresses in the unrotated configuration
else
    stressGLO_REF = [] ;
end
DATA.NameFile_msh = NameFile_msh;
DATA.NameFile_res = NameFile_res;
SNAP = DefaultField(SNAP,'VONMISES_CAUCHY_STRESS',[]) ;

if ~isempty(SNAP.VONMISES_CAUCHY_STRESS)   && isempty(MESH.AUXILIAR)
    VON_mises_COARSE= zeros(size(CN,1),size(SNAP.VONMISES_CAUCHY_STRESS,2));
    for ielem = 1:size(CN,1)
        % VONloc  = reshape(SNAP.VONMISES_CAUCHY_STRESS(:,itime),DATA.MESH.ngaus_STRESS,[]) ;
        VONloc = SNAP.VONMISES_CAUCHY_STRESS(DATA.MESH.IndexECMpoints_per_element{ielem},:);
        VON_mises_COARSE(ielem,:) = max(VONloc,[],1)  ;
    end
    VON_mises_COARSE = VON_mises_COARSE; %(:,STEPS_PRINT);
else
    VON_mises_COARSE = [] ;
end

SNAP = DefaultField(SNAP,'InternalVarStrain',[]) ;
if ~isempty(SNAP.InternalVarStrain)
    InternalVarStrain= zeros(size(CN,1),size(SNAP.InternalVarStrain,2));
    for ielem = 1:size(CN,1)
        VONloc  =   SNAP.InternalVarStrain(DATA.MESH.IndexECMpoints_per_element{ielem},:); % reshape(SNAP.InternalVarStrain(:,itime),DATA.MESH.ngaus_STRESS,[]) ;
        InternalVarStrain(:,itime) = max(VONloc,[],1)  ;
    end
    InternalVarStrain = InternalVarStrain(:,STEPS_PRINT) ;
else
    InternalVarStrain = [] ;
end





DATALOC = DefaultField(DATALOC,'FineScaleDomainsToShow',[])  ;
DATA.FineScaleDomainsToShow = DATALOC.FineScaleDomainsToShow  ;
DATALOC = DefaultField(DATALOC,'ReconstructionStressesUsingCoarseScaleDISP',1) ;
DATA.ReconstructionStressesUsingCoarseScaleDISP = DATALOC.ReconstructionStressesUsingCoarseScaleDISP  ;
OTHERinputs.MaterialType_coarse_GID = MaterialType_coarse_GID;

strainCOARSE_history = PostProcess_EIFE_vectBUBcorot(COOR,CN, PROPMAT,MaterialType,...
    TypeElement,DATA,MESH.TRANSF_COORD,d,stressGLO_REF,VON_mises_COARSE,InternalVarStrain,OTHERinputs,...
    QrotTIME,QrotINI,LboolCall);




