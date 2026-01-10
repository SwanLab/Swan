function [INPUTS_PARAMETERS,INFO_BASIC_DEFORMATIONAL_MODES] = Q8_for_beamBEHAVIOR_trainingTESTS(AMPLITUDE_alltests,include_warping,SURFACES,DATAOUT)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/01_BEAMQ8_ELAST.mlx

if nargin == 0
    load('tmp1.mat')
end

INPUTS_PARAMETERS = [] ;
INFO_BASIC_DEFORMATIONAL_MODES = [] ; 

% BEAM MODES, FIXED SURFACE = 1;
AMPLITUDE = AMPLITUDE_alltests*eye(7) ;

if include_warping == 0
    AMPLITUDE = AMPLITUDE(:,1:6);
end
for itest = 1:size(AMPLITUDE,2)
    INPUTS_PARAMETERS(itest).DIRICHLET =  DirichBound_BEAM_generic(DATAOUT.t0,DATAOUT.tEND,AMPLITUDE(:,itest),SURFACES)  ;
    INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = 'RIGID_BODY_WARPING_ON_SURFACES'; %
    
    
    
end



if include_warping == 1
    
    AMPLITUDE_warp = zeros(7,1);  AMPLITUDE_warp(7) = AMPLITUDE_alltests ;
    SURFACES = [SURFACES(2),SURFACES(1)] ;
    itest = length(INPUTS_PARAMETERS ) +1  ;
    INPUTS_PARAMETERS(itest).DIRICHLET =  DirichBound_BEAM_generic(DATAOUT.t0,DATAOUT.tEND,AMPLITUDE_warp,SURFACES)  ;
    INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = 'RIGID_BODY_WARPING_ON_SURFACES'; %
end

INFO_BASIC_DEFORMATIONAL_MODES.INDEX_FUNDAMENTAL_MODES = 1:(length(INPUTS_PARAMETERS)*length(DATAOUT.DATA_STEPS)) ;  % THESE ARE THE INDEXES OF THE

INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_FUNDAMENTAL_MODES = length(INPUTS_PARAMETERS) ; 
% SNAPSHOTS CORRESPONDING TO THE FUNDAMENTAL MODES (BEAM MODES IN THIS CASE)

%
% STUDY_ONLY_COMPLEMENTARY_MODES = 1 ;
%
%
% if   STUDY_ONLY_COMPLEMENTARY_MODES  == 1
%
%%%%%% COMPLEMENTARY MODES **************************************+
nstrain = 6;
MESH = DATAOUT.MESHcoarse;
MESH.DATA = [];
[~,~,~,~,~,~,GEOproperties,~,~,...
    ~,~] = ...
    GeometricMatricesFun(MESH,nstrain)  ;
ndim = size(MESH.COOR,2) ;
nDOFS = size(MESH.COOR,1)*ndim ;
ModesFaceGLO_all = [] ;
for isurfaceLOC = 1:length(SURFACES)
    isurface= SURFACES(isurfaceLOC) ;  % Index whose surface is to be rotated + warped
    rnodLOC=  MESH.NODES_FACES{isurface} ;  % Nodes pertaining to the surface
    DOFrLOC =  small2large(rnodLOC,ndim) ;  % DOFs associatted to the surface. ndim is the number of spatial dimensions
    COORface = MESH.COOR(rnodLOC,:) ;   % Spatial coordinates of the nodes
    Centroid = GEOproperties.FACES{isurface}.CENTROID;  % Centroid of the surface
    % 4. Coordinate relative to the centroid
    COORfaceREL = bsxfun(@minus,COORface',Centroid')' ;  %
    BasisRrb = ConstructBasisRigidBody(COORfaceREL) ; % Rigid body modes
    Rwarping = ConstructWarpingMode(COORfaceREL) ; % Warping mode
    if include_warping == 1
        ModesFace = [BasisRrb,Rwarping] ;
    else
        ModesFace =  BasisRrb ;
    end
    % GLOBAL COORDINATES
    ModesFaceGLO = zeros( (nDOFS),size(ModesFace,2)) ;
    ModesFaceGLO(DOFrLOC,:) = ModesFace ;
    ModesFaceGLO_all =[ModesFaceGLO_all,ModesFaceGLO] ;
end

% NOW WE COMPUTE THE ORTHOGONAL COMPLEMENT

AMPLITUDE_LINEAR_TESTS = eye(nDOFS) - ModesFaceGLO_all*((ModesFaceGLO_all'*ModesFaceGLO_all)\ModesFaceGLO_all') ;
[UU,SS,VV] = SVDT(AMPLITUDE_LINEAR_TESTS) ;
ncompl = nDOFS-size(ModesFaceGLO_all,2) ;
AMPLITUDE = AMPLITUDE_alltests*UU(:,1:ncompl) ;
ini  = length( INPUTS_PARAMETERS) ;
for itestLOC = 1:size(AMPLITUDE,2)
    itest = ini + itestLOC ;
    INPUTS_PARAMETERS(itest).DIRICHLET.AMPLITUDE = AMPLITUDE(:,itestLOC) ;
    INPUTS_PARAMETERS(itest).DIRICHLET.TIMEFUN = @(t) t/DATAOUT.tEND ;  %
    INPUTS_PARAMETERS(itest).DIRICHLET.INTERVAL =[DATAOUT.t0,DATAOUT.tEND];  %
    INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = ''; %
end



INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_BASIC_DEFORMATIONAL_MODES = length(INPUTS_PARAMETERS) ; 

disp(['Q8 EIF ELEMENT FOR  representing beam behavior'])
disp(['TOTAL NUMBER OF TESTS = ',num2str(length(INPUTS_PARAMETERS) ),' = ', 'number of deformational modes'])


% else
%     %---------------------------------------------------------------
%     AMPLITUDE = eye(24)*AMPLITUDE_alltests ;
%     ini  = length( INPUTS_PARAMETERS) ;
%     for itestLOC = 1:length(AMPLITUDE)
%         itest = ini + itestLOC ;
%         INPUTS_PARAMETERS(itest).DIRICHLET.AMPLITUDE = AMPLITUDE(:,itestLOC) ;
%         INPUTS_PARAMETERS(itest).DIRICHLET.TIMEFUN = @(t) t/DATAOUT.tEND ;  %
%         INPUTS_PARAMETERS(itest).DIRICHLET.INTERVAL =[DATAOUT.t0,DATAOUT.tEND];  %
%         INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = ''; %
%     end
%
% end

