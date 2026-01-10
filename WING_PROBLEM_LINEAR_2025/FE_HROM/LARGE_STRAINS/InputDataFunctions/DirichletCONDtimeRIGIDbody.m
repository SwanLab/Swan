function [U,a,DOFrLOC] = DirichletCONDtimeRIGIDbody(LOCVAR,DATA,ndim,MESH,GEOproperties)
%--------------------------------------------------------------------------
% [U,a,DOFrLOC] = DirichletCONDtimeRIGIDbody(LOCVAR,DATA,ndim,MESH,GEOproperties)
%
% PURPOSE:
%   Defines rigid-body type Dirichlet boundary conditions (translation +
%   rotation) for a single surface. The prescribed displacement field is
%   expressed in a space–time separated form:
%
%       dR(t) = U * a(t)
%
%   where:
%       - DOFrLOC : Local degrees of freedom (DOFs) on the surface
%       - U       : Spatial pattern(s) of rigid-body displacement (basis)
%       - a(t)    : Corresponding time-dependent amplitudes
%
%   Translations are defined at the surface centroid. Rotations are defined
%   relative to either the surface centroid (local) or an arbitrary global
%   center of rotation. The function supports both global and local axes.
%
% INPUT:
%   LOCVAR       - Struct describing prescribed motion (TRANSLATION,
%                  ROTATION, INTERVAL, ISLOCAL flags, etc.).
%   DATA         - Global data structure (time steps, etc.).
%   ndim         - Number of spatial dimensions (2D/3D).
%   MESH         - Mesh structure (contains nodes per surface, coordinates).
%   GEOproperties- Geometric properties (centroids, normals, tangents, etc.).
%
% OUTPUT:
%   U        - Spatial displacement modes (size = nDOF x nModes).
%   a        - Time-dependent amplitudes (size = nModes x nSteps).
%   DOFrLOC  - Vector of DOF indices constrained by this rigid body motion.
%
% NOTES:
%   - Motion is decomposed via SVD to obtain a compact space–time form.
%   - Handles compound rigid-body motions (global + local).
%   - Archetypal input: 
%       INPUTS_ROTATION.m
%   - Documentation:
%       DOCS/IMPLEMENTATION_STATIC.pdf
%
% JAHO - Date not specified
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end

% --- Nodes and DOFs on prescribed surface --------------------------------
isurface  = LOCVAR.NUMBER_SURFACE;               % Surface ID
rnodLOC   = MESH.NODES_FACES{isurface};          % Node list on this surface
DOFrLOC   = small2large(rnodLOC,ndim);           % DOFs on surface
COORface  = MESH.COOR(rnodLOC,:);                % Node coordinates
Centroid  = GEOproperties.FACES{isurface}.CENTROID;

% Relative coordinates with respect to centroid
COORfaceREL = bsxfun(@minus,COORface',Centroid')';

% --- Initialize ----------------------------------------------------------
nloads      = length(LOCVAR.PRESCRIBED_DISP);    % Number of load states
ntimesteps  = length(DATA.STEPS);
U = cell(1,nloads); 
a = cell(nloads,1);

% --- Loop over load states ------------------------------------------------
for iload = 1:nloads
    TRANSLATION = LOCVAR.PRESCRIBED_DISP(iload).TRANSLATION;
    ROTATION    = LOCVAR.PRESCRIBED_DISP(iload).ROTATION; % Rotation rel. to centroid
    INTERVAL    = LOCVAR.PRESCRIBED_DISP(iload).INTERVAL;
    
    LOC_DISPF = LOCVAR.PRESCRIBED_DISP(iload);
    LOC_DISPF = DefaultField(LOC_DISPF,'ISLOCAL',0);
    
    % --- Local/global axes handling ---------------------------------------
    if LOC_DISPF.ISLOCAL == 0
        % Motion defined in global axes
        ROTATION.RotationVectorGlobal = [];
    else
        % Motion defined in local surface axes
        % Build local frame (normal + tangents at reference Gauss point)
        igausREF    = 1;
        UnitNormal  = GEOproperties.FACES{isurface}.UnitNormalAtGaussPoint(:,igausREF);
        UnitTangent1= GEOproperties.FACES{isurface}.UnitTangentAtGaussPoint{1}(:,igausREF);
        UnitTangent2= GEOproperties.FACES{isurface}.UnitTangentAtGaussPoint{2}(:,igausREF);
        RotMatFace  = [UnitNormal,UnitTangent1,UnitTangent2];
        
        % Align local axes with principal inertia directions
        RotMatFace = AxesPrincipalInertia(RotMatFace,...
            GEOproperties.FACES{isurface}.COORrelA_global',...
            GEOproperties.FACES{isurface}.GeometricMassMatrix);
        
        % Convert local prescribed vectors to global frame
        ROTATION.RotationVectorGlobal = RotMatFace * ROTATION.ROTATION_VECTOR_MAXIMUM(:);
        TRANSLATION.AMPLITUDE         = RotMatFace * TRANSLATION.AMPLITUDE(:);
    end
    
    % --- Additional global rotation/translation ---------------------------
    PRESCRIBED_DISP_LOC = LOCVAR.PRESCRIBED_DISP(iload);
    PRESCRIBED_DISP_LOC = DefaultField(PRESCRIBED_DISP_LOC,'ROTATION_GLO',[]);
    PRESCRIBED_DISP_LOC = DefaultField(PRESCRIBED_DISP_LOC,'TRANSLATION_GLO',[]);
    ROTATION_GLO        = PRESCRIBED_DISP_LOC.ROTATION_GLO;
    TRANSLATION_GLO     = PRESCRIBED_DISP_LOC.TRANSLATION_GLO;
    
    if ~isempty(ROTATION_GLO)
        ROTATION_GLO = DefaultField(ROTATION_GLO,'CENTER',[]);
        if isempty(ROTATION_GLO.CENTER.COOR)
            % Use centroid of a reference surface as center of rotation
            isurfGLO = PRESCRIBED_DISP_LOC.ROTATION_GLO.CENTER.ISCENTROID;
            CENTERglo= GEOproperties.FACES{isurfGLO}.CENTROID;
        else
            CENTERglo= ROTATION_GLO.CENTER.COOR;
        end
    end
    
    % --- Displacement evolution over time --------------------------------
    D = zeros(length(DOFrLOC),ntimesteps);
    for istep = 1:ntimesteps
        TIMELOC = DATA.STEPS(istep);
        Dloc    = D(:,istep); % Initially zero
        
        % Rigid-body motion relative to surface centroid
        Dloc = RigidBodyMotionALL(TRANSLATION,INTERVAL,TIMELOC,...
                                   rnodLOC,Dloc,ndim,ROTATION,COORfaceREL);
        
        % Global rotation/translation (compound motion)
        if ~isempty(ROTATION_GLO)
            dLOCAL = reshape(Dloc,ndim,[]);
            COORfaceREL_new = bsxfun(@minus,COORface',CENTERglo')';
        else
            COORfaceREL_new = [];
            dLOCAL = [];
        end
        
        if ~isempty(ROTATION_GLO) || ~isempty(TRANSLATION_GLO)
            Dloc = RigidBodyMotionALL_compound(TRANSLATION_GLO,INTERVAL,TIMELOC,...
                                               rnodLOC,Dloc,ndim,ROTATION_GLO,...
                                               COORfaceREL_new,dLOCAL);
        end
        D(:,istep) = Dloc;
    end
    
    % --- Space–time separation via SVD -----------------------------------
    if any(any(D))
        [Uloc,S,V] = RSVDT(D);          % Reduced SVD of displacement matrix
        aloc = bsxfun(@times,V',S);     % Amplitudes
        U{iload} = Uloc;
        a{iload} = aloc;
    else
        % No motion -> zero displacement
        U{iload} = zeros(size(D,1),1);
        a{iload} = zeros(1,size(D,2));
    end
end

% --- Assemble outputs ----------------------------------------------------
U = cell2mat(U);
a = cell2mat(a);

% BEFORE BEING COMMENTED BY CHATGPT 5, 24th August 2025

% function      [U,a,DOFrLOC] =    DirichletCONDtimeRIGIDbody(LOCVAR,DATA,ndim,MESH,GEOproperties) ;
% % DIRICHLET BOUNDARY CONDITIONS, RIGID BODY TYPE (ONE SINGLE SURFACE)
% % See archetypal input file
% % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/INPUTS_ROTATION.m
% % See DOCS/IMPLEMENTATION_STATIC.pdf
% %
% % iloadstate = 1;  % Loading stage
% % LOCVAR.PRESCRIBED_DISP = [] ;
% % disp_max = 0;
% % disp_x = disp_max*cosd(DATAINPUT.PARAM); disp_y = disp_max*sind(DATAINPUT.PARAM) ;
% %
% %
% % LOCVAR.PRESCRIBED_DISP(iloadstate).TRANSLATION.AMPLITUDE = [disp_x,disp_y] ;  % % Translation (centroid face)
% % LOCVAR.PRESCRIBED_DISP(iloadstate).TRANSLATION.TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
% %
% % LOCVAR.PRESCRIBED_DISP(iloadstate).ROTATION.X.MAXANGLE = 0 ;     % % Rotation around the x-axix
% % LOCVAR.PRESCRIBED_DISP(iloadstate).ROTATION.X.TIMEFUN = @(t) 0 ;  ;     % % Translation (centroid face)
% %
% % LOCVAR.PRESCRIBED_DISP(iloadstate).ROTATION.Y.MAXANGLE = 0 ;     % % Rotation around the x-axix
% % LOCVAR.PRESCRIBED_DISP(iloadstate).ROTATION.Y.TIMEFUN = @(t) 0 ;     % % Translation (centroid face)
% %
% % LOCVAR.PRESCRIBED_DISP(iloadstate).ROTATION.Z.MAXANGLE = 0 ;     % % Rotation around the x-axix
% % LOCVAR.PRESCRIBED_DISP(iloadstate).ROTATION.Z.TIMEFUN = @(t) t/tEND ;     % % Translation (centroid face)
% %
% if nargin == 0
%     load('tmp1.mat')
% end
% 
% %%%%
% isurface= LOCVAR.NUMBER_SURFACE ;
% % 1. Nodes comprising the surface
% rnodLOC=  MESH.NODES_FACES{isurface};
% DOFrLOC =  small2large(rnodLOC,ndim) ;
% COORface = MESH.COOR(rnodLOC,:) ;
% 
% 
% Centroid = GEOproperties.FACES{isurface}.CENTROID;
% 
% 
% % 4. Coordinate relative to the centroid
% COORfaceREL = bsxfun(@minus,COORface',Centroid')' ;
% 
% 
% nloads = length( LOCVAR.PRESCRIBED_DISP);
% % Uloc = zeros(length(DOFrLOC),nloads) ;
% % aloc = zeros(nloads,length(DATA.STEPS)) ;
% ntimesteps = length(DATA.STEPS) ;
% 
% U = cell(1,nloads) ;
% a = cell(nloads,1) ;
% 
% for  iload = 1:nloads
%     TRANSLATION= LOCVAR.PRESCRIBED_DISP(iload).TRANSLATION ;
%     ROTATION= LOCVAR.PRESCRIBED_DISP(iload).ROTATION ; % Rotation relative to the centroid of the surface
%     
%     INTERVAL =  LOCVAR.PRESCRIBED_DISP(iload).INTERVAL ;
%     
%     LOC_DISPF  = LOCVAR.PRESCRIBED_DISP(iload) ;
%     
%     LOC_DISPF = DefaultField(LOC_DISPF,'ISLOCAL',0) ;
%     
%     if LOC_DISPF.ISLOCAL ==0
%         % Global (axes X, Y, Z)
%         ROTATION.RotationVectorGlobal  = [] ;
%     else
%         % Local coordinates
%         % ------------------
%         % Rotation matrix associated to the surface
%         
%         % All unit normals
%         % FACES{iface}.UnitNormalAtGaussPoint
%         
%         
%         % Unit normal of the first Gauss point
%         % ------------------------------------ (reference)
%         igausREF  = 1 ;
%         UnitNormal = GEOproperties.FACES{isurface}.UnitNormalAtGaussPoint(:,igausREF) ;
%         UnitTangent1 = GEOproperties.FACES{isurface}.UnitTangentAtGaussPoint{1}(:,igausREF) ;
%         UnitTangent2 = GEOproperties.FACES{isurface}.UnitTangentAtGaussPoint{2}(:,igausREF) ;
%         % Rotation matrix
%         RotMatFace = [UnitNormal,UnitTangent1,UnitTangent2] ;
%         
%         % ROTATIONS
%         % ----------------
%         % Rotation vector
%         
%         %
%         RotMatFace = AxesPrincipalInertia(RotMatFace,...
%             GEOproperties.FACES{isurface}.COORrelA_global',...
%             GEOproperties.FACES{isurface}.GeometricMassMatrix) ;
%         
%         RotationVectorGlobal = RotMatFace*ROTATION.ROTATION_VECTOR_MAXIMUM(:) ;
%         
%         ROTATION.RotationVectorGlobal = RotationVectorGlobal ;
%         
%         
%         
%         
%         
%         
%         
%         
%         
%         TRANSLATION.AMPLITUDE = RotMatFace*TRANSLATION.AMPLITUDE(:) ;
%         
%         
%     end
%     
%     
%     % 3.Global rotation (for defining rotation)
%     %*******************************************
%     PRESCRIBED_DISP_LOC = LOCVAR.PRESCRIBED_DISP(iload) ;
%     PRESCRIBED_DISP_LOC = DefaultField(PRESCRIBED_DISP_LOC,'ROTATION_GLO',[]) ;
%     PRESCRIBED_DISP_LOC = DefaultField(PRESCRIBED_DISP_LOC,'TRANSLATION_GLO',[]) ;
%     ROTATION_GLO = PRESCRIBED_DISP_LOC.ROTATION_GLO ; % Rotation relative to a given point (CENTER)
%     TRANSLATION_GLO = PRESCRIBED_DISP_LOC.TRANSLATION_GLO ;
%     
%     if ~isempty(ROTATION_GLO)
%         ROTATION_GLO = DefaultField(ROTATION_GLO,'CENTER',[]) ;
%         if isempty(ROTATION_GLO.CENTER.COOR)
%             isurfGLO = PRESCRIBED_DISP_LOC.ROTATION_GLO.CENTER.ISCENTROID  ; % Number of surface whose centroid is the center of rotation
%             CENTERglo  =  GEOproperties.FACES{isurfGLO}.CENTROID; % Coordinates center of rotation
%         else
%             CENTERglo =ROTATION_GLO.CENTER.COOR ;
%         end
%     end
%     
%     D= zeros(length(DOFrLOC),ntimesteps)   ;
%     
%     for istep = 1:ntimesteps
%         
%         TIMELOC = DATA.STEPS(istep) ;
%         
%         
%         Dloc =  D(:,istep)  ; % This is zero initially
%         % RIGID BODY MOTION RELATIVE TO THE CENTROID OF THE FACE UNDER
%         % CONSIDERATION
%         Dloc = RigidBodyMotionALL(TRANSLATION,INTERVAL,TIMELOC,rnodLOC,Dloc,ndim,ROTATION,COORfaceREL) ;
%         % D(:,istep)  = Dloc ;
%         %
%         
%         %% GLOBAL ROTATION
%         % --------------------
%         if ~isempty(ROTATION_GLO)
%             % Global rotation
%             % ------------------
%             %  dLOCAL = D(:,istep) ; %
%             dLOCAL = reshape(Dloc,ndim,[]) ;
%             %   COORface_NEW = dLOCAL + COORface' ; %
%             COORfaceREL_new = bsxfun(@minus,COORface',CENTERglo')' ;
%             %  Dloc = zeros(size(Dloc));
%         else
%             COORfaceREL_new = [] ;
%         end
%         if ~isempty(ROTATION_GLO) ||      ~isempty(TRANSLATION_GLO)
%             Dloc = RigidBodyMotionALL_compound(TRANSLATION_GLO,INTERVAL,TIMELOC,rnodLOC,Dloc,ndim,ROTATION_GLO,COORfaceREL_new,dLOCAL) ;
%             
%         end
%         
%         D(:,istep)  = Dloc ;
%         
%         
%         
%         
%         
%     end
%     
%     % Next we apply the SVD
%     if any(any(D))
%         [Uloc,SS,VV] =  RSVDT(D) ;
%         aloc = bsxfun(@times,VV',SS) ;
%         U{iload} =  Uloc ;
%         a{iload} =  aloc ;
%     else
%         
%         U{iload}  = zeros(size(D,1),1) ;
%         a{iload} = zeros(1,size(D,2)) ;
%     end
%     
%     
%     %     U{izone}(:,iload) = displ{iload} ;    % Uniform Boundary Conditions
%     %     FactorSteps =  LOCVAR.PRESCRIBED_DISP(iload).TIMEFUN(DATA.STEPS) ;
%     %     LIMITS_interval =  LOCVAR.PRESCRIBED_DISP(iload).INTERVAL ;
%     %     FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
%     %     a{izone}(iload,:) = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
% end
% 
% 
% U = cell2mat(U) ;
% a = cell2mat(a) ;
% 
% 
% 
% 
