function Ftrac = NeumannCONDtime(NEUMANN,DATA,ndim,MESH,NstT_W_N_boundaries,GEOproperties)
%--------------------------------------------------------------------------
% function Ftrac = NeumannCONDtime(NEUMANN, DATA, ndim, MESH, ...
%                                   NstT_W_N_boundaries, GEOproperties)
%
% PURPOSE:
%   Constructs the global traction vector `Ftrac` due to Neumann boundary
%   conditions (e.g., surface loads or point forces) using a space-time
%   separated representation:
%
%       Ftrac = Ftrac.U * Ftrac.a(t)
%
%   This separation allows efficient handling of time-dependent loads in
%   reduced-order models and large-scale simulations.
%
% INPUTS:
%   - NEUMANN : Structure array with boundary conditions. Each entry may
%               specify either:
%                  - NUMBER_SURFACE with FORCE_PER_UNIT_SURFACE or
%                    GENERALIZED_POINT_FORCE
%                  - NUMBER_NODE with FORCE (nodal force)
%   - DATA    : Structure containing global simulation information.
%               Must include field `DATA.STEPS`, a vector of time steps.
%   - ndim    : Number of spatial dimensions (2 or 3).
%   - MESH    : Mesh structure with fields:
%                  - COOR: Nnodes x ndim matrix of nodal coordinates
%                  - CNb, Indexes_faces_bnd_element: boundary connectivity
%                  - NormalBoundaryElementsFace, TangentBoundaryElementsFace
%   - NstT_W_N_boundaries : Cell array containing interpolation matrices
%                           for boundary integrals (one per boundary face).
%   - GEOproperties : Geometry-specific information (not used directly).
%
% OUTPUT:
%   - Ftrac : Structure with fields:
%              - Ftrac.U : Spatial patterns of the applied loads (sparse matrix)
%              - Ftrac.a : Temporal evolution of the load amplitudes (matrix)
%
% REMARKS:
%   - The routine handles both distributed surface tractions and
%     generalized point forces (e.g., Dirac-like forces).
%   - Forces can be defined in global or local coordinates (`ISLOCAL`).
%   - Time dependency is introduced via user-defined `TIMEFUN` functions
%     inside each force specification.
%   - The force distribution over surfaces is assumed uniform unless local
%     coordinate transformations are specified.
%
% REFERENCES:
%   - See implementation details in:
%       DOCS/01_LargeStrainsCode.pdf, page 26
%   - For generalized point forces:
%       /.../FIBREGY_PROJECT_2022/03_CYLINDRICAL_TOWER/ImplemetationGeneralizedForces.mlx
%
% HISTORY:
%   - Initial version: JAHO, 23-Nov-2020
%   - Update for generalized point forces: 11-May-2022
%   - Patch for empty Ftrac: 31-May-2024
%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal. Determine Ftrac    in a space-time separated fashion
%
%  Ftrac = Ftrac.U*Ftrac.a(t)
%  See          DOCS/01_LargeStrainsCode.pdf, page 26
%
% JAHO- 23-NOV-2O2O
% 11-May-2022:  Adding "generalized point forces" over surfaces, see
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/FIBREGY_PROJECT_2022/03_CYLINDRICAL_TOWER/ImplemetationGeneralizedForces.mlx
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

nzones = length(NEUMANN) ; % number of zones with prescribed forces
U = cell(nzones,1) ;  % This is for the space/time decomposition. U contains the spatial patterns
a = cell(nzones,1) ; % while contains the temporal patterns (functions )
nnodeALL = size(MESH.COOR,1) ;
for  izone = 1:length(NEUMANN)
    %%%%
    if ~isempty(NEUMANN(izone).NUMBER_SURFACE)
        iface = NEUMANN(izone).NUMBER_SURFACE ; % Number of face
        if iface > length(MESH.Indexes_faces_bnd_element)
            error(['Surface ',num2str(iface),' has not been defined'])
        end
        IND_CNbLOC = MESH.Indexes_faces_bnd_element{iface} ;  %
        CNbLOC = MESH.CNb(IND_CNbLOC,:)  ; % Connectivities elements of the surface
        nnodeEB = size(CNbLOC,2) ;
        NstT_W = NstT_W_N_boundaries{iface} ;
        %
        if isfield(NEUMANN(izone),'FORCE_PER_UNIT_SURFACE')
            FIELDLOC = 'FORCE_PER_UNIT_SURFACE' ;
        elseif isfield(NEUMANN(izone),'GENERALIZED_POINT_FORCE')
            FIELDLOC = 'GENERALIZED_POINT_FORCE' ;
          elseif isfield(NEUMANN(izone),'GAUSSIAN_TRACTION')  
             FIELDLOC = 'GAUSSIAN_TRACTION' ;
                 elseif isfield(NEUMANN(izone),'UNIFORM_TRACTION')  
             FIELDLOC = 'UNIFORM_TRACTION' ;
        else
            error('Option not implemented ')
        end
        NORMALS = MESH.NormalBoundaryElementsFace{iface} ;
        TANGENTS = MESH.TangentBoundaryElementsFace{iface} ;
    elseif ~isempty(NEUMANN(izone).NUMBER_NODE)
        FIELDLOC = 'FORCE' ;
        NODES = NEUMANN(izone).NUMBER_NODE ;
    else
        error('Option not implemented')
    end
    nloads=  length(NEUMANN(izone).(FIELDLOC)) ;
    U{izone} = sparse(nnodeALL*ndim,nloads) ;
    a{izone} = zeros(nloads,length(DATA.STEPS)) ;
    
    
    for iload = 1:nloads
        FORCE_INFO =NEUMANN(izone).(FIELDLOC)(iload);
        % We assume that the traction vector (force per unit surface)
        % is uniform over the surface.
        FactorSteps = FORCE_INFO.TIMEFUN(DATA.STEPS) ;
        LIMITS_interval =  FORCE_INFO.INTERVAL ;
        %  FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
        %  mistake....
        
        FactorSteps = FactorSteps.*(DATA.STEPS > LIMITS_interval(1)) ;
        switch FIELDLOC
            case {'FORCE_PER_UNIT_SURFACE','GENERALIZED_POINT_FORCE','GAUSSIAN_TRACTION','UNIFORM_TRACTION'}
                FORCE_INFO = DefaultField(FORCE_INFO,'ISLOCAL',0);
                ISLOCAL = FORCE_INFO.ISLOCAL ;
                
        end
        a{izone}(iload,:) = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
        
        
        switch FIELDLOC
            case {'FORCE_PER_UNIT_SURFACE','FORCE'}
                Floc = sparse(nnodeALL*ndim,1) ;
                for idim = 1:ndim
                    switch FIELDLOC
                        case 'FORCE'
                            AmplitudeIDIR    =  FORCE_INFO.AMPLITUDE ;
                            DOFs = small2large(NODES,ndim)  ;
                            for idim = 1:ndim
                                DOFSloc = DOFs(idim:ndim:end) ;
                                Floc(DOFSloc) = AmplitudeIDIR(idim)  ;
                            end
                        case 'FORCE_PER_UNIT_SURFACE'
                            if ISLOCAL == 0
                                AmplitudeIDIR    =  FORCE_INFO.AMPLITUDE(idim) ;
                                TnodLOC= AmplitudeIDIR*ones(size(CNbLOC))' ;
                                TnodVECT  =  TnodLOC(:) ;
                            else
                                % KW:NormalsLocalCoordinateForces
                                % Local coordinates to global coordinates
                                if ndim == 2
                                    Telem =NORMALS(idim,:)*FORCE_INFO.AMPLITUDE(1) + TANGENTS(idim,:)*FORCE_INFO.AMPLITUDE(2) ;
                                else
                                    Telem =NORMALS(idim,:)*FORCE_INFO.AMPLITUDE(1) + TANGENTS{1}(idim,:)*FORCE_INFO.AMPLITUDE(2) + ...
                                        TANGENTS{3}(idim,:)*FORCE_INFO.AMPLITUDE(3);
                                end
                                TnodeELEM = repmat(Telem,nnodeEB,1) ;
                                TnodVECT = TnodeELEM(:);
                            end
                            
                            Fdim = NstT_W*TnodVECT ;
                            INDEXES = idim:ndim:length(Floc) ;
                            Floc(INDEXES) = Fdim ;
                    end
                end
                
            case 'GENERALIZED_POINT_FORCE'
                
                
                
                Floc =  PointForceGeneralized(nnodeALL,ndim,ISLOCAL,FORCE_INFO,NstT_W,MESH,iface) ;
                
            case 'GAUSSIAN_TRACTION'
                
                
                
                Floc =  GaussianTraction(nnodeALL,ndim,ISLOCAL,FORCE_INFO,NstT_W,MESH,iface) ;
                
            case 'UNIFORM_TRACTION'
                 Floc =  UniformTraction(nnodeALL,ndim,ISLOCAL,FORCE_INFO,NstT_W,MESH,iface) ;
                
            otherwise
                
                
        end
        
        U{izone}(:,iload) =  Floc ;
    end
    
end

Ftrac.U = cell2mat(U') ;
Ftrac.a = cell2mat(a) ;

if isempty(Ftrac.U)
    % 31-May-2024
    % -------------
    Ftrac.U = zeros(DATA.MESH.ndof,1) ;
    Ftrac.a  = ones(size(DATA.STEPS)) ;
end

PLOT_AMPL = 0;
if  PLOT_AMPL == 1
    
    figure(200)
    hold on
    for  iii  = 1:size(Ftrac.a,1)
        
        plot(Ftrac.a(iii,:))
        
    end
    
end

%
% dR.U = U ;
% %dR.a = a ;

% DOFr  =cell2mat(DOFr) ;
% [DOFr,III ]= sort(DOFr)
%
% dR.U = blkdiag(U{:}) ;
% dR.a = cell2mat(a) ;
%
% dR.U = dR.U(III,:) ;
%



%


%
% % Removed repeated condions
% % ---------------------------
% for idim = 1:ndim
%     [rnod{idim}, AAA] = unique(rnod{idim}) ;
%     uPRES{idim} = uPRES{idim}(AAA) ;
% end
%
% % Degrees of freedom and prescribed displacements
% DOFr = [] ; dR = [] ;
% for idim = 1:ndim
%     DOFr = [DOFr ; (rnod{idim}-1)*ndim+idim];
%     dR = [dR ; uPRES{idim}];
% end