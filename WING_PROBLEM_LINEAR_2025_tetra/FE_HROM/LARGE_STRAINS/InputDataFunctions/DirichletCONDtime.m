function [DOFr,dR] = DirichletCONDtime(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
%--------------------------------------------------------------------------
% [DOFr, dR] = DirichletCONDtime(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
%
% PURPOSE:
%   Defines the prescribed Dirichlet boundary conditions in a time-dependent
%   manner, returning the list of constrained degrees of freedom (DOFr) and
%   the associated space–time separated representation of the prescribed 
%   displacement field.
%
%   The displacement field is expressed as:
%       dR(t) = U * a(t)
%   where:
%       - DOFr = vector of prescribed DOF indices
%       - U    = spatial pattern(s) of the prescribed displacements
%       - a(t) = time-dependent amplitude functions
%
%   Supports both standard Dirichlet conditions (fixed displacements on 
%   surfaces/points) and rigid-body motions (translations + rotations).
%
% INPUT:
%   DIRICHLET    - Struct array with prescribed displacement info per zone.
%                  Each zone may be defined by a surface or point set, and
%                  contains prescribed displacement amplitudes and time laws.
%   DATA         - Global data structure (time steps, file names, etc.).
%   ndim         - Number of spatial dimensions (2D/3D).
%   MESH         - Mesh structure (contains connectivity, nodes, face lists).
%   GEOproperties- Geometric parameters (centers, references for rigid motion).
%   DATALOC      - Local data structure (used e.g. to activate rigid body motion).
%
% OUTPUT:
%   DOFr   - Vector of prescribed DOF indices (assembled across all zones).
%   dR     - Struct containing separated representation of displacement:
%               dR.U  : Block-diagonal matrix with spatial patterns of BCs
%               dR.a  : Time-dependent amplitudes
%               dR.RIGID_BODY_MOTION : (if activated) structure with rigid
%                                       body motion decomposition
%
% NOTES:
%   - Multiple load states (nloads) per zone are supported.
%   - Automatically removes duplicate DOFs across zones.
%   - Ensures that the same DOFs are prescribed across all load states.
%   - If rigid body motion is activated, the prescribed displacements are
%     recomputed in a rotated reference frame and reduced via SVD.
%
% REFERENCES:
%   See DOCS/01_LargeStrainsCode.pdf, page 17
%   Archetypal input files:
%       INPUTS_ROTATION.m   (rigid body motion)
%       INPUTS_translation.m (standard translations)
%
% JAHO - 22-Nov-2020, Commented by ChatGPT-5, on 24th August 2025
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp2.mat')
elseif nargin == 5
    DATALOC = [] ;
end

nzones = length(DIRICHLET);        % Number of zones with prescribed BCs
DOFr = cell(nzones,1);             % Cell array: DOFs per zone
U = cell(nzones,1); a = cell(nzones,1); % Displacement shapes and amplitudes

for izone = 1:nzones
    
    LOCVAR = DIRICHLET(izone);
    LOCVAR = DefaultField(LOCVAR,'IS_RIGID_BODY',0); % Flag for rigid-body BC
    
    if LOCVAR.IS_RIGID_BODY == 1
        % --- Rigid body motion case ---------------------------------------
        % Calls dedicated function that returns rigid-body displacement shape
        % and amplitudes for prescribed motion.
        [Uloc,aloc,DOFrLOCC] = DirichletCONDtimeRIGIDbody(LOCVAR,DATA,ndim,MESH,GEOproperties);
        U{izone} = sparse(Uloc);
        a{izone} = aloc;
        DOFr{izone} = DOFrLOCC;
        
    else
        % --- Standard displacement BCs -----------------------------------
        % BCs can be defined by:
        %   • Surface number (NUMBER_SURFACE)
        %   • Specific mesh points (NUMBER_POINT_MESH)
        
        if isfield(DIRICHLET(izone),'NUMBER_SURFACE') && ~isempty(DIRICHLET(izone).NUMBER_SURFACE)
            MESH = DefaultField(MESH,'NODES_FACES',[]);
            if isempty(MESH.NODES_FACES)
                % Extract nodes on the prescribed surface from mesh file
                rnodLOC = ListOfNodesFACES(DATA.NameFileMeshDATA,DIRICHLET(izone).NUMBER_SURFACE,ndim);
            else
                % Use precomputed node list
                rnodLOC = MESH.NODES_FACES{DIRICHLET(izone).NUMBER_SURFACE};
            end
        elseif ~isempty(DIRICHLET(izone).NUMBER_POINT_MESH)
            % BC defined by explicit node list
            rnodLOC = DIRICHLET(izone).NUMBER_POINT_MESH;
        else
            error('Option not implemented'); % JAHO, 27-Dec-2023
        end
        
        % Candidate DOFs for this zone
        DOFloc = small2large(rnodLOC,ndim);
        
        % --- Loop over load states ---------------------------------------
        nloads = length(DIRICHLET(izone).PRESCRIBED_DISP);
        DOFrLOC = cell(1,nloads); displ = cell(1,nloads);
        for iload = 1:nloads
            PRESCRIBED_DISPLACEMENT = DIRICHLET(izone).PRESCRIBED_DISP(iload);
            for idim = 1:ndim
                if ~isempty(PRESCRIBED_DISPLACEMENT.AMPLITUDE{idim})
                    % DOFs in this direction
                    DOFloc_dim = DOFloc(idim:ndim:end);
                    % Constant prescribed value for all these DOFs
                    dddd = [PRESCRIBED_DISPLACEMENT.AMPLITUDE{idim}] * ones(length(DOFloc_dim),1);
                    % Append
                    displ{iload} = [displ{iload}; dddd];
                    DOFrLOC{iload} = [DOFrLOC{iload}; DOFloc_dim];
                end
            end
            % Check consistency of prescribed DOFs across load states
            if iload > 1
                if ~isempty(setdiff(DOFrLOC{iload},DOFrLOC{iload-1}))
                    error('Ill-defined Boundary Conditions... Load states with distinct DOFs');
                end
            end
        end
        
        % --- Remove repeated DOFs across zones ----------------------------
        if izone > 1
            [REPEATED,iOLD,iNEW] = intersect(cell2mat(DOFr(1:izone-1)),DOFrLOC{1});
            if ~isempty(REPEATED)
                % Keep only non-repeated entries
                IND = setdiff(1:length(DOFrLOC{1}),iNEW);
                DOFrLOC = DOFrLOC{1}(IND);
                for iload = 1:nloads
                    displ{iload} = displ{iload}(IND);
                end
            else
                DOFrLOC = DOFrLOC{1};
            end
        else
            DOFrLOC = DOFrLOC{1};
        end
        
        % --- Assemble space-time separation -------------------------------
        DOFr{izone} = DOFrLOC;
        U{izone} = zeros(length(DOFrLOC),nloads);
        a{izone} = zeros(nloads,length(DATA.STEPS));
        
        for iload = 1:nloads
            % Spatial shape (uniform across nodes of this zone)
            U{izone}(:,iload) = displ{iload};
            
            % Temporal factor from user-provided time law
            FactorSteps = DIRICHLET(izone).PRESCRIBED_DISP(iload).TIMEFUN(DATA.STEPS);
            LIMITS_interval = DIRICHLET(izone).PRESCRIBED_DISP(iload).INTERVAL;
            % Apply activation window
            FactorSteps = FactorSteps .* (DATA.STEPS >= LIMITS_interval(1));
            a{izone}(iload,:) = FactorSteps .* (DATA.STEPS <= LIMITS_interval(2));
        end
        U{izone} = sparse(U{izone});
    end
end

% --- Global assembly -----------------------------------------------------
DOFr = cell2mat(DOFr);          % Concatenate all prescribed DOFs
[DOFr,III] = sort(DOFr);        % Sort DOFs for consistency
dR.U = blkdiag(U{:});           % Block-diagonal assembly of spatial shapes
dR.a = cell2mat(a);             % Amplitudes across all zones
dR.U = dR.U(III,:);             % Reorder to match sorted DOFs

% --- Optional rigid body motion decomposition ----------------------------
DATALOC = DefaultField(DATALOC,'RB_MOTION',[]);
if ~isempty(DATALOC.RB_MOTION)
    % Compute rigid-body motion matrices
    RIGID_BODY_MOTION = RigidBodyMotionTIME(DATA,ndim,MESH,GEOproperties,DATALOC);
    
    % Recompute dR.U and dR.a in rotated reference frame
    D = zeros(size(dR.U,1),size(dR.a,2));
    for itimestep = 1:size(dR.a,2)
        uBAR = dR.U*dR.a(:,itimestep); % Original prescribed disp
        uBAR_rb = RIGID_BODY_MOTION.U(DOFr,:)*RIGID_BODY_MOTION.a(:,itimestep);
        uBARnew = (uBAR - uBAR_rb);
        uBARnew = reshape(uBARnew,ndim,[]);
        uBARnew = RIGID_BODY_MOTION.RotREFP{itimestep}'*uBARnew;
        D(:,itimestep) = uBARnew(:);
    end
    
    disp('Recomputing dR.U and dR.a');
    TOLrel = 1e-10;
    [U,S,V] = RSVDTrel(D,TOLrel);   % Reduced SVD of recomputed displacements
    dR.U = U;
    dR.a = bsxfun(@times,V',S);
else
    RIGID_BODY_MOTION = [];
end

% Attach rigid body motion structure to output
dR.RIGID_BODY_MOTION = RIGID_BODY_MOTION;

%---------------------------------------------------------------------------------
% Before being commented by ChatGPT (24th August 2025)

% function [DOFr,dR] = DirichletCONDtime(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
% % Goal. Determine DOFr and    dR(t)   (prescribed DOFs and value of the corresponding
% % displacements)  in a space-time separated fashion
% %  DOFr  =[DOFrZONE{1} ;  DOFrZONE{2}  ; ...DOFrZONE{nzone}]
% %  dR(t) = U*a(t)
% %  See          DOCS/01_LargeStrainsCode.pdf, page 17
% %
% % JAHO- 22-NOV-2O2O
% %--------------------------------------------------------------------------------------
% if nargin == 0
%     load('tmp.mat')
% elseif nargin == 5
%     DATALOC = [] ;
% end
% 
% nzones = length(DIRICHLET) ; % number of zones with prescribed DOFs
% DOFr = cell(nzones,1); U = cell(nzones,1) ; a = cell(nzones,1) ;
% for  izone = 1:length(DIRICHLET)
%     
%     LOCVAR=  DIRICHLET(izone) ;
%     LOCVAR = DefaultField(LOCVAR,'IS_RIGID_BODY',0) ;
%     
%     if LOCVAR.IS_RIGID_BODY == 1
%         %  See archetypal input file
%         % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/INPUTS_ROTATION.m
%         [Uloc,aloc,DOFrLOCC] =    DirichletCONDtimeRIGIDbody(LOCVAR,DATA,ndim,MESH,GEOproperties) ;
%         U{izone} = sparse(Uloc) ;
%         a{izone} = aloc ;
%         DOFr{izone} = DOFrLOCC ;
%         
%    
%         
%     else
%         
%         % See archetypal input file
%         % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/LARGE_DISPLACEMENTS/RIGID_BODY_MOTION/INPUTS_translation.m
%         
%         %%%%
%         
%         if  isfield(DIRICHLET(izone),'NUMBER_SURFACE') && ~isempty(DIRICHLET(izone).NUMBER_SURFACE)
%             MESH = DefaultField(MESH,'NODES_FACES',[]) ; 
%             if isempty(MESH.NODES_FACES)
%             rnodLOC=  ListOfNodesFACES(DATA.NameFileMeshDATA,DIRICHLET(izone).NUMBER_SURFACE,ndim) ;
%             else
%                 rnodLOC = MESH.NODES_FACES{DIRICHLET(izone).NUMBER_SURFACE} ; % 14-Apr-2025
%             end
%         elseif ~isempty(DIRICHLET(izone).NUMBER_POINT_MESH)
%             rnodLOC = DIRICHLET(izone).NUMBER_POINT_MESH ;
%         else
%             error('Option not implemented')  % JAHO, 27-Dec-2023
%         end
%         DOFloc = small2large(rnodLOC,ndim) ;  % Candidates for being DOFr
%         
%         % Determining the DOFr corresponding to this zone
%         % Notice that the DOFr must be the same for all load states
%         % ----------------------------------------------------------
%         DOFrLOC =  cell(1,length(DIRICHLET(izone).PRESCRIBED_DISP)) ;
%         displ =  cell(1,length(DIRICHLET(izone).PRESCRIBED_DISP)) ;
%         nloads=  length(DIRICHLET(izone).PRESCRIBED_DISP) ;
%         for iload = 1:nloads
%             PRESCRIBED_DISPLACEMENT = DIRICHLET(izone).PRESCRIBED_DISP(iload);
%             for idim = 1:ndim
%                 if ~isempty(PRESCRIBED_DISPLACEMENT.AMPLITUDE{idim})
%                     DOFloc_dim = DOFloc(idim:ndim:end) ;
%                     dddd =[ PRESCRIBED_DISPLACEMENT.AMPLITUDE{idim} ]*ones(length(DOFloc_dim),1) ;
%                     displ{iload} = [ displ{iload}; dddd] ;
%                     DOFrLOC{iload} = [DOFrLOC{iload};DOFloc_dim] ;  % DOFr corresponding
%                 end
%             end
%             if iload >1
%                 % Checking that the DOFs are the same
%                 if ~isempty(setdiff(DOFrLOC{iload},DOFrLOC{iload-1}))
%                     error('Ill-defined Boundary Conditions... Load states with distinct DOFs')
%                 end
%             end
%         end
%         
%         % Checking that there are no repeated indexes  (and remove them if indeed there are repetitions)
%         if izone >1
%             [REPEATED,iOLD,iNEW] = intersect(cell2mat(DOFr(1:izone-1)),DOFrLOC{1}) ;
%             if ~isempty(REPEATED)
%                 IND =  setdiff(1:length(DOFrLOC{1}),iNEW) ;
%                 DOFrLOC = DOFrLOC{1}(IND) ;
%                 for iload =1:nloads
%                     displ{iload} = displ{iload}(IND) ;
%                 end
%             else
%                 DOFrLOC = DOFrLOC{1} ;
%             end
%         else
%             DOFrLOC = DOFrLOC{1} ;
%         end
%         %
%         DOFr{izone} = DOFrLOC;
%         U{izone} = zeros(length(DOFrLOC),nloads) ;
%         a{izone} = zeros(nloads,length(DATA.STEPS)) ;
%         
%         
%         for  iload = 1:nloads
%             U{izone}(:,iload) = displ{iload} ;    % Uniform Boundary Conditions
%             FactorSteps =  DIRICHLET(izone).PRESCRIBED_DISP(iload).TIMEFUN(DATA.STEPS) ;
%             LIMITS_interval =  DIRICHLET(izone).PRESCRIBED_DISP(iload).INTERVAL ;
%             FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
%             a{izone}(iload,:) = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
%         end
%         U{izone} = sparse( U{izone}) ;
%         
%     end
%     
%     
% end
% 
% %dR.U = U ;
% %dR.a = a ;
% 
% DOFr  =cell2mat(DOFr) ;
% [DOFr,III ]= sort(DOFr);
% 
% dR.U = blkdiag(U{:}) ;
% dR.a = cell2mat(a) ;
% 
% dR.U = dR.U(III,:) ;
% 
% 
% %%% RIGID BODY MOTION
% % -----------------------------
% % Compute, for each time step
% % 1) Translation  (3xtime_steps)
% % 2) Relative coordinates
% % 3) Rotation (3x3xtime_steps)
% 
% DATALOC = DefaultField(DATALOC,'RB_MOTION',[]) ;
% 
% if ~isempty(DATALOC.RB_MOTION)
%     RIGID_BODY_MOTION = RigidBodyMotionTIME(DATA,ndim,MESH,GEOproperties,DATALOC)   ;
%     
%     % We have to redefine dR.U and dR.a ...
%     D = zeros(size(dR.U,1),size(dR.a,2)) ;
%     for itimestep = 1:size(dR.a,2)
%         uBAR = dR.U*dR.a(:,itimestep) ; % Original prescribed displacement
%         uBAR_rb  = RIGID_BODY_MOTION.U(DOFr,:)*RIGID_BODY_MOTION.a(:,itimestep) ; % Rigid body motion
%         uBARnew = (uBAR-uBAR_rb) ; %
%         uBARnew = reshape(uBARnew,ndim,[]) ;
%         uBARnew = RIGID_BODY_MOTION.RotREFP{itimestep}'*uBARnew ;
%         D(:,itimestep) =uBARnew(:);
%         
%         %
%     end
%     disp('Recomputing dR.U and dR.a')
%     TOLrel = 1e-10;
%     [U,S,V] = RSVDTrel(D,TOLrel) ;
%     dR.U = U ;
%     a = bsxfun(@times,V',S) ;
%     dR.a = a;
%     
%     
% else
%     RIGID_BODY_MOTION = [] ;
% end
% 
% dR.RIGID_BODY_MOTION = RIGID_BODY_MOTION;
% 
