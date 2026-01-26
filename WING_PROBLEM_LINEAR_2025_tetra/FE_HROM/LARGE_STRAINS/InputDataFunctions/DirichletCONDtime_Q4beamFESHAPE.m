function [DOFr,dR] = DirichletCONDtime_Q4beamFESHAPE(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
% Goal. Determine DOFr and    dR(t) .
% Case boundary conditions Q4-beam 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/18_BEAMQ4.mlx
 %
% JAHO- 20-DEC-2O23, Balmes, Barcelona 
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
      %  DATA.EnableInverseMappingBoundaryConditions = 0; 

elseif nargin == 5
    DATALOC = [] ;
end


% BOUNDARY NODES THAT WILL BE SUBJECT TO
% THE FE-SHAPE-LIKE BOUNDARY CONDITIONS
MESHcoar  = DATALOC.MESHcoarse_FOR_DIRICHLET_BOUND_COND ; % COARSE MESH, CONTAINING THE CORNER NODES
FACES_BND  = DATALOC.LabelEntitiesDefiningBoundary ; % FACES DEFINING THE BOUNDARIES OF THE FINE MESH DOMAIN
NODES_FACES = MESH.NODES_FACES(FACES_BND) ;  %
NODES_FACES = NODES_FACES(:) ;               %
rnodLOC = unique(cell2mat(NODES_FACES)) ;  % LIST OF NODES WITH PRESCRIBED DISPLACEMENTS (BOUNDARY)
% ----------------------------------------------------------
COORbnd = MESH.COOR(rnodLOC,:) ;  % COORDINATES
DATA = DefaultField(DATA,'EnableInverseMappingBoundaryConditions',0); % =

% SHAPE FUNCTIONS
% ---------------
% See example in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/Fun3D/Dpoly3D_P5.m
%
switch  MESHcoar.TypeElement
    case {'Quadrilateral','Hexahedra'}
        
        [nnodeBND, ndim] = size(MESHcoar.COOR) ;  % Number of nodes of the mesh upon which displacements are imposed
        % Number of dimensions
        nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
        
        if  ceil(nnodes_DIM) == nnodes_DIM
            
            if nnodes_DIM == 2 || DATA.EnableInverseMappingBoundaryConditions ==0
                ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
                DATAshape = ShapeFunCoefficients(MESHcoar.COOR,ORDER_POLYNOMIALS) ;
                DATAlocSHAPE.DATAshape  = DATAshape;
                xLIM = [] ;
                DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
                [Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;
            else
                 Nshape = InverseMapping3Dinterpolation(MESHcoar,DATA,COORbnd) ; 
                
               
            end
        else
            % Serendipity element
            ORDER_POLYNOMIALS = (ceil(nnodes_DIM)-1)*ones(1,ndim) ;
            DATAshape = ShapeFunCoefficientsSEREN(MESHcoar.COOR,ORDER_POLYNOMIALS) ;
            DATAlocSHAPE.DATAshape  = DATAshape;
            xLIM = [] ;
            DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
            [Nshape, ~,~ ]=    ShapeFunctionFEseren(xLIM,COORbnd,DATAlocSHAPE) ;
        end
        
        
        
        
    case {'Triangle'}
        
        [nnodeBND, ndim] = size(MESHcoar.COOR) ;
        % Order of polynomial
        % nnodes_DIM = (nnodeBND)^(1/ndim)  ;  % Number of nodes per dimension
        % ORDER_POLYNOMIALS = (nnodes_DIM-1)*ones(1,ndim) ;
        DATAshape = ShapeFunCoefficientsTRI(MESHcoar.COOR) ;
        DATAlocSHAPE.DATAshape  = DATAshape;
        xLIM = [] ;
        DATAlocSHAPE.ORDER_POLYNOMIALS = [];
        [Nshape, ~,~ ]=    ShapeFunctionFEtri(xLIM,COORbnd,DATAlocSHAPE) ;
        
    otherwise
        error('Element not implemented yet')
end

%


% Displacement pattern
DOFr = small2large(rnodLOC,ndim) ;  % Prescribed DOFS
DIRICHLET.AMPLITUDE = DIRICHLET.AMPLITUDE(:); 
U = zeros(length(DOFr),1) ;   % Space pattern
for idim = 1:ndim
    U(idim:ndim:end) = Nshape*DIRICHLET.AMPLITUDE(idim:ndim:end) ;
end
a = DIRICHLET.TIMEFUN(DATA.STEPS) ;  % Time pattern

dR.U = U ;
dR.a = a ;
%
%      nloads=  1 ;
% %---------------------------
% % Decomposition space.time
% U = {zeros(length(DOFrLOC),nloads)} ; % Nodal amplitude ---space pattern
% a = {zeros(nloads,length(DATA.STEPS))} ; % Time factor  --- time functions
%
%
% for  iload = 1:nloads
%     U{izone}(:,iload) = displ{iload} ;    % Uniform Boundary Conditions
%     FactorSteps =  DIRICHLET(izone).PRESCRIBED_DISP(iload).TIMEFUN(DATA.STEPS) ;
%     LIMITS_interval =  DIRICHLET(izone).PRESCRIBED_DISP(iload).INTERVAL ;
%     FactorSteps = FactorSteps.*(DATA.STEPS >= LIMITS_interval(1)) ;
%     a{izone}(iload,:) = FactorSteps.*(DATA.STEPS <= LIMITS_interval(2)) ;
% end
% U{izone} = sparse( U{izone}) ;

% end

%
% end

%dR.U = U ;
%dR.a = a ;

% DOFr  =cell2mat(DOFr) ;
% [DOFr,III ]= sort(DOFr);
%
% dR.U = blkdiag(U{:}) ;
% dR.a = cell2mat(a) ;
%
% dR.U = dR.U(III,:) ;






% % Determining the DOFr corresponding to this zone
% % Notice that the DOFr must be the same for all load states
% % ----------------------------------------------------------
% DOFrLOC =  cell(1,length(DIRICHLET(izone).AMPLITUDE)) ;
% displ =  cell(1,length(DIRICHLET(izone).AMPLITUDE)) ;
% nloads=  length(DIRICHLET(izone).AMPLITUDE) ;
% for iload = 1:nloads
%     PRESCRIBED_DISPLACEMENT = DIRICHLET(izone).PRESCRIBED_DISP(iload);
%     for idim = 1:ndim
%         if ~isempty(PRESCRIBED_DISPLACEMENT.AMPLITUDE{idim})
%             DOFloc_dim = DOFloc(idim:ndim:end) ;
%             dddd =[ PRESCRIBED_DISPLACEMENT.AMPLITUDE{idim} ]*ones(length(DOFloc_dim),1) ;
%             displ{iload} = [ displ{iload}; dddd] ;
%             DOFrLOC{iload} = [DOFrLOC{iload};DOFloc_dim] ;  % DOFr corresponding
%         end
%     end
%     if iload >1
%         % Checking that the DOFs are the same
%         if ~isempty(setdiff(DOFrLOC{iload},DOFrLOC{iload-1}))
%             error('Ill-defined Boundary Conditions... Load states with distinct DOFs')
%         end
%     end
% end

% % Checking that there are no repeated indexes  (and remove them if indeed there are repetitions)
% if izone >1
%     [REPEATED,iOLD,iNEW] = intersect(cell2mat(DOFr(1:izone-1)),DOFrLOC{1}) ;
%     if ~isempty(REPEATED)
%         IND =  setdiff(1:length(DOFrLOC{1}),iNEW) ;
%         DOFrLOC = DOFrLOC{1}(IND) ;
%         for iload =1:nloads
%             displ{iload} = displ{iload}(IND) ;
%         end
%     else
%         DOFrLOC = DOFrLOC{1} ;
%     end
% else
%     DOFrLOC = DOFrLOC{1} ;
% end
%

%%% RIGID BODY MOTION
% -----------------------------
% Compute, for each time step
% 1) Translation  (3xtime_steps)
% 2) Relative coordinates
% 3) Rotation (3x3xtime_steps)

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

dR.RIGID_BODY_MOTION = [];

