function [DOFr,dR,DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = DirichletCONDtime_FESHAPEperiod(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC)
% Goal. Determine DOFr and    dR(t) .
% Case boundary conditions expressed as standard boundary shape functions +
% periodic boundaryb conditions
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
%
% JAHO- 17-nov-2O23, Green's Aribau, Barcelona
%--------------------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
INFO_PERIODIC_CONDITIONS = [] ;

% BOUNDARY NODES THAT WILL BE SUBJECT TO
% THE FE-SHAPE-LIKE BOUNDARY CONDITIONS
MESHcoar  = DATALOC.MESHcoarse_FOR_DIRICHLET_BOUND_COND ; % COARSE MESH, CONTAINING THE CORNER NODES
DATALOC  = DefaultField(DATALOC,'LabelEntitiesDefiningBoundary',[1,2,3,4]) ;
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
USE_PERIODIC = 1;
if USE_PERIODIC == 0
    % JUST FOR TESTING PURPOSE
    
    DOFr = small2large(rnodLOC,ndim) ;  % Prescribed DOFS
    
    U = zeros(length(DOFr),1) ;   % Space pattern
    for idim = 1:ndim
        U(idim:ndim:end) = Nshape*DIRICHLET.AMPLITUDE(idim:ndim:end) ;
    end
    a = DIRICHLET.TIMEFUN(DATA.STEPS) ;  % Time pattern
    
    dR.U = U ;
    dR.a = a ;
    
    dR.RIGID_BODY_MOTION = [];
    
else
    %     DISP_CONDITIONS.A = Abub ; %
    % DISP_CONDITIONS.DOFr = bDOFS(hDOFS) ; % SLAVE DOFS
    % DISP_CONDITIONS.DOFm = bDOFS(gDOFS) ; % MASTER DOFS
    % DISP_CONDITIONS.DOFl = lDOFS ; % remaining DOFS
    DOFr = [] ; dR= [] ;
    DATALOC = DefaultField(DATALOC,'TypeFunctionDisplacementInterfaces','QUADRILATERAL_LINEAR') ;
    switch DATALOC.TypeFunctionDisplacementInterfaces
        case 'QUADRILATERAL_LINEAR'
            [ DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = PeriodicQ4(DIRICHLET,DATA,ndim,MESH,GEOproperties,DATALOC,Nshape,rnodLOC,MESHcoar);
        otherwise
            error('Option not implemented')
    end
    
end

