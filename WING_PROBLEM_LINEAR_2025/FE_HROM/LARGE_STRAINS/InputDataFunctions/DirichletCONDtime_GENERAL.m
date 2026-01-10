function [DOFr,dR,OTHEROUTPUT] = DirichletCONDtime_GENERAL(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS)
%--------------------------------------------------------------------------
%  DirichletCONDtime_GENERAL
%
%  General interface for applying **Dirichlet boundary conditions** (prescribed
%  displacements) in the context of finite element problems, including those
%  using EIFEM with co-rotational or homogenization-based features.
%
%  This function dispatches the appropriate boundary condition handler based on
%  the geometry, training setup, and special constraints such as:
%    - FE shape function-based displacement modes
%    - Plate kinematics from hexahedral elements
%    - Periodicity constraints (1D, full face)
%    - Rigid body motion and warping constraints
%
%  INPUTS:
%    - DIRICHLET       : structure array specifying which surfaces/nodes are constrained
%    - DATA            : problem definition, including RVE flags, time parameters, etc.
%    - ndim            : spatial dimension (2 or 3)
%    - MESH            : geometry, connectivities, etc.
%    - GEOproperties   : computed geometric data such as centroids, normals, etc.
%    - OTHER_INPUTS    : optional structure including:
%         > MESHcoarse_FOR_DIRICHLET_BOUND_COND
%         > TypeFunctionDisplacementTRAINING
%         > RVE_like_boundary_conditions
%         > BOUNDARY_CONDITIONS_SHAPE_FUNCTIONS_ONLY_SELECTED_SURFACES
%
%  OUTPUTS:
%    - DOFr        : global DOFs constrained by Dirichlet conditions
%    - dR          : structure with imposed displacements (can include time functions)
%    - OTHEROUTPUT : optionally includes:
%         > DISP_CONDITIONS: substructure for periodic or rigid setups
%         > INFO_PERIODIC_CONDITIONS: bookkeeping for periodic BCs
%
%  SUPPORTED OPTIONS (via OTHER_INPUTS.TypeFunctionDisplacementTRAINING):
%    - '':                          standard shape function-based boundary modes
%    - 'FE_SHAPE_functions_only_selected_surfaces'
%    - 'Q4beam':                    specific for beam-type training elements
%    - 'PLATE_ELEMENT_FROM_LINEAR_HEXAHEDRA_ELEMENT'
%    - 'MIXED_PLATE_ELEMENT_FROM_LINEAR_HEXAHEDRA_ELEMENT'
%    - 'PERIODICITY_*':            various periodic BCs (full, 1D, with/without corners)
%    - 'RIGID_BODY_WARPING_ON_SURFACES'
%    - 'RIGID_BODY_ON_SURFACES':   with automatic detection of large/small rotations
%    - or leave empty for `DirichletCONDtime` (standard)
%
%  SPECIAL CASES:
%    - If `RVE_like_boundary_conditions == 1`, Dirichlet values are derived from
%      a macro deformation tensor (`DIRICHLET.MACRODEF`) — see `DirichletCONDtime_homogZERO`.
%
%  REMARKS:
%    - Automatically detects if large or small rotations should be applied in rigid modes
%    - Modular design to allow plug-in of new constraint types without modifying this dispatcher
%
%  AUTHOR:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 29-Jan-2023 to 28-Oct-2024
%    Comments by ChatGPT4, 13-May-2025
%  SEE ALSO:
%    - DirichletCONDtime_FESHAPE, DirichletCONDtime_FESHAPEsurf, DirichletCONDtime_rigid
%    - DirichletCONDtime_Period1D, DirichletCONDtime_homogZERO
%    - Corresponding documentation in README_HOMOG.mlx, README_PLATES.mlx, etc.
%
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
    OTHER_INPUTS.TypeFunctionDisplacementTRAINING = ''  ;
end
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'TypeFunctionDisplacementTRAINING','') ;  %
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'BOUNDARY_CONDITIONS_SHAPE_FUNCTIONS_ONLY_SELECTED_SURFACES',[]) ;  %
OTHER_INPUTS = DefaultField(OTHER_INPUTS,'MESHcoarse_FOR_DIRICHLET_BOUND_COND',[]) ;  %



% 'PLATE_ELEMENT_FROM_LINEAR_HEXAHEDRA_ELEMENT'
OTHEROUTPUT = [] ;
OTHEROUTPUT.INFO_PERIODIC_CONDITIONS = [] ;
if ~isempty(OTHER_INPUTS.MESHcoarse_FOR_DIRICHLET_BOUND_COND)
    % --------------------------------------------------------------------------------------------------------------------------
    switch OTHER_INPUTS.TypeFunctionDisplacementTRAINING
        case ''
            % BOUNDARY CONDITIONS WITH PRESCRIBED DISPLACEMENTS EQUAL TO STANDARD
            % SHAPE FUNCTIONS
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
            % JAHO, 29-Jan-2023
            [DOFr,dR] = DirichletCONDtime_FESHAPE(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
        case 'FE_SHAPE_functions_only_selected_surfaces'
            
            [DOFr,dR] = DirichletCONDtime_FESHAPEsurf(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
        case 'Q4beam'
            [DOFr,dR] = DirichletCONDtime_Q4beamFESHAPE(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
        case {'PERIODICITY_ALL_FACES_WITH_CORNERS','PERIODICITY_ALL_FACES_WITH_NOCORNERS'}
            [DOFr,dR,DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = ...
                DirichletCONDtime_FESHAPEperiod(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
            OTHEROUTPUT.DISP_CONDITIONS = DISP_CONDITIONS;
            OTHEROUTPUT.INFO_PERIODIC_CONDITIONS = INFO_PERIODIC_CONDITIONS ;
            
            
        case 'RIGID_BODY_WARPING_ON_SURFACES'
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/01_BEAMQ8_ELAST.mlx
            [DOFr,dR] = DirichletCONDtime_rigid_WARPING(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
        otherwise
            error('Option not implemented')
    end
else
    switch OTHER_INPUTS.TypeFunctionDisplacementTRAINING
        case 'Periodicity_Quadrilateral_With_Corners_Homogenization'
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/13_HOMOG.mlx
            % 7-Oct-2025
            
           [DOFr,dR,DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = ...
                DirichletCONDtime_PeriodHomogQ4(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
            OTHEROUTPUT.DISP_CONDITIONS = DISP_CONDITIONS;
            OTHEROUTPUT.INFO_PERIODIC_CONDITIONS = INFO_PERIODIC_CONDITIONS ;
        
        
        case 'PLATE_ELEMENT_FROM_LINEAR_HEXAHEDRA_ELEMENT'
            % JAHO, 16-feb-2023
            % PLATE-LIKE BOUNDARY CONDITIONS
            % sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/02_PLATES/README_PLATES.mlx
            [DOFr,dR] = DirichletCONDtime_PLATESfromFEshape(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
        case 'MIXED_PLATE_ELEMENT_FROM_LINEAR_HEXAHEDRA_ELEMENT'
            % JAHO, 22-Aug-2023
            % MIXED TRACTION/DISPLACEMENTS BOUNDARY CONDITIONS (LINEAR)
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/103_EIFEM_PLATES/02_CONV_SIZE_TRAC.mlx
            [DOFr,dR,OTHEROUTPUT] = MIXEDCONDtime_PLATESfromFEshape(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
        case 'RIGID_BODY_WARPING_ON_SURFACES'
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/01_BEAMQ8_ELAST.mlx
            [DOFr,dR] = DirichletCONDtime_rigid_WARPING(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
            %              case 'FE_SHAPE_functions_only_selected_surfaces'
            %
            %               [DOFr,dR] = DirichletCONDtime_FESHAPEsurf(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
         case 'RIGID_BODY_ON_SURFACES'
            % The snippet below was introduced to reconcile large/small
            % rotations implementation (confusion diagnosed in 
            % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/10_AUXETIC_2D.mlx
             if isfield(DIRICHLET(1).PRESCRIBED_DISP,'ROTATION')
                 % Large rotations
                 [DOFr,dR] = DirichletCONDtime_rigid(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
             else        
                 % Small rotations
                 [DOFr,dR] = DirichletCONDtime_rigid2(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
             end
           
            
        case 'PERIODICITY_1D'
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/01_meta_1D.mlx
            [DOFr,dR,DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = ...
                DirichletCONDtime_Period1D(DIRICHLET,DATA,ndim,MESH,OTHER_INPUTS) ;
            OTHEROUTPUT.DISP_CONDITIONS = DISP_CONDITIONS;
            OTHEROUTPUT.INFO_PERIODIC_CONDITIONS = INFO_PERIODIC_CONDITIONS ;
        case 'PERIODICITY_1D_PLUS_ZERO_CONSTRAINTS'
             [DOFr,dR,DISP_CONDITIONS,INFO_PERIODIC_CONDITIONS] = ...
                DirichletCONDtime_Period1Dpz(DIRICHLET,DATA,ndim,MESH,OTHER_INPUTS) ;
            OTHEROUTPUT.DISP_CONDITIONS = DISP_CONDITIONS;
            OTHEROUTPUT.INFO_PERIODIC_CONDITIONS = INFO_PERIODIC_CONDITIONS ;
        otherwise
            if OTHER_INPUTS.RVE_like_boundary_conditions == 1
                % Boundary with conditions defined by the disp. gradient
                % DIRICHLET.MACRODEF
                [DOFr,dR] = DirichletCONDtime_homogZERO(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
            else
                % Standard
                [DOFr,dR] = DirichletCONDtime(DIRICHLET,DATA,ndim,MESH,GEOproperties,OTHER_INPUTS) ;
            end
    end
end