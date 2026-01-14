function  [Vall,MESH,DATA,INFOPLOTMESHBND,FluctuationFacesGLO] = FictInterfaceDISP_raw(MESH,DATAcommon,DATA,Mintf,DATAoffline)
%--------------------------------------------------------------------------
% FUNCTION: FictInterfaceDISP_raw
%
% PURPOSE:
%   Returns the interface displacement modes `Vall` for a given subdomain mesh, based
%   on the type of interface specified in `DATAcommon.TypeFunctionDisplacementInterfaces`.
%   Unlike `FictInterfaceDISP.m`, this version **does not compute the fluctuation modes**.
%
% DESCRIPTION:
%   This function dispatches the construction of kinematic interface modes (e.g.,
%   rigid body, linear, or higher-order interpolated displacements) depending on the
%   type of boundary/interface geometry. These interface modes are crucial in multiscale
%   or domain decomposition methods (e.g., EIFEM), where reduced-order bases or coupling
%   strategies require accurate interfacial representations.
%
% INPUTS:
%   - MESH        : Finite element mesh structure of the subdomain.
%   - DATAcommon  : Structure with high-level modeling choices and options,
%                   including field `.TypeFunctionDisplacementInterfaces` that controls logic.
%   - DATA        : Structure containing physical and numerical modeling data.
%   - Mintf       : (Optional) Interface mass matrix or relevant geometric structure (can be empty).
%   - DATAoffline : Structure with flags and settings for offline computations.
%
% OUTPUTS:
%   - Vall                : Matrix whose columns are the prescribed displacement modes
%                           on the interface (rigid body, polynomial, user-defined, etc.).
%   - MESH                : Possibly updated mesh structure.
%   - DATA                : Possibly updated data structure.
%   - INFOPLOTMESHBND     : Optional structure with boundary visualization info (only for some modes).
%   - FluctuationFacesGLO : Empty here (in contrast to full version), placeholder for interface fluctuation modes.
%
% SUPPORTED TYPES in `DATAcommon.TypeFunctionDisplacementInterfaces`:
%   - 'GIVEN_BY_USER'                     : Load modes defined manually in `FictInterfaceModes_USER`.
%   - 'QUADRILATERAL_LINEAR'             : Linear interface displacements (2D quadrilateral elements).
%   - 'QUADRILATERAL_LINEAR_UNCOUPLED'   : Same as above, but uncoupled mode computation.
%   - 'QUADRILATERAL_QUADRATIC'          : Higher-order quadratic interface displacements.
%   - 'PLATE_QUAD_LINEAR'                : Plate bending kinematics for shell-like interfaces.
%   - 'TRIANGULAR_LINEAR'                : 2D triangular linear boundary elements.
%   - '1D_LINEAR'                        : 1D interface (beam or bar elements).
%   - 'HEXAHEDRA_LINEAR' / 'HEXAHEDRA_8nodes'        : 3D brick elements, linear.
%   - 'HEXAHEDRA_QUADRATIC' / 'HEXAHEDRA_27nodes'    : 3D brick, quadratic.
%   - 'HEXAHEDRA_20nodes'                : 3D brick, serendipity formulation.
%
% NOTE:
%   - All switch options call specialized helper functions like `FictInterfaceDISP_QUADlin_ONLY_GID`,
%     `FictInterfaceDISP_HEXAquad_onlyGID`, etc., depending on element type.
%   - INFOPLOTMESHBND is only non-empty for some visualization-enabled types.
%
% SEE ALSO:
%   FictInterfaceDISP, FictInterfaceDISP_QUADlin_ONLY_GID, FictInterfaceModes_USER,
%   FictInterfaceDISP_HEXAquad_onlyGID, GeometricVarDOMAINS
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (CIMNE/UPC)
%   Last updated: 12-Feb-2023, Barcelona
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------



% Similar to  FictInterfaceDISP.m, but without computing fluctuations
% JAHO 12-FEB-2023
%5------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
FluctuationFacesGLO =[] ;
 

INFOPLOTMESHBND = [] ; 
FluctuationFacesGLO = [] ; 
% ---------------------------------------------------------------------------------------------
DATAcommon = DefaultField(DATAcommon,'TypeFunctionDisplacementInterfaces','QUADRILATERAL_LINEAR') ;




Mintf = [] ; 
FluctuationFacesGLO = [] ; 
 switch DATAcommon.TypeFunctionDisplacementInterfaces
     case 'GIVEN_BY_USER'
     % Created on 2-April-2025, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/08_GENERAL_MODES.mlx
     [Vall,DATA] = FictInterfaceModes_USER(MESH,DATAcommon,DATA,DATAoffline);
     
    case 'QUADRILATERAL_LINEAR'
    [Vall,Mintf,MESH,DATA,FluctuationFacesGLO] = FictInterfaceDISP_QUADlin_ONLY_GID([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ; 
    
        case '1D_LINEAR'
            % See 
    [Vall,Mintf,MESH,DATA,FluctuationFacesGLO] = FictInterfaceDISP_1D([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ; 
    
    case 'QUADRILATERAL_LINEAR_UNCOUPLED'
    [Vall,MESH,DATA] = FictInterfaceDISP_QUADlin_UNCOUPLED(MESH,DATAcommon,DATA,Mintf,DATAoffline) ; 
    
    
    
    
     case 'QUADRILATERAL_QUADRATIC'
 %   [Vall,Mintf,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP_QUADq([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ; 
      [Vall,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP_QUADquad_onlyGID(MESH,DATAcommon,DATA,DATAoffline) ; 
    case 'PLATE_QUAD_LINEAR'
     [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_PLATEquadLIN([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ; 
      case 'TRIANGULAR_LINEAR'
        [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_TRIlin([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ; 
       case {'HEXAHEDRA_LINEAR','HEXAHEDRA_8nodes'}  % Alternative names
    [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_HEXAlin([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ;    
        
 case {'HEXAHEDRA_QUADRATIC','HEXAHEDRA_27nodes'}  % Alternative names
   %  disp('Rather use ')
   
 %   [Vall,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP_HEXAquad([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ;   
     [Vall,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP_HEXAquad_onlyGID([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ;    
    
    
    case 'HEXAHEDRA_20nodes'
    [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_HEXA20NODES([],[],MESH,DATAcommon,DATA,Mintf,[],[],DATAoffline) ;    

    
    
    otherwise
        error('Option not implemented yet')
end

