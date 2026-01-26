function [weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
%%
% =========================================================================
% ComputeElementShapeFun — Element shape functions, Gauss points & derivatives
% =========================================================================
% PURPOSE
%   Return Gauss weights/points, shape functions, and their derivatives for a
%   given finite element type, number of nodes per element, and integration
%   purpose (stiffness vs RHS).
%
% SIGNATURE
%   [weig, posgp, shapef, dershapef] = ComputeElementShapeFun(TypeElement, nnodeE, TypeIntegrand)
%
% INPUTS
%   TypeElement    (char) : 'Linear' | 'Triangle' | 'Quadrilateral' | 'Tetrahedra' | 'Hexahedra'
%   nnodeE         (int)  : # nodes per element (depends on TypeElement)
%   TypeIntegrand  (char) : 'K'   → stiffness/stress integration
%                           'RHS' → loads/mass (right-hand side) integration
%
% OUTPUTS
%   weig        : 1×ngaus vector of Gauss weights
%   posgp       : ndim×ngaus Gauss-point coordinates in parent domain
%   shapef      : ngaus×nnodeE shape-function values at Gauss points
%   dershapef   : ndim×nnodeE×ngaus derivatives of shape functions
%
% WHAT THE FUNCTION DOES
%   • Dispatches to the appropriate element-specific quadrature/shape routine:
%       - 'Linear'        : Linear2NInPoints, Linear3NInPoints
%       - 'Triangle'      : Triangular3NInPoints
%       - 'Quadrilateral' : Quadrilateral4NInPoints, Quadrilateral8NInPoints, Quadrilateral9NInPoints
%       - 'Tetrahedra'    : Tetrahedra4NInPoints
%       - 'Hexahedra'     : Hexahedra8NInPoints, Hexahedra27NInPoints
%   • Each routine returns (weig, posgp, shapef, dershapef) tailored to the
%     requested integration type ('K' or 'RHS').
%
% SUPPORTED ELEMENTS & NODE COUNTS
%   Linear(2,3), Triangle(3), Quadrilateral(4,8,9),
%   Tetrahedra(4), Hexahedra(8,27).
%
% PRACTICAL NOTES
%   • Use 'K' when assembling stiffness/stress-related terms; use 'RHS' for
%     mass or external-force integrals (may employ different quadrature).
%   • The returned dershapef is with respect to the parent element coordinates.
%     Mapping to physical space requires Jacobians built elsewhere.
%   • Unsupported (TypeElement, nnodeE) pairs raise an error.
%
% .mlx references: (none)
%
% Written by Joaquín A. Hernández (JAHO), UPC/CIMNE
% Contact: jhortega@cimne.upc.edu
% Automatically commented by ChatGPT — 07-Nov-2025
% =========================================================================

switch TypeElement
    case 'Linear'
        if nnodeE ==2
            [weig,posgp,shapef,dershapef] = Linear2NInPoints(TypeIntegrand) ;
        elseif   nnodeE ==3 
             [weig,posgp,shapef,dershapef] = Linear3NInPoints(TypeIntegrand) ;
        else
            error('Option not implemented')
        end
    case 'Triangle'
        if nnodeE ==3
            [weig,posgp,shapef,dershapef] = Triangular3NInPoints(TypeIntegrand) ;
        else
            error('Option not implemented')
        end
        
    case 'Quadrilateral'
        if nnodeE ==4
            [weig,posgp,shapef,dershapef] = Quadrilateral4NInPoints(TypeIntegrand) ;
        elseif nnodeE ==9 
            [weig,posgp,shapef,dershapef] = Quadrilateral9NInPoints(TypeIntegrand) ;
        elseif nnodeE == 8
            [weig,posgp,shapef,dershapef] = Quadrilateral8NInPoints(TypeIntegrand) ;
        else
            error('Option not implemented')
        end
    case 'Tetrahedra'
        if nnodeE ==4 
             [weig,posgp,shapef,dershapef] = Tetrahedra4NInPoints(TypeIntegrand) ; 
        else
             error('Option not implemented')
        end    
        
        
    case 'Hexahedra'
        if nnodeE ==8
            [weig,posgp,shapef,dershapef] = Hexahedra8NInPoints(TypeIntegrand) ;
        elseif nnodeE ==27
            [weig,posgp,shapef,dershapef] = Hexahedra27NInPoints(TypeIntegrand) ;
        else
            error('Option not implemented')
        end
    otherwise
        error('Option not implemented')
end