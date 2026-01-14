function [shapef,dershapef,COORnodes] = ComputeElementShapeFun_givenpoint(TypeElement,nnodeE,TypeIntegrand,xi) ;
% This function returns, for each "TypeElement" and 'TypeIntegrand'
% (K/RHS)*
% weig = Vector of Gauss weights (1xngaus)
% posgp: Position of Gauss points  (ndim x ngaus)
% shapef: Array of shape functions (ngaus x nnodeE)
% dershape: Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)
%
% *) TypeIntegrand = 'K': For integration of Conductance Matrix
% *) TypeIntegrand = 'RHS': For integration of flux vectors

switch TypeElement
    case 'Linear'
        if nnodeE ==2
             error('Option not implemented')
            [shapef,dershapef] = Linear2NInPoints_givenpoint(TypeIntegrand,xi) ;
        elseif   nnodeE ==3 
             error('Option not implemented')
             [shapef,dershapef] = Linear3NInPoints_givenpoint(TypeIntegrand,xi) ;
        else
            error('Option not implemented')
        end
    case 'Triangle'
        if nnodeE ==3
             error('Option not implemented')
            [shapef,dershapef] = Triangular3NInPoints_givenpoint(TypeIntegrand,xi) ;
        else
            error('Option not implemented')
        end
        
    case 'Quadrilateral'
        if nnodeE ==4
             
            [shapef,dershapef,COORnodes] = Quadrilateral4NInPoints_givenpoint(xi) ;
        elseif nnodeE ==9            
            [shapef,dershapef,COORnodes] = Quadrilateral9NInPoints_givenPoints(xi) ;
        else 
              [shapef,dershapef] = QuadrilateralGen_givenPoints(nnodeE,xi) ;
        end
    case 'Hexahedra'
        if nnodeE ==8
             error('Option not implemented')
            [shapef,dershapef] = Hexahedra8NInPoints_givenpoint(xi) ;
        elseif nnodeE ==27
            [shapef,dershapef] = Hexahedra27NInPoints_givenpoint(xi)  ;
        else
            error('Option not implemented')
        end
    otherwise
        error('Option not implemented')
end