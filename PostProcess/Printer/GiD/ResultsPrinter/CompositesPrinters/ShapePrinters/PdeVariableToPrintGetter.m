classdef PdeVariableToPrintGetter < handle
       
    properties (Access = private)
      physicalProblem  
    end
    
    methods (Access = public)
        
        function obj = PdeVariableToPrintGetter(cParams)
            obj.init(cParams)            
        end
        
        function v = compute(obj)
            p = obj.physicalProblem;
            v.stress = p.variables.stress;
            v.strain = p.variables.strain;
            v.u      = obj.splitDisplacement(p.variables.d_u,p.dof.nunkn);
            v.quad   = p.element.quadrature;            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.physicalProblem = cParams.physicalProblem;
        end
        
    end
    
    methods (Access = private, Static)
        
        function uM = splitDisplacement(u,nu)
            nnode = round(length(u)/nu); 
            nodes = 1:nnode;
            uM = zeros(nnode,nu);
            for idim = 1:nu
                dofs = nu*(nodes-1)+idim;
                uM(:,idim) = u(dofs);
            end
        end       
        
    end    
    
end