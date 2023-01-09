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
            if isempty(p.variables)
                p.computeChomog();
            end
            v.stress = p.variables.stress;
            v.strain = p.variables.strain;
%             v.u      = obj.splitDisplacement(p.variables.d_u,p.getDimensions().ndimf);
            v.u      = p.variables.d_u;
            v.quad   = p.getQuadrature();
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