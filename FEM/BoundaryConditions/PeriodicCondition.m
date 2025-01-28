classdef PeriodicCondition < BoundaryCondition
    
    properties (Access = public)
        fun
        domain
        direction
        type = 'Periodic';

        dofs
        values
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = PeriodicCondition(mesh, s)
            % P1
            % fun = P1Function.create(mesh, mesh.ndim); % not necessarily mesh.ndim
            % pl_dofs = s.domain(mesh.coord);
            % fun.fValues(pl_dofs,s.direction) = s.value;
            % 
            % obj.fun    = fun;
            % obj.domain = s.domain;
            % obj.mesh   = mesh;
            % obj.direction = s.direction;
            % obj.dofs = obj.getDofs();
            % obj.values = obj.getValues();
        end

        function dofs = getDofs(obj)
            % ndimf = obj.fun.ndimf;
            % nodes = find(obj.domain(obj.mesh.coord));
            % dofs = ndimf*(nodes - 1) + obj.direction;
            % dofs = dofs(:);
        end

        function v = getValues(obj)
            % dofs = obj.domain(obj.mesh.coord);
            % vals = obj.fun.fValues(dofs, obj.direction);
            % v = vals(:);
        end

        function Ct = computeLinearConditionsMatrix(obj)
            % % dir_dofs = sort(dirich.getDofs());
            % dir_dofs = obj.getDofs();
            % nDofs = obj.fun.nDofs;
            % nDirich = length(dir_dofs);
            % Ct = full(sparse(1:nDirich, dir_dofs, 1, nDirich, nDofs));
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end