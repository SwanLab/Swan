classdef PointLoad < BoundaryCondition
    
    properties (Access = public)
        fun
        domain
        direction
        type = 'Neumann';

        dofs
        values
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = PointLoad(mesh, s)
            % P1
            fun = LagrangianFunction.create(mesh, mesh.ndim,'P1'); % not necessarily mesh.ndim
            pl_dofs = s.domain(mesh.coord);
            fun.fValues(pl_dofs,s.direction) = s.value;
            
            % % P2
            % fun2 = P2Function.create(mesh, mesh.ndim); % not necessarily mesh.ndim
            % pl_nods = fun2.getNodesFromCondition(s.domain);
            % fun2.fValues(pl_nods,s.direction) = s.value;
            
            obj.fun    = fun;
            obj.domain = s.domain;
            obj.mesh   = mesh;
            obj.direction = s.direction;
            obj.dofs = obj.getDofs();
            obj.values = obj.getValues();
        end

        function dofs = getDofs(obj)
            ndimf = obj.fun.ndimf;
            nodes = find(obj.domain(obj.mesh.coord));
            dofs = ndimf*(nodes - 1) + obj.direction;
            dofs = dofs(:);
        end

        function v = getValues(obj)
            dofs = obj.domain(obj.mesh.coord);
            vals = obj.fun.fValues(dofs, obj.direction);
            v = vals(:);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end