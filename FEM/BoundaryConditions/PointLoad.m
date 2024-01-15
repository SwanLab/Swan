classdef PointLoad < BoundaryCondition
    
    properties (Access = public)
        fun
        domain
        type = 'Neumann';
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = PointLoad(mesh, s)
            % P1
            fun = P1Function.create(mesh, mesh.ndim); % not necessarily mesh.ndim
            pl_dofs = s.domain(mesh.coord);
            fun.fValues(pl_dofs,s.direction) = s.value;
            
            % % P2
            % fun2 = P2Function.create(mesh, mesh.ndim); % not necessarily mesh.ndim
            % pl_nods = fun2.getNodesFromCondition(s.domain);
            % fun2.fValues(pl_nods,s.direction) = s.value;
            
            obj.fun    = fun;
            obj.domain = s.domain;
            obj.mesh   = mesh;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end