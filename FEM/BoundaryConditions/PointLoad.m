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
        
        function obj = PointLoad(mesh, domain, direction, value)
            fun = P1Function.create(mesh, mesh.ndim); % not necessarily mesh.ndim
            pl_dofs = domain(mesh.coord);
            fun.fValues(pl_dofs,direction) = value;
            
            obj.fun    = fun;
            obj.domain = domain;
            obj.mesh   = mesh;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end