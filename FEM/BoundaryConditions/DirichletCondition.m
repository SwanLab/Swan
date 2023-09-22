classdef DirichletCondition < BoundaryCondition
    
    properties (Access = public)
        fun
        domain
        type = 'Dirichlet';
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = DirichletCondition(mesh, domain, direction, value)
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