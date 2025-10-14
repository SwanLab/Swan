classdef DirichletCondition < BoundaryCondition
    
    properties (Access = public)
        fun
        domain
        direction
        type = 'Dirichlet';

        dofs
        values
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = DirichletCondition(mesh, s)
            % P1
            if isfield(s,'ndim')
                ndim = s.ndim;
            else
                ndim = mesh.ndim;
            end
            fun = LagrangianFunction.create(mesh, ndim,'P1');         
            % pl_dofs = find(s.domain(mesh.coord));
            % eval_values = s.value(mesh.coord);
            % fun.fValues(pl_dofs,s.direction) = eval_values(pl_dofs);
            %pl_dofs = s.domain(mesh.coord);
            [~,l2g] = mesh.createSingleBoundaryMesh();
            fun.fValues(l2g,s.direction) = s.value;

            obj.fun    = fun;
            obj.domain = s.domain;
            obj.mesh   = mesh;
            obj.direction = s.direction;
            obj.dofs = obj.getDofs();
            obj.values = obj.getValues();
        end

        function dofs = getDofs(obj)
            ndimf = obj.fun.ndimf;
            nodesLog = obj.domain(obj.mesh.coord);
            if islogical(nodesLog)
                nodes = find(nodesLog);
            else
                nodes = nodesLog;
            end            
            dofs = ndimf*(nodes - 1) + obj.direction;
            dofs = dofs(:);
        end

        function v = getValues(obj)
            dofs = obj.domain(obj.mesh.coord);
            vals = obj.fun.fValues(dofs, obj.direction);
            v = vals(:);
        end


        function Ct = computeLinearConditionsMatrix(obj)
            % dir_dofs = sort(dirich.getDofs());
            dir_dofs = obj.getDofs();
            nDofs = obj.fun.nDofs;
            nDirich = length(dir_dofs);
            Ct = full(sparse(1:nDirich, dir_dofs, 1, nDirich, nDofs));
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end

  
    
    
end