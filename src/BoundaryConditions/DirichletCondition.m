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
        
        function obj = DirichletCondition(mesh, s, varargin)
            % P1
            if isfield(s,'ndim')
                ndim = s.ndim;
            else
                ndim = mesh.ndim;
            end

            if ~isempty(varargin)
                order = varargin{1};
            else
                order = 'P1';
            end
            fun = LagrangianFunction.create(mesh, ndim,order);
            % pl_dofs = find(s.domain(mesh.coord));
            % eval_values = s.value(mesh.coord);
            % fun.fValues(pl_dofs,s.direction) = eval_values(pl_dofs);
            dofs = fun.getDofsFromCondition(s.domain);
            nodes = unique(ceil(dofs/fun.ndimf));
            fun.fValues(nodes,s.direction) = s.value;
            
%             pl_dofs = s.domain(mesh.coord);
%             fun.fValues(pl_dofs,s.direction) = s.value;

            obj.fun    = fun;
            obj.domain = s.domain;
            obj.mesh   = mesh;
            obj.direction = s.direction;
%             obj.dofs = obj.getDofs();
%             obj.values = obj.getValues();

                  obj.dofs = dofs;
            obj.values = fun.fValues(nodes,s.direction);
        end

        function dofs = getDofs(obj)
            dofs = obj.dofs;
%             ndimf = obj.fun.ndimf;
%             nodesLog = obj.domain(obj.mesh.coord);
%             if islogical(nodesLog)
%                 nodes = find(nodesLog);
%             else
%                 nodes = nodesLog;
%             end            
%             dofs = ndimf*(nodes - 1) + obj.direction;
%             dofs = dofs(:);
        end

        function v = getValues(obj)
              dofs  = obj.dofs;
              nodes = unique(ceil(dofs/obj.fun.ndimf));
              vals = obj.fun.fValues(nodes, obj.direction);
              v = vals(:);
%             dofs = obj.domain(obj.mesh.coord);
%             vals = obj.fun.fValues(dofs, obj.direction);
%             v = vals(:);
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