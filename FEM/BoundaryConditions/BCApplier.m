classdef BCApplier < handle
    
    % Goal: group dirichlet and neumann conditions
    % to allow multiple boundary conditions at the same time
    % Use: BCApplier.computeLinearConditionsMatrix()
    properties (Access = public)
        dirichlet_dofs, dirichlet_vals
    end
    
    properties (Access = private)
        mesh
        dirichletInput, dirichletFun
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = BCApplier(cParams)
            obj.init(cParams)
            obj.createDirichletFun();
        end
        
        function Ct = computeLinearConditionsMatrix(obj)
            dir_dofs = obj.dirichlet_dofs;
            nDofs = obj.dirichletFun.nDofs;
            nDirich = length(dir_dofs);
            Ct = full(sparse(1:nDirich, dir_dofs, 1, nDirich, nDofs));
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.dirichletInput = cParams.boundaryConditions.dirichletFun;
        end

        function createDirichletFun(obj)
            ndofs  = obj.dirichletInput(1).fun.nDofs;
            ndimf  = obj.dirichletInput(1).fun.ndimf;
            dirich = P1Function.create(obj.mesh, ndimf);
            dir_fV = [];
            dir_dofs = [];
            dir_vals = [];
            for i = 1:numel(obj.dirichletInput)
                values = obj.dirichletInput(i).getValues();
                dofs   = obj.dirichletInput(i).getDofs();
    
                fV = dirich.fValues(:); 
                fV(dofs) = values;
                fV = reshape(fV, [ndimf ndofs/ndimf])';
                dir_dofs = [dir_dofs; dofs];
                dir_vals = [dir_vals; values];
            end
            obj.dirichlet_dofs = dir_dofs;
            obj.dirichlet_vals = dir_vals;
            dirich.fValues = fV;
            obj.dirichletFun = dirich;
        end
        
    end
    
end