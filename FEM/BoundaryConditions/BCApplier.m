classdef BCApplier < handle
    
    % Goal: group dirichlet and neumann conditions
    % to allow multiple boundary conditions at the same time
    % Use: BCApplier.computeLinearConditionsMatrix()
    properties (Access = public)
        dirichlet_dofs, dirichlet_vals
        dirichletFun
        periodic_leader, periodic_follower
    end
    
    properties (Access = private)
        mesh
        dirichletInput
        periodicInput
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = BCApplier(cParams)
            obj.init(cParams)
            obj.createDirichletFun();
            obj.createPeriodicConditions();
        end
        
        function Ct = computeLinearConditionsMatrix(obj)
            dir_dofs = obj.dirichlet_dofs;
            nDofs = obj.dirichletFun.nDofs;
            nDirich = length(dir_dofs);
            Ct = full(sparse(1:nDirich, dir_dofs, 1, nDirich, nDofs));
        end

        function Ct = computeLinearPeriodicConditionsMatrix(obj)
            per_lead = obj.periodic_leader;
            per_fllw = obj.periodic_follower;
            nDofs = obj.dirichletFun.nDofs;
            nPer = length(per_lead);
            Ct = full(sparse([(1:nPer)', (1:nPer)'], [per_lead, per_fllw], [ones(size(per_lead,1),1), -ones(size(per_lead,1),1)], nPer, nDofs));
        end

        function Ct = computeSingleDirichletPeriodicCondition(obj)
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.dirichletInput = cParams.boundaryConditions.dirichletFun;
            obj.periodicInput  = cParams.boundaryConditions.periodicFun;
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
    
                % fV = dirich.fValues(:); % wrong, it needs to be overwritten
                % fV(dofs) = values;
                % fV = reshape(fV, [ndimf ndofs/ndimf])';
                % dirich.fValues = fV;
                dir_dofs = [dir_dofs; dofs];
                dir_vals = [dir_vals; values];
            end
            obj.dirichlet_dofs = dir_dofs;
            obj.dirichlet_vals = dir_vals;
            obj.dirichletFun = dirich;
        end

        function createPeriodicConditions(obj)
            obj.periodic_leader   = [];
            obj.periodic_follower = [];
            if ~isequal(obj.periodicInput, [])
                mR = MasterSlaveRelator(obj.mesh.coord);
                MS = mR.getRelation();
                obj.periodic_leader   = obj.computePeriodicNodes(MS(:,1));
                obj.periodic_follower = obj.computePeriodicNodes(MS(:,2));
            end
        end
        
        function perDof = computePeriodicNodes(obj,perNodes)
            nunkn = obj.dirichletFun.ndimf;
            nlib = size(perNodes,1);
            perDof = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                indDof = nlib*(iunkn - 1) + [1:nlib];
                perDof(indDof,1) = obj.nod2dof(obj.dirichletFun.ndimf, perNodes,iunkn);
            end
        end

        function idof = nod2dof(obj, ndimf, inode, iunkn)
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

    end
    
end