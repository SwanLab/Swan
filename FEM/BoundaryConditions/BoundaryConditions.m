classdef BoundaryConditions < handle
    
    properties (Access = public)
        dirichletFun, dirichlet_dofs, dirichlet_vals, dirichlet_domain
        pointloadFun, pointload_dofs, pointload_vals
        periodic_leader, periodic_follower

        iVoigt, nVoigt
    end
    
    properties (Access = private)
        mesh
        dirichletInput
        pointloadInput
        periodicInput
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = BoundaryConditions(cParams)
            obj.init(cParams)
            obj.createDirichletFun();
            obj.createPointloadFun();
            obj.createPeriodicConditions();
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.dirichletInput = cParams.dirichletFun;
            obj.pointloadInput  = cParams.pointloadFun;
            obj.periodicInput  = cParams.periodicFun;
        end

        function createPointloadFun(obj)
            if ~isequal(obj.pointloadInput, [])
                opndofs  = obj.pointloadInput(1).fun.nDofs;
                ndimf  = obj.pointloadInput(1).fun.ndimf;
                pl = P1Function.create(obj.mesh, ndimf);
                pl_dofs = [];
                pl_vals = [];
                pl_domain = @(coor) obj.pointloadInput(1).domain(coor);
                for i = 1:numel(obj.pointloadInput)
                    values = obj.pointloadInput(i).getValues();
                    dofs   = obj.pointloadInput(i).getDofs();
        
                    % fV = dirich.fValues(:); % wrong, it needs to be overwritten
                    % fV(dofs) = values;
                    % fV = reshape(fV, [ndimf ndofs/ndimf])';
                    % dirich.fValues = fV;
                    pl_dofs = [pl_dofs; dofs];
                    pl_vals = [pl_vals; values];
                    pl_domain = @(coor) pl_domain(coor) | obj.pointloadInput(i).domain(coor);
                end
                obj.pointload_dofs = pl_dofs;
                obj.pointload_vals = pl_vals;
                obj.pointloadFun = pl;
            end
        end

        function createDirichletFun(obj)
            if ~isequal(obj.dirichletInput, [])
            opndofs  = obj.dirichletInput(1).fun.nDofs;
            ndimf  = obj.dirichletInput(1).fun.ndimf;
            dirich = P1Function.create(obj.mesh, ndimf);
            dir_fV = [];
            dir_dofs = [];
            dir_vals = [];
            dir_domain = @(coor) obj.dirichletInput(1).domain(coor);
            for i = 1:numel(obj.dirichletInput)
                values = obj.dirichletInput(i).getValues();
                dofs   = obj.dirichletInput(i).getDofs();
    
                % fV = dirich.fValues(:); % wrong, it needs to be overwritten
                % fV(dofs) = values;
                % fV = reshape(fV, [ndimf ndofs/ndimf])';
                % dirich.fValues = fV;
                dir_dofs = [dir_dofs; dofs];
                dir_vals = [dir_vals; values];
                dir_domain = @(coor) dir_domain(coor) | obj.dirichletInput(i).domain(coor);
            end
            obj.dirichlet_dofs = dir_dofs;
            obj.dirichlet_vals = dir_vals;
            obj.dirichlet_domain = dir_domain;
            obj.dirichletFun = dirich;
            end
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