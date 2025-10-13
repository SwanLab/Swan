classdef BoundaryConditions < handle
    
    properties (Access = public)
        dirichletFun, dirichlet_dofs, dirichlet_vals, dirichlet_domain
        tractionFun
        periodic_leader, periodic_follower

        iVoigt, nVoigt
        free_dofs
    end
    
    properties (Access = private)
        mesh
        dirichletInput
        periodicInput
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = BoundaryConditions(cParams)
            obj.init(cParams)
            obj.createDirichletFun();
            obj.createPeriodicConditions();
        end

        function updatePeriodicConditions(obj,MS)
            obj.periodic_leader   = obj.computePeriodicNodes(MS(:,1));
            obj.periodic_follower = obj.computePeriodicNodes(MS(:,2));
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.dirichletInput = cParams.dirichletFun;
            obj.tractionFun = cParams.pointloadFun;
            obj.periodicInput  = cParams.periodicFun;
        end

%         function createPointloadFun(obj)
%             [dofs,vals,domain,fun] = obj.createBCFun(obj.pointloadInput);
%             obj.traction_dofs = dofs;
%             obj.traction_vals = vals;
%             obj.traction_domain = domain;
%             obj.tractionFun = fun;
%         end

        function createDirichletFun(obj)
            [dofs,vals,domain,fun] = obj.createBCFun(obj.dirichletInput);
            obj.dirichlet_dofs = dofs;
            obj.dirichlet_vals = vals;
            obj.dirichlet_domain = domain;
            obj.dirichletFun = fun;
            if ~isempty(fun)
                obj.free_dofs = setdiff(1:fun.nDofs,dofs);
            else
                obj.free_dofs = [];
            end
        end

        function [dofs,vals,domain,bcFun] = createBCFun(obj,input)
            if ~isequal(input, [])
                ndimf  = input(1).fun.ndimf;
                bcFun = LagrangianFunction.create(obj.mesh, ndimf,'P1');
                dofs = [];
                vals = [];
                domain = @(coor) input(1).domain(coor);
                for i = 1:numel(input)
                    values_i = input(i).getValues();
                    dofs_i   = input(i).getDofs();

                    dofs = [dofs; dofs_i];
                    vals = [vals; values_i];
                    domain = @(coor) pl_domain(coor) | input(i).domain(coor);
                    bcFun.setFValues(obj.addValues(bcFun,dofs_i,values_i));
                end
            else
                dofs = [];
                vals = [];
                domain = [];
                bcFun = [];
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
        
        function fV = addValues(obj,dirich,dofs,values)
            ndofs = dirich.nDofs;
            ndimf = dirich.ndimf;

            fVdofs = zeros(ndofs,1);
            fVdofs(dofs) = values;
            fVdofs = reshape(fVdofs,[ndimf ndofs/ndimf])';
            fV = dirich.fValues + fVdofs;
        end
    end
    
end