classdef BoundaryConditions < handle

    properties (GetAccess = public)
        dirichlet
        dirichlet_values
        free
        neumann
        neumann_values
        masterSlave
        periodic_free
        periodic_constrained
    end

    properties (Access = private)
        dim
        scale
        dirichletInput
        pointloadInput

        ndofs
    end
    
    methods (Access = public)
        
        function obj = BoundaryConditions(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            [dirID, dirVals]     = obj.formatInputData(obj.dirichletInput);
            [neuID, neuVals]     = obj.formatInputData(obj.pointloadInput);
            obj.dirichlet        = dirID;
            obj.dirichlet_values = dirVals;
            obj.neumann          = neuID;
            obj.neumann_values   = neuVals;
            obj.free             = obj.computeFreeDOF();
        end

        function red = fullToReducedMatrix(obj, mat)
            switch obj.scale
                case 'MACRO'
                    red = obj.reduceMatrixDirichlet(mat);
                case 'MICRO'
                    red = obj.reduceMatrixPeriodic(mat);
            end
        end

        function red = fullToReducedVector(obj, vec)
            switch obj.scale
                case 'MACRO'
                    red = obj.reduceVectorDirichlet(vec);
                case 'MICRO'
                    red = obj.reduceVectorPeriodic(vec);
            end
        end
        
        function full = reducedToFullVector(obj, vec)
            switch obj.scale
                case 'MACRO'
                    full = obj.expandVectorDirichlet(vec);
                case 'MICRO'
                    full = obj.expandVectorPeriodic(vec);
            end
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.dim            = cParams.dim;
            obj.scale          = cParams.scale;
            obj.dirichletInput = cParams.bc.dirichlet;
            obj.pointloadInput = cParams.bc.pointload;
            if isfield(cParams, 'ndofs') && ~isempty(cParams.ndofs)
                obj.ndofs = cParams.ndofs;
            else
                obj.ndofs = obj.dim.ndofs;
            end
            obj.initPeriodicMasterSlave(cParams);
        end

        function initPeriodicMasterSlave(obj, cParams)
            switch obj.scale
                case 'MICRO'
                    if isfield(cParams.bc, 'masterSlave')
                        obj.masterSlave = cParams.bc.masterSlave;
                    end
                    MS = obj.masterSlave;
                    if isempty(MS)
                        mesh = cParams.mesh;
                        MS = obj.computeMasterSlave(mesh.coord);
                        obj.masterSlave = MS;
                    end
                    obj.periodic_free = obj.computePeriodicNodes(MS(:,1));
                    obj.periodic_constrained = obj.computePeriodicNodes(MS(:,2));
            end
        end
        
        function perDof = computePeriodicNodes(obj,perNodes)
            nunkn = obj.dim.ndimf;
            nlib = size(perNodes,1);
            perDof = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                indDof = nlib*(iunkn - 1) + [1:nlib];
                perDof(indDof,1) = obj.nod2dof(perNodes,iunkn);
            end
        end

        function free = computeFreeDOF(obj)
            ndof  = obj.ndofs;
            cnstr = [obj.periodic_constrained;obj.dirichlet];
            free  = setdiff(1:ndof,cnstr);
        end

        function [dofs, vals] = formatInputData(obj, data)
            dofs = [];
            vals = [];
            if ~isempty(data)
                inod = data(:,1);
                iunk = data(:,2);
                vals = data(:,3);
                dofs = obj.nod2dof(inod,iunk);
            end
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimf;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

        
        function Ared = reduceMatrixDirichlet(obj,A)
%             fr = obj.computeGlobalFree();
            fr = obj.free';
            Ared = A(fr,fr);
        end
        
        function b_red = reduceVectorDirichlet(obj,b)
            fr = obj.free';
            b_red = b(fr);
        end

        function Ared = reduceMatrixPeriodic(obj,A)
            MS = obj.masterSlave;
            vF = obj.free;
            vP = obj.computePeriodicNodes(MS(:,1));
            vQ = obj.computePeriodicNodes(MS(:,2));
            vI = setdiff(vF,vP);
            
            A_II = A(vI,vI);
            A_IP = A(vI,vP) + A(vI,vQ); %Grouping P and Q nodal values
            A_PI = A(vP,vI) + A(vQ,vI); % Adding P  and Q equation
            A_PP = A(vP,vP) + A(vP,vQ) + A(vQ,vP) + A(vQ,vQ); % Adding and grouping
            
            Ared = [A_II, A_IP; A_PI, A_PP];
        end
        
        function b_red = reduceVectorPeriodic(obj,b)
            vF = obj.free;
            vP = obj.periodic_free;
            vQ = obj.periodic_constrained;
            vI = setdiff(vF,vP);
            b_I = b(vI);
            b_P = b(vP) + b(vQ);
            b_red = [b_I; b_P];
        end

        function b = expandVectorDirichlet(obj,bfree)
            dir = obj.dirichlet;
            uD  = obj.dirichlet_values;
            fr  = obj.free;
            nsteps = length(bfree(1,:));
            ndof = sum(obj.ndofs);
            uD = repmat(uD,1,nsteps);
            
            b = zeros(ndof,nsteps);
            b(fr,:) = bfree;
            if ~isempty(dir)
                b(dir,:) = uD;
            end
        end

        function b = expandVectorPeriodic(obj,bfree)
            vF = obj.free;
            vP = obj.periodic_free;
            vC = obj.periodic_constrained;
            vI = setdiff(vF,vP);
            b = zeros(obj.ndofs,1);
            b(vI) = bfree(1:1:size(vI,2));
            b(vP) = bfree(size(vI,2)+1:1:size(bfree,1));
            b(vC) = b(vP);
        end

        function MS = computeMasterSlave(obj, coord)
           mR = MasterSlaveRelator(coord);
           MS = mR.getRelation();
        end

    end
end
