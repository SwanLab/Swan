classdef BoundaryConditionsStokes < handle

    properties (GetAccess = public)
        dirichlet
        dirichlet_values
        free
        freeFields
        masterSlave
        periodic_free
        periodic_constrained
        neumann
        neumann_values
    end

    properties (Access = private)
        mesh
        ndimf
        ndofs
        scale
        dirichletInput
        pointloadInput
        nodesEdges
    end
    
    methods (Access = public)
        
        function obj = BoundaryConditionsStokes(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            [dirID, dirVals]     = obj.formatInputData(obj.ndimf, obj.dirichletInput);
            [neuID, neuVals]     = obj.formatInputData(obj.ndimf, obj.pointloadInput);
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

        function changeBoundaryConditions(obj,neumannDOFs,neumannValues)
            obj.neumann        = neumannDOFs;
            obj.neumann_values = neumannValues;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            if isfield(cParams,'boundaryType')
                switch cParams.boundaryType
                    case 'Periodic'
                        obj.scale = 'MICRO';
                    otherwise
                        obj.scale = 'MACRO';
                end
            else
                obj.scale = cParams.scale;
            end
            obj.ndofs          = cParams.ndofs; % Stokes
            obj.ndimf          = cParams.bc{1}.ndimf; % Elastic
            obj.initPeriodicMasterSlave(cParams);
            obj.initDirichletInput(cParams);
        end

        function initDirichletInput(obj, s)
            nfields = numel(s.bc);
            dirich  = [];
            neumann    = [];
            dirichVals = [];
            neumnnVals = [];
            free       = [];
            globalNdof = 0;
            for i = 1:nfields
                bc = s.bc{i};
%                 cornerDirichlet = bc.dirichlet;
%                 edgesDirichlet = zeros(length(obj.nodesEdges),obj.ndimf);
%                 rows = repmat(obj.nodesEdges, [size(edgesDirichlet,2), 1]);
%                 cols = repmat(1:size(edgesDirichlet,2), [length(obj.nodesEdges), 1]);
%                 cols = cols(:);
%                 newValues = zeros(size(rows,1), 1);
%                 edgesDirichlet = [rows, cols, newValues];
%                 totalDirichlet = [cornerDirichlet; edgesDirichlet];
%                 totalDirichlet = sortrows(totalDirichlet, 1);
                
                totalDirichlet = bc.dirichlet;

                obj.dirichletInput = totalDirichlet;
                obj.pointloadInput = bc.pointload;

%                 inD = bc.dirichlet;
                inD = totalDirichlet;
                inN = bc.pointload;
                [idxD, valD] = obj.formatInputData(bc.ndimf,inD);
                [idxN, valN] = obj.formatInputData(bc.ndimf,inN);
                idxD = idxD + globalNdof;
                idxN = idxN + globalNdof;

                dirich     = [dirich; idxD];
                dirichVals = [dirichVals; valD];
                neumann    = [neumann; idxN];
                neumnnVals = [neumnnVals; valN];

                firstDof = globalNdof + 1;
                lastDof  = firstDof + bc.ndofs - 1;
                obj.freeFields{i} = setdiff(firstDof:lastDof,idxD);
                globalNdof = globalNdof+ bc.ndofs;
            end
            obj.dirichlet        = dirich;
            obj.dirichlet_values = dirichVals;
            obj.neumann          = neumann;
            obj.neumann_values   = neumnnVals;
            obj.free             = obj.computeFreeDOF();
        end

        function initPeriodicMasterSlave(obj, cParams)
            switch obj.scale
                case 'MICRO'
                    if isfield(cParams.bc, 'masterSlave')
                        obj.masterSlave = cParams.bc.masterSlave;
                    end
                    MS = obj.masterSlave;
                    if isempty(MS)
                        obj.mesh = cParams.mesh;
                        MS = obj.computeMasterSlave(obj.mesh.coord);
                        obj.masterSlave = MS;
                    end
                    obj.periodic_free = obj.computePeriodicNodes(MS(:,1));
                    obj.periodic_constrained = obj.computePeriodicNodes(MS(:,2));
            end
        end
        
        function perDof = computePeriodicNodes(obj,perNodes)
            nunkn = obj.ndimf;
            nlib = size(perNodes,1);
            perDof = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                indDof = nlib*(iunkn - 1) + [1:nlib];
                perDof(indDof,1) = obj.nod2dof(obj.ndimf, perNodes,iunkn);
            end
        end

        function free = computeFieldFree(obj, ndofs, dirich)
            free  = setdiff(1:ndofs,dirich);
        end

        function free = computeFreeDOF(obj)
            ndof  = obj.ndofs;
            cnstr = [obj.periodic_constrained;obj.dirichlet];
            free  = setdiff(1:ndof,cnstr);
        end

        function [dofs, vals] = formatInputData(obj, ndimf, data)
            dofs = [];
            vals = [];
            if ~isempty(data)
                inod = data(:,1);
                iunk = data(:,2);
                vals = data(:,3);
                dofs = obj.nod2dof(ndimf, inod,iunk);
            end
        end

        function idof = nod2dof(obj, ndimf, inode, iunkn)
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
            nsteps = size(bfree,2);
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
