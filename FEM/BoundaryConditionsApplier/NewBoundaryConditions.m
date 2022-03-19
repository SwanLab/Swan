classdef NewBoundaryConditions < handle

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
        dofsInElem
        applierType
        dirichletInput
        pointloadInput
    end
    
    methods (Access = public)
        
        function obj = NewBoundaryConditions(cParams)
            obj.init(cParams);
        end

%         function updateDirichlet(obj)
%         end

        function compute(obj)
            [dirID, dirVals]     = obj.formatInputData(obj.dirichletInput);
            [neuID, neuVals]     = obj.formatInputData(obj.pointloadInput);
            obj.dirichlet{1}        = dirID;
            obj.dirichlet_values{1} = dirVals;
            obj.neumann             = neuID;
            obj.neumann_values      = neuVals;
            obj.free{1}             = obj.computeFreeDOF();
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
            obj.applierType    = cParams.type;
            obj.dofsInElem     = cParams.dofsInElem;
            obj.dirichletInput = cParams.bc.dirichlet;
            obj.pointloadInput = cParams.bc.pointload;
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
                       mesh.computeMasterSlaveNodes();
                       MS = mesh.masterSlaveNodes;
                    end
                    obj.periodic_free = obj.computePeriodicNodes(MS(:,1));
                    obj.periodic_constrained = obj.computePeriodicNodes(MS(:,2));
            end
        end
        
        function periodic_dof = computePeriodicNodes(obj,periodic_nodes)
            nunkn = obj.dim.ndimField;
            nlib = size(periodic_nodes,1);
            periodic_dof = zeros(nlib*nunkn,1);
            for iunkn = 1:nunkn
                index_glib = nlib*(iunkn - 1) + [1:nlib];
                periodic_dof(index_glib,1) = obj.nod2dof(periodic_nodes,iunkn);
            end
        end

        function free = computeFreeDOF(obj)
            ndof  = obj.dim.ndof;
            cnstr = [obj.periodic_constrained;obj.dirichlet{1}];
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
            ndimf = obj.dim.ndimField;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

        
        function Ared = reduceMatrixDirichlet(obj,A)
%             fr = obj.computeGlobalFree();
            fr = obj.free{1}';
            Ared = A(fr,fr);
        end
        
        function b_red = reduceVectorDirichlet(obj,b)
            fr = obj.free{1}';
            b_red = b(fr);
        end

        function Ared = reduceMatrixPeriodic(obj,A)
            MS = obj.computeMasterSlave();
            vF = obj.free{1};
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
            vF = obj.free{1};
            vP = obj.periodic_free;
            vQ = obj.periodic_constrained;
            vI = setdiff(vF,vP);
            b_I = b(vI);
            b_P = b(vP) + b(vQ);
            b_red = [b_I; b_P];
        end

        function b = expandVectorDirichlet(obj,bfree)
            dir = obj.dirichlet{1};
            uD = obj.dirichlet_values{1};
            fr = obj.free{1};
            nsteps = length(bfree(1,:));
            ndof = sum(obj.dim.ndof);
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
            vI = setdiff(vF{1},vP);
            
            b = zeros(obj.dim.ndof,1);
            b(vI) = bfree(1:1:size(vI,2));
            b(obj.periodic_free) = bfree(size(vI,2)+1:1:size(bfree,1));
            b(obj.periodic_constrained) = b(obj.periodic_free);
        end

        function MS = computeMasterSlave(obj)
            MS = obj.masterSlave;
            if isempty(MS)
               mesh = cParams.mesh;
               mesh.computeMasterSlaveNodes();
               MS = mesh.masterSlaveNodes;
            end
        end

    end
end
