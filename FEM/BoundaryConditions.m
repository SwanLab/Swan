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
        dirichletInput
        pointloadInput

    end
    
    methods (Access = public)
        
        function obj = BoundaryConditions(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            [dirID, dirVals]     = obj.formatInputData(obj.dirichletInput);
            [neuID, neuVals]     = obj.formatInputData(obj.pointloadInput);
            obj.dirichlet{1}        = dirID;
            obj.dirichlet_values{1} = dirVals;
            obj.neumann             = neuID;
            obj.neumann_values      = neuVals;
            obj.free{1}             = obj.computeFreeDOF();
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
%             free = setdiff(1:obj.dim.ndof(ifield),obj.constrained{ifield});
            ndof  = obj.dim.ndof;
            cnstr = [obj.periodic_constrained;obj.dirichlet{1}];
            free  = setdiff(1:ndof,cnstr);
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.dim            = cParams.dim;
            obj.dirichletInput = cParams.bc.dirichlet;
            obj.pointloadInput = cParams.bc.pointload;
            if isfield(cParams.bc, 'masterSlave')
                obj.masterSlave = cParams.bc.masterSlave;
            else
%                obj.mesh.computeMasterSlaveNodes();
%                obj.masterSlave = obj.mesh.masterSlaveNodes;
            end
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

    end
end
