classdef BoundaryConditions < handle

    properties (GetAccess = public)
        dirichlet
        dirichlet_values
        free
        in_elem
    end

    properties (Access = private)
        dim
        dirichletInput
        pointloadInput
        globalConnec
    end
    
    methods (Access = public)
        
        function obj = BoundaryConditions(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            [dirID, dirVals]     = obj.formatInputData(obj.dirichletInput);
            obj.dirichlet{1}        = dirID;
            obj.dirichlet_values{1} = dirVals;
            obj.free{1}             = obj.computeFreeDOF();
            obj.in_elem{1}          = obj.compute_idx();
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.dim            = cParams.dim;
            obj.globalConnec   = cParams.globalConnec;
            obj.dirichletInput = cParams.bc.dirichlet;
            obj.pointloadInput = cParams.bc.pointload;
        end

        function [dofs, vals] = formatInputData(obj, data)
            inod = data(:,1);
            iunk = data(:,2);
            dofs = obj.nod2dof(inod,iunk);
            vals = data(:,3);
        end

        function idof = nod2dof(obj, inode, iunkn)
            nunkn = obj.dim.nunkn;
            idof(:,1)= nunkn*(inode - 1) + iunkn;
        end
        
        function free = computeFreeDOF(obj)
%             free = setdiff(1:obj.dim.ndof(ifield),obj.constrained{ifield});
            ndof  = obj.dim.ndof;
            cnstr = obj.dirichlet;
            free  = setdiff(1:ndof,cnstr{1});
        end

        function dof_elem = compute_idx(obj)
            connec = obj.globalConnec;
            nunkn  = obj.dim.nunkn;
            nnode  = obj.dim.nnode;
            dof_elem  = zeros(nnode*nunkn,size(connec,1));
            for inode = 1:nnode
                for iunkn = 1:nunkn
                    idof_elem = obj.nod2dof(inode,iunkn);
                    global_node = connec(:,inode);
                    idof_global = obj.nod2dof(global_node,iunkn);
                    dof_elem(idof_elem,:) = idof_global;
                end
            end
            
        end

    end
end
