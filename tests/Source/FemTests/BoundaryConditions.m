classdef BoundaryConditions < handle

    properties (GetAccess = public)
        dirichlet
        dirichlet_values
        free
        dofsInElem
        neumann
        neumann_values
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
            [neuID, neuVals]     = obj.formatInputData(obj.pointloadInput);
            obj.dirichlet{1}        = dirID;
            obj.dirichlet_values{1} = dirVals;
            obj.neumann             = neuID;
            obj.neumann_values      = neuVals;
            obj.free{1}             = obj.computeFreeDOF();
            obj.dofsInElem{1}       = obj.computeDofsInElem();
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
            vals = data(:,3);
            dofs = obj.nod2dof(inod,iunk);
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimField;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end
        
        function free = computeFreeDOF(obj)
%             free = setdiff(1:obj.dim.ndof(ifield),obj.constrained{ifield});
            ndof  = obj.dim.ndof;
            cnstr = obj.dirichlet;
            free  = setdiff(1:ndof,cnstr{1});
        end

        function dofsElem = computeDofsInElem(obj)
            connec = obj.globalConnec;
            ndimf  = obj.dim.ndimField;
            nnode  = obj.dim.nnode;
            dofsElem  = zeros(nnode*ndimf,size(connec,1));
            for inode = 1:nnode
                for iunkn = 1:ndimf
                    idofElem   = obj.nod2dof(inode,iunkn);
                    globalNode = connec(:,inode);
                    idofGlobal = obj.nod2dof(globalNode,iunkn);
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            
        end

    end
end
