classdef ForcesComputer < handle

    properties (Access = private)
        dim
        mesh
        boundaryConditions
    end
    
    methods (Access = public)

        function obj = ForcesComputer(cParams)
            obj.init(cParams);
        end

        function Fext = compute(obj)
            Fsup       = obj.computeSuperficialFext;
            Fvol       = obj.computeVolumetricFext;
            FextSupVol = obj.assembleVector(Fsup, Fvol);
            Fpoint     = obj.computePunctualFext();
            Fext = FextSupVol +  Fpoint;
        end

        function R = computeReactions(obj, K)
            bc      = obj.boundaryConditions;
            dirich  = bc.dirichlet{1};
            dirichV = bc.dirichlet_values{1};
            if ~isempty(dirich)
                R = -K(:,dirich)*dirichV;
            else
                R = zeros(sum(obj.dim.ndof(:)),1);
            end

        end

    end

    methods (Access = private)

        function init(obj, s)
            obj.dim                = s.dim;
            obj.mesh               = s.mesh;
            obj.boundaryConditions = s.BC;
        end

        function Fs = computeSuperficialFext(obj)
            d = obj.dim;
            nnode = d.nnode;
            ndimf = d.ndimField;
            nelem = d.nelem;
            Fs = zeros(nnode*ndimf,1,nelem);
        end
        
        function Fv = computeVolumetricFext(obj)
            d = obj.dim;
            nnode = d.nnode;
            ndimf = d.ndimField;
            nelem = d.nelem;
            Fv = zeros(nnode*ndimf,1,nelem);
        end

        function b = assembleVector(obj, Fsup, Fvol)
            forces = squeeze(Fsup + Fvol);
            ndofPerElem = obj.dim.ndofPerElement;
            ndof        = obj.dim.ndof;
            dofsInElemCell = obj.boundaryConditions.dofsInElem;
            dofsInElem = cell2mat(dofsInElemCell);
            b = zeros(ndof,1);
            for iDof = 1:ndofPerElem
                dofs = dofsInElem(iDof,:);
                c = forces(iDof,:);
                badd = obj.computeAddVectorBySparse(dofs, c);
                % badd = obj.computeAddVectorByAccumarray(dofs, c);
                b = b + badd;
            end
        end

        function Bdof = computeAddVectorBySparse(obj,dofs, c)
           ndof = obj.dim.ndof;
           Bdof = sparse(dofs,1,c',ndof,1);
        end

        function Bdof = computeAddVectorByAccumarray(obj,dofs,c)
           ndof = obj.dim.ndof;
           Bdof = accumarray(dofs',c',[ndof 1]);
        end

        function Fp = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            neumann       = obj.boundaryConditions.neumann;
            neumannValues = obj.boundaryConditions.neumann_values;
            Fp = zeros(obj.dim.ndof,1);
            if ~isempty(neumann)
                Fp(neumann) = neumannValues;
            end
        end

    end
    
end