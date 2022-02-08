classdef ForcesComputer < handle

    properties (Access = private)
        dim
        mesh
        dof
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

    end

    methods (Access = private)

        function init(obj, s)
            obj.dim  = s.dim;
            obj.mesh = s.mesh;
            obj.dof  = s.dof;
        end

        function Fs = computeSuperficialFext(obj)
            d = obj.dim;
            nnode = d.nnode;
            nunkn = d.nunkn;
            nelem = d.nelem;
            Fs = zeros(nnode*nunkn,1,nelem);
        end
        
        function Fv = computeVolumetricFext(obj)
            d = obj.dim;
            nnode = d.nnode;
            nunkn = d.nunkn;
            nelem = d.nelem;
            Fv = zeros(nnode*nunkn,1,nelem);
        end

        function b = assembleVector(obj, Fsup, Fvol)
            bElem_cell = {Fsup + Fvol};
            nfields = 1;
            d = obj.dim;
            for iField = 1:nfields
                bElem = bElem_cell{iField,1};
                b = zeros(obj.dim.ndof(iField),1);
                nDofPerElem = d.ndofPerElement;
                nDof = d.ndof(iField);
                nGaus = size(bElem,2);
                for iDof = 1:nDofPerElem
                    for igaus = 1:nGaus
                        c = squeeze(bElem(iDof,igaus,:));
                        idof_elem = obj.dof.in_elem{iField}(iDof,:);
                        b = b + sparse(idof_elem,1,c',nDof,1);
                    end
                end
                b_global{iField,1} = b;
            end
            b=cell2mat(b_global);
        end

        function Fp = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            Fp = zeros(obj.dim.ndof,1);
            if ~isempty(obj.dof.neumann)
                Fp(obj.dof.neumann) = obj.dof.neumann_values;
            end
        end

    end
    
end