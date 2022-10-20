classdef LHSintegrator_StiffnessElastic < LHSintegrator

    properties (Access = private)
        geometry
        material
        field
    end

    methods (Access = public)

        function obj = LHSintegrator_StiffnessElastic(cParams)
            obj.mesh = cParams.mesh;
            obj.field = cParams.field;
            obj.material = cParams.material;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrixField(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            dvolu  = obj.mesh.computeDvolume(obj.field.quadrature);
            ngaus  = obj.field.quadrature.ngaus;
            nelem  = size(obj.material.C,3);
            npe    = obj.field.dim.ndofsElem;
            lhs = zeros(npe,npe,nelem);
            Bcomp = obj.createBComputer();
            Cmat = obj.material.C(:,:,:,1);
            for igaus = 1:ngaus
                Bmat = Bcomp.compute(igaus);
%                 Cmat = obj.material.C(:,:,:,igaus);
                dV(1,1,:) = dvolu(igaus,:)';
                Bt   = permute(Bmat,[2 1 3]);
                BtC  = pagemtimes(Bt,Cmat);
                BtCB = pagemtimes(BtC, Bmat);
                lhs = lhs + bsxfun(@times, BtCB, dV);
            end
        end

    end

    methods (Access = private)

        function Bcomp = createBComputer(obj)
            s.dim          = obj.field.dim;
            s.geometry     = obj.field.geometry;
            Bcomp = BMatrixComputer(s);
        end

        function LHS = assembleMatrixField(obj, lhs)
            s.dim          = obj.field.dim;
            s.globalConnec = obj.field.connec;
            s.nnodeEl      = obj.field.dim.nnodeElem;
            assembler = Assembler(s);
            LHS = assembler.assembleFields(lhs, obj.field, obj.field);
        end
    end

end