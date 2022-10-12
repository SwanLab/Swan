classdef LHSintegrator_Stiffness < LHSintegrator

    properties (Access = private)
        geometry
        field
    end

    methods (Access = public)

        function obj = LHSintegrator_Stiffness(cParams)
            obj.mesh  = cParams.mesh;
            obj.field = cParams.field;
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrixField(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            f = obj.field;
            dvolu = obj.mesh.computeDvolume(f.quadrature);
            ngaus = size(dvolu,1);
            nelem = obj.mesh.nelem;
            ndpe  = f.dim.ndofsElem;
            lhs = zeros(ndpe,ndpe,nelem);
            Bcomp = obj.createBComputer();
            for igaus = 1:ngaus
                Bmat = Bcomp.compute(igaus);
                nvoigt = size(Bmat,1);
                for istre = 1:nvoigt
                    BmatI = Bmat(istre,:,:);
                    BmatJ = permute(BmatI,[2 1 3]);
                    dNdN = bsxfun(@times,BmatJ,BmatI);
                    dv(1,1,:) = dvolu(igaus, :);
                    inc = bsxfun(@times,dv,dNdN);
                    lhs = lhs + inc;
                end
            end
        end

    end

    methods (Access = private)

        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interpolation;
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function Bcomp = createBComputer(obj)
            s.dim          = obj.field.dim;
            s.geometry     = obj.field.geometry;
            s.globalConnec = [];
            Bcomp = BMatrixComputer(s);
        end

        function lhs = assembleMatrixField(obj, Ae)
            s.dim          = obj.field.dim;
            s.globalConnec = obj.field.connec;
            s.nnodeEl      = obj.field.dim.nnodeElem;
            assembler = Assembler(s);
            lhs = assembler.assembleFields(Ae, obj.field, obj.field);
        end

    end

end