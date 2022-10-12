classdef LHSintegratorAnisotropicStiffness < LHSintegrator
    
    properties (Access = private)
        CAnisotropic
        Celas
        alphaDeg
        field
    end

    methods (Access = public)
        
        function obj = LHSintegratorAnisotropicStiffness(cParams)
            obj.mesh  = cParams.mesh;
            obj.field = cParams.field;
            obj.initAnisotropicTensor(cParams);
        end

        function LHS = compute(obj)
            obj.assemblyCMatrix();
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
            Cmat   = obj.Celas;
            Bcomp = obj.createBComputer();
            for igaus = 1:ngaus
                Bmat = Bcomp.compute(igaus);
                dV(1,1,:) = dvolu(igaus,:)';
                Bt   = permute(Bmat,[2 1 3]);
                BtC  = pagemtimes(Bt,Cmat);
                BtCB = pagemtimes(BtC, Bmat);
                lhs = lhs + bsxfun(@times, BtCB, dV);
            end
        end
        
   end
    
   methods (Access = private)
       
       function initAnisotropicTensor(obj,cParams)
           CLocal = cParams.CAnisotropic;
           obj.alphaDeg = cParams.aniAlphaDeg;
           obj.CAnisotropic = obj.rotateAnisotropicMatrix(CLocal);
       end

        function Bcomp = createBComputer(obj)
            s.dim          = obj.field.dim;
            s.geometry     = obj.field.geometry;
            s.globalConnec = [];
            Bcomp = BMatrixComputer(s);
        end

        function CGlobal = rotateAnisotropicMatrix(obj,CLocal)
            R = [cosd(obj.alphaDeg),-sind(obj.alphaDeg)
                sind(obj.alphaDeg), cosd(obj.alphaDeg)];
            CGlobal = R*CLocal*R';
        end

        function assemblyCMatrix(obj)
            nelem = size(obj.mesh.connec,1);
            C = repmat(obj.CAnisotropic, [1 1 nelem]);
            obj.Celas = C;
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