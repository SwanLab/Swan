classdef LHSintegratorAnisotropicStiffness < handle

    properties (Access = private)
        CAnisotropic
        Celas
        alphaDeg
        fun
        mesh
        quadrature
        quadratureOrder
    end

    methods (Access = public)

        function obj = LHSintegratorAnisotropicStiffness(cParams)
            obj.init(cParams);
            obj.initAnisotropicTensor(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            obj.assemblyCMatrix();
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            dNdx  = obj.fun.computeCartesianDerivatives(obj.quadrature);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nGaus = obj.quadrature.ngaus;
            nElem = size(dVolu,2);
            nNodE = size(dNdx,2);
            nDofE = nNodE*obj.fun.ndimf;
            lhs = zeros(nDofE,nDofE,nElem);
            Cmat   = obj.Celas;
            Bcomp = obj.createBComputer(dNdx);
            for iGaus = 1:nGaus
                Bmat = Bcomp.compute(iGaus);
                dV(1,1,:) = dVolu(iGaus,:)';
                Bt   = permute(Bmat,[2 1 3]);
                BtC  = pagemtimes(Bt,Cmat);
                BtCB = pagemtimes(BtC, Bmat);
                lhs = lhs + bsxfun(@times, BtCB, dV);
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.fun      = cParams.fun;
            obj.mesh     = cParams.mesh;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = obj.fun.order;
            end
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = quad;
        end

        function Bcomp = createBComputer(obj, dNdx)
            s.fun  = obj.fun;
            s.dNdx = dNdx;
            Bcomp = BMatrixComputerFun(s);
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = obj.fun; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs);
        end

        function initAnisotropicTensor(obj,cParams)
            CLocal = cParams.CAnisotropic;
            obj.alphaDeg = cParams.aniAlphaDeg;
            obj.CAnisotropic = obj.rotateAnisotropicMatrix(CLocal);
        end

%         function Bcomp = createBComputer(obj)
%             s.dim          = obj.field.dim;
%             s.geometry     = obj.field.geometry;
%             Bcomp = BMatrixComputer(s);
%         end

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

%         function lhs = assembleMatrixField(obj, Ae)
%             s.dim          = obj.field.dim;
%             s.globalConnec = obj.field.connec;
%             s.nnodeEl      = obj.field.dim.nnodeElem;
%             assembler = Assembler(s);
%             lhs = assembler.assembleFields(Ae, obj.field, obj.field);
%         end
    end

end