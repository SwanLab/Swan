classdef LHSintegrator_Laplacian < handle

    properties (Access = private)
        fun
        mesh
        quadrature
%         material
    end

    methods (Access = public)

        function obj = LHSintegrator_Laplacian(cParams)
            %             obj.init(cParams);
            %             obj.createInterpolation();
            %             obj.createGeometry();
            obj.initLaplacian(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function LHSe = computeElementalLHS(obj)
            shapes = obj.fun.computeCartesianDerivatives(obj.quadrature);
            ngaus = size(shapes,4);
            dVolu = obj.mesh.computeDvolume(obj.quadrature);
            nDofs = size(shapes,1) * size(shapes,2);
            nElem = size(shapes,3);

%             material = obj.material;
%             Cmat = material.mu;
            
            lhs = zeros(nDofs, nDofs, nElem);
            for igaus = 1:ngaus
                dNdx = shapes(:,:,:,igaus);
                dV(1,1,:) = dVolu(igaus,:)';
                Bmat = obj.computeB(dNdx);
                Bt   = permute(Bmat,[2 1 3]);
%                 BtC  = pagemtimes(Bt,Cmat);
%                 BtCB = pagemtimes(BtC, Bmat);
                BtB = pagemtimes(Bt,Bmat);
                lhs = lhs + bsxfun(@times, BtB, dV);
            end
            LHSe = lhs;
        end

    end

    methods (Access = private)

        function initLaplacian(obj, cParams)
            obj.mesh  = cParams.mesh;
            obj.fun = cParams.fun;
%             obj.material = cParams.material;
        end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC'); % ehhh
            obj.quadrature = q;
        end

        function B = computeB(obj,dNdx)
            nunkn = size(dNdx,1);
            nnode = size(dNdx,2);
            nelem = size(dNdx,3);
            B = zeros(4,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)   = dNdx(1,i,:);
                B(2,j+1,:) = dNdx(1,i,:);
                B(3,j,:)   = dNdx(2,i,:);
                B(4,j+1,:) = dNdx(2,i,:);
            end
        end

        function LHS = assembleMatrix(obj, lhs)
            s.fun    = obj.fun; % !!!
            assembler = AssemblerFun(s);
            LHS = assembler.assemble(lhs);
        end

    end

end