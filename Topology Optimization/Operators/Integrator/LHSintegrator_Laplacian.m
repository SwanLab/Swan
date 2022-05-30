classdef LHSintegrator_Laplacian < LHSintegrator

    properties (Access = private)
        field
        material
    end

    methods (Access = public)

        function obj = LHSintegrator_Laplacian(cParams)
            %             obj.init(cParams);
            %             obj.createQuadrature();
            %             obj.createInterpolation();
            %             obj.createGeometry();
            obj.initLaplacian(cParams);
        end

        function lhs = compute(obj)
            lhs = obj.computeElementalLHS();
            %             LHS = obj.assembleMatrix(lhs);
        end

    end

    methods (Access = protected)

        function LHSe = computeElementalLHS(obj)
            f = obj.field;
            ndofs = f.dim.ndofsElem;
            nelem = obj.mesh.nelem;

%             material = obj.material;
%             Cmat = material.mu;
            
            geom = f.geometry;
            shape = geom.dNdx;
            ngaus = size(shape,4);
            dvolu = geom.dvolu;
            lhs = zeros(ndofs, ndofs, nelem);
            for igaus = 1:ngaus
                dNdx = shape(:,:,:,igaus);
                dV(1,1,:) = dvolu(:,igaus);
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
            obj.mesh     = cParams.mesh;
            obj.field    = cParams.field;
%             obj.material = cParams.material;
        end

        function B = computeB(obj,dNdx)
            f = obj.field;
            nunkn = f.dim.ndimf;
            nnode = f.dim.nnodeElem;
            nelem = obj.mesh.nelem;
            B = zeros(4,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)  = dNdx(1,i,:);
                B(2,j+1,:)= dNdx(1,i,:);
                B(3,j,:)  = dNdx(2,i,:);
                B(4,j+1,:)= dNdx(2,i,:);
            end
        end

    end

end