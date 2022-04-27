classdef RHSintegrator_ShapeDerivative < handle

    properties (Access = private)
        npnod
        mesh
        globalConnec
        fNodal
        xGauss
        fGauss
        quadOrder
        quadrature
    end

    methods (Access = public)

        % Via Integrator_Simple + Integrator
        function obj = RHSintegrator_ShapeDerivative(cParams)
            obj.init(cParams);
        end

        function rhs = compute(obj)
            obj.computeQuadrature();
            obj.computeGaussPoints();
            obj.computeFgauss();
            rhs = obj.integrateFgauss();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.npnod        = cParams.npnod;
            obj.fNodal       = cParams.fNodal;
            obj.quadOrder    = cParams.quadOrder;
            obj.globalConnec = cParams.globalConnec;
        end
        
        function computeQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end
        
        function computeGaussPoints(obj)
            q = obj.quadrature;
            xG = repmat(q.posgp,[1,1,obj.mesh.nelem]);
            obj.xGauss = xG;
        end

        function computeFgauss(obj)
            s.fNodes = obj.fNodal;
            s.connec = obj.globalConnec;
            s.type   = obj.mesh.type;
            f = FeFunction(s);
            fG = f.interpolateFunction(obj.xGauss);
            fG = permute(fG,[2 3 1]);
            obj.fGauss = fG;
        end
        
        function rhs = integrateFgauss(obj)
            rhsElem = obj.computeElementalRHS();
            rhs = obj.assembleIntegrand(rhsElem);
        end
        
        function rhs = computeElementalRHS(obj)
            % integrateWithShapeDerivative@RHSintegrator
            fG      = obj.fGauss;
            dV      = obj.computeDvolume();
            grad    = obj.computeGrad();
            nnode   = size(grad,2);
            ndim    = size(grad,1);
            nelem   = obj.mesh.nelem;
            int = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                for idime = 1:ndim
                    fI     = squeezeParticular(fG(idime,igaus,:),1);
                    fdV    = (fI.*dV(igaus));
                    dShape = squeeze(grad(idime,:,:,igaus));
                    intI = bsxfun(@times,dShape,fdV);
                    int = int + intI;
                end
            end
            rhs = transpose(int);
        end

        function f = assembleIntegrand(obj,rhsElem)
            integrand = rhsElem;
            ndofs = obj.npnod;
            connec = obj.globalConnec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end
        
        function dV = computeDvolume(obj)
            q = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);
        end

        function grad = computeGrad(obj)
            m.type = obj.mesh.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(obj.quadrature,int);
            grad = g.dNdx;
        end

    end

end