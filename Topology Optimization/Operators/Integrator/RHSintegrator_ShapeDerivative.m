classdef RHSintegrator_ShapeDerivative < handle

    properties (Access = private)
        npnod
        mesh
        globalConnec
        fType
        fNodal
        fGauss
        xGauss
        quadOrder
        quadrature
    end

    methods (Access = public)

        % Via Integrator_Simple + Integrator
        function obj = RHSintegrator_ShapeDerivative(cParams)
            obj.init(cParams);
        end

        function rhs = compute(obj)
            rhsElem = obj.computeElementalRHS();
            rhs = obj.assembleIntegrand(rhsElem);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.npnod        = cParams.npnod;
            obj.quadOrder    = cParams.quadOrder;
            obj.globalConnec = cParams.globalConnec;
            obj.quadrature   = obj.computeQuadrature();
            if isequal(cParams.fType,'Nodal')
                obj.fNodal   = cParams.fNodal;
                obj.computeFgaussFromFnodal();
            elseif isequal(cParams.fType,'Gauss')
                obj.xGauss   = cParams.xGauss;
                obj.fGauss   = cParams.fGauss;
            end
        end
        
        function q = computeQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
        end
        
        function computeFgaussFromFnodal(obj)
            obj.computeGaussPoints();
            obj.computeFgauss();
        end

        function computeGaussPoints(obj)
            q = obj.quadrature;
            xG = repmat(q.posgp,[1,1,obj.mesh.nelem]);
            obj.xGauss = xG;
        end

        function computeFgauss(obj)
            s.fValues = obj.fNodal;
            s.mesh   = obj.mesh;
            f = P1Function(s);
            fG = f.evaluate(obj.xGauss);
            fG = permute(fG,[2 3 1]);
            obj.fGauss = fG;
        end
        
        function rhsC = computeElementalRHS(obj)
            % integrateWithShapeDerivative@RHSintegrator
            fG      = obj.fGauss;
            dV      = obj.computeDvolume();
            dN      = obj.computeGrad();
            nnode   = size(dN,2);
            ndim    = size(dN,1);
            nelem   = obj.mesh.nelem;
            int = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                for idime = 1:ndim
                    for inode = 1:nnode
                    fI     = squeezeParticular(fG(idime,igaus,:),1);
                    fdV    = fI.*dV(igaus,:);
                   %dShape = squeeze(grad(idime,:,:,igaus));
                    %intI = bsxfun(@times,dShape,fdV);                   
                    dShape = squeeze(dN(idime,inode,:,igaus))';                    
                    intI = dShape.*fdV;
                    int(inode,:) = int(inode,:) + intI;
                    end
                end
            end
            rhsC = transpose(int);
        end

        function grad = computeGrad(obj)
            m.type = obj.mesh.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(obj.quadrature,int);
            grad = g.dNdx;
         %   grad = int.deriv;
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
        
        function shapes = computeShapeFunctions(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            shapes = permute(int.shape,[1 3 2]);
        end

    end

end